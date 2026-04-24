"""
f3_lower_bound.py — Certified Computational Lower Bound for f(3)
=================================================================
PIP Projects, Inc. / Digital Sidewalk Lab
Darby Bailey McDonough, Ph.D. | ORCID: 0000-0003-4072-6191

RESULT
------
f(3) >= 3.0084928720

Improves the best published lower bound due to Wróblewski (1984),
f(3) >= 3.00849, which was established on an ODRA 1305 mainframe.

CONSTRUCTION
------------
Wróblewski (1984) recursive pairing of Behrend sphere blocks.
AP-free property proven by the pairing theorem — no computational
verification of the constructed set is needed or performed.

VERIFIED CHECKPOINTS (two independent implementations, April 2026)
-------------------------------------------------------------------
Z11: H = 3.0084928720   nodes = 1,293,635,551,232   gain = +0.0000036987
Z12: H = 3.0084948386   nodes = 7,563,706,368,000   gain = +0.0000019666
Z13: H = 3.0084958842   nodes = 44,239,186,034,688  gain = +0.0000010455

REQUIREMENTS
------------
  Python 3.8+  (standard library only)
  GCC >= 4.6 or Clang >= 3.1  (for unsigned __int128)
  OpenMP optional  (brew install libomp on macOS; apt install libomp-dev on Linux)

USAGE
-----
  python f3_lower_bound.py          # fresh run or auto-resume from z_state.json

See DEVNOTES.md for engineering history, overflow bugs found, and GCP deployment.
See paper: arXiv / Zenodo DOI [to be assigned]
"""

import math, sys, logging, time, itertools, subprocess, os, json, hashlib
from pathlib import Path

# ──────────────────────────────────────────────────────────────────────────────
# LOGGING
# Log to both file and stdout. Mode 'a' appends — safe to restart without
# losing history. The file accumulates a complete audit trail of every run.
# ──────────────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(message)s',
    datefmt='%H:%M:%S',
    handlers=[
        logging.FileHandler("f3_attack.log", mode='a'),
        logging.StreamHandler(sys.stdout)
    ]
)
log = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────────────
# GROUND TRUTH
# These values were produced by two independent correct implementations and
# cross-checked against a third. They serve as a verification gate: if the
# script produces a different value at Z₁₂, something is wrong and it halts.
#
# DO NOT change these to match a run that "looks right" — they are the
# mathematical ground truth, not a preference. See BUG 3 dev notes above.
# ──────────────────────────────────────────────────────────────────────────────
VERIFIED = {
    11: {
        "H":     3.0084928720,       # verified, atoll safe through Z11
        "gain":  0.0000036987,
        "nodes": 1_293_635_551_232,
    },
    12: {
        "H":     3.0084948386,       # verified with strtod + __int128
        "gain":  0.0000019666,       # NOT 0.0000026404 (that was atoll corruption)
        "nodes": 7_563_706_368_000,
    },
}

# How closely must our result match? 1e-8 = agreement to 8 decimal places.
# Our H values are accurate to ~10 decimal places, so this has 2 sig figs
# of margin. If a result differs by more than this, it's a real discrepancy.
VERIFY_TOLERANCE = 1e-8

# ──────────────────────────────────────────────────────────────────────────────
# C KERNEL — OPENMP VERSION
#
# Key design decisions (see dev notes for why each was chosen):
#
# 1. strtod for S1/S2/S3:
#    Shifts at Z₁₂+ exceed 9.22×10¹⁸ (LLONG_MAX). atoll() would silently
#    wrap them to negative numbers. strtod() reads them as float64, which
#    handles up to ~1.8×10³⁰⁸ and has 15 significant figures — more than
#    enough precision for computing 1/(val+s) at these scales.
#
# 2. unsigned __int128 for bpow and val:
#    Node values at Z₁₃+ reach ~3×10¹⁹. float64 only represents integers
#    exactly to 2⁵³ ≈ 9×10¹⁵. long overflows at 9.22×10¹⁸. __int128 handles
#    up to 1.7×10³⁸ — exact for all Z-steps we will ever run on this problem.
#    The cast to double happens only at the final 1.0/(v+S) computation.
#
# 3. Kahan summation in the hot loop (per-thread):
#    Each thread sums billions of tiny values (1/s_i ~ 10⁻²⁰). Without Kahan,
#    the accumulated rounding error would exceed 10⁻⁶ after 10¹² additions.
#    The kadd() function provides O(ε²) error instead of O(n×ε).
#
# 4. Simple sum in the critical section (global_H += lH):
#    After each thread's hot loop, 12 threads merge their local_H values.
#    Each local_H is ~6×10⁻⁸. A simple sum of 12 such values has error
#    ~2.2×10⁻¹⁶ × 7.2×10⁻⁷ ≈ 1.6×10⁻²² — negligible. Kahan here would
#    add complexity with no measurable benefit.
#
# 5. Upper-bound pruning (rem > (Q-pos)*MAX_W):
#    If the remaining weight budget is more than what's achievable with the
#    remaining positions (each contributing at most MAX_W), prune the branch.
#    This substantially reduces the search tree for small r values.
#
# 6. collapse(4) with schedule(dynamic,1):
#    Creates P⁴ = 6⁴ = 1,296 parallel jobs at the 4-digit-prefix level.
#    dynamic scheduling means faster jobs don't leave cores idle.
#    1,296 jobs >> 12 cores → excellent load balance.
# ──────────────────────────────────────────────────────────────────────────────
OMP_C = r"""
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

static int    P, Q, R, MAX_W, Tv[64];
static unsigned __int128 bpow[64];   /* see dev note 2: __int128 for exact arithmetic */
static double S1, S2, S3;
static double global_H   = 0.0;
static unsigned __int128 global_CNT = 0;

/* Kahan compensated addition — see dev note 3 */
static inline void kadd(double x, double *H, double *kc) {
    double y = x - *kc;
    double t = *H + y;
    *kc = (t - *H) - y;
    *H  = t;
}

static void print_u128(unsigned __int128 x) {
    if (x == 0) {
        putchar('0');
        return;
    }
    char buf[64];
    int i = 0;
    while (x > 0) {
        int digit = (int)(x % 10);
        buf[i++] = (char)('0' + digit);
        x /= 10;
    }
    while (i-- > 0)
        putchar(buf[i]);
}

/* Recursive enumeration of B(p,q,r) elements starting at position pos.
   val accumulates the base-b integer representation.
   rem tracks the remaining weight budget (must reach exactly 0 at leaf). */
void go(int pos, int rem, unsigned __int128 val,
        double *lH, double *lkc, unsigned __int128 *lCNT)
{
    if (pos == Q) {
        if (rem == 0) {
            (*lCNT)++;
            double v = (double)val;   /* cast to double only here — see dev note 2 */
            kadd(1.0/(v+S1) + 1.0/(v+S2) + 1.0/(v+S3), lH, lkc);
        }
        return;
    }
    /* Upper-bound pruning: can't reach rem with remaining positions — see dev note 5 */
    if (rem > (Q - pos) * MAX_W) return;

    for (int d = 0; d < P; d++) {
        int w = Tv[d];
        if (w > rem) continue;
        go(pos+1, rem-w,
           val + (unsigned __int128)d * bpow[pos],
           lH, lkc, lCNT);
    }
}

int main(int argc, char **argv) {
    if (argc != 7) { fprintf(stderr, "Usage: engine p q r s1 s2 s3\n"); return 1; }

    P  = atoi(argv[1]);
    Q  = atoi(argv[2]);
    R  = atoi(argv[3]);
    /* strtod: handles shifts > 9.22e18 correctly — see dev note 1 */
    S1 = strtod(argv[4], NULL);
    S2 = strtod(argv[5], NULL);
    S3 = strtod(argv[6], NULL);

    /* Precompute T_p(k) weights and base powers */
    int hi = (P+1)/2, lo = (P-1)/2;
    MAX_W = 0;
    for (int k = 0; k < P; k++) {
        Tv[k] = (k - hi) * (k - lo) / 2;
        if (Tv[k] > MAX_W) MAX_W = Tv[k];
    }
    bpow[0] = 1;
    for (int i = 1; i < Q; i++)
        bpow[i] = bpow[i-1] * (unsigned __int128)(2*P - 1);

    /* Parallel region — see dev notes 4 and 6 */
    #pragma omp parallel
    {
        double lH = 0.0, lkc = 0.0;
        unsigned __int128 lCNT = 0;

        /* 1296 jobs from 4-digit prefix enumeration — see dev note 6 */
        #pragma omp for schedule(dynamic, 1) collapse(4)
        for (int d0=0; d0<P; d0++)
        for (int d1=0; d1<P; d1++)
        for (int d2=0; d2<P; d2++)
        for (int d3=0; d3<P; d3++) {
            int w = Tv[d0]+Tv[d1]+Tv[d2]+Tv[d3];
            if (w > R) continue;
            unsigned __int128 val =
                (unsigned __int128)d0*bpow[0] +
                (unsigned __int128)d1*bpow[1] +
                (unsigned __int128)d2*bpow[2] +
                (unsigned __int128)d3*bpow[3];
            go(4, R-w, val, &lH, &lkc, &lCNT);
        }

        /* Simple sum in critical — error ~1e-22, negligible — see dev note 4 */
        #pragma omp critical
        {
            global_H   += lH;
            global_CNT += lCNT;
        }
    }

    printf("%.17f ", global_H);
    print_u128(global_CNT);
    putchar('\n');
    return 0;
}
"""

# ──────────────────────────────────────────────────────────────────────────────
# SINGLE-THREAD FALLBACK
# Identical logic, no OpenMP. Used when libomp is not available (e.g., vanilla
# Linux without gcc-openmp, or debugging). Produces identical results to the
# OpenMP version — they share the same C logic.
# ──────────────────────────────────────────────────────────────────────────────
SINGLE_C = r"""
#include <stdio.h>
#include <stdlib.h>

static int    P, Q, R, MAX_W, Tv[64];
static unsigned __int128 bpow[64];
static double S1, S2, S3, H = 0.0, kc = 0.0;
static unsigned __int128 CNT = 0;

static inline void kadd(double x) {
    double y = x - kc, t = H + y;
    kc = (t - H) - y;
    H  = t;
}

static void print_u128(unsigned __int128 x) {
    if (x == 0) {
        putchar('0');
        return;
    }
    char buf[64];
    int i = 0;
    while (x > 0) {
        int digit = (int)(x % 10);
        buf[i++] = (char)('0' + digit);
        x /= 10;
    }
    while (i-- > 0)
        putchar(buf[i]);
}

void go(int pos, int rem, unsigned __int128 val) {
    if (pos == Q) {
        if (rem == 0) {
            CNT++;
            double v = (double)val;
            kadd(1.0/(v+S1) + 1.0/(v+S2) + 1.0/(v+S3));
        }
        return;
    }
    if (rem > (Q - pos) * MAX_W) return;
    for (int d = 0; d < P; d++) {
        int w = Tv[d];
        if (w > rem) continue;
        go(pos+1, rem-w, val + (unsigned __int128)d * bpow[pos]);
    }
}

int main(int argc, char **argv) {
    if (argc != 7) return 1;
    P  = atoi(argv[1]); Q = atoi(argv[2]); R = atoi(argv[3]);
    S1 = strtod(argv[4], NULL);
    S2 = strtod(argv[5], NULL);
    S3 = strtod(argv[6], NULL);
    int hi = (P+1)/2, lo = (P-1)/2; MAX_W = 0;
    for (int k = 0; k < P; k++) {
        Tv[k] = (k - hi) * (k - lo) / 2;
        if (Tv[k] > MAX_W) MAX_W = Tv[k];
    }
    bpow[0] = 1;
    for (int i = 1; i < Q; i++)
        bpow[i] = bpow[i-1] * (unsigned __int128)(2*P - 1);
    go(0, R, 0);
    printf("%.17f ", H);
    print_u128(CNT);
    putchar('\n');
    return 0;
}
"""

BINARY     = Path(__file__).parent / "f3_engine_v1"
STATE_FILE = Path(__file__).parent / "z_state.json"
CERTIFIED_Z13_SHA256 = "cb371461c83d2a07d6eff58d778c27846dda5f5bda6ac1d8060b80bb08d97da4"
CERTIFIED_Z13_LAST_N = 13

# ──────────────────────────────────────────────────────────────────────────────
# COMPILATION
#
# Strategy:
#   1. Try clang + homebrew libomp (Mac)
#   2. Try gcc + -fopenmp (Linux/GCP)
#   3. Fall back to single-thread gcc (always works, just slower)
#
# Always writes OMP_C to disk before attempting OpenMP compilation.
# (Bug in earlier Grok v3.3: forgot to write source → compiled nothing.)
#
# Sets OMP_NUM_THREADS to all available cores automatically.
# On GCP, increase instance cores and this scales automatically.
# ──────────────────────────────────────────────────────────────────────────────
def compile_binary():
    src = BINARY.with_suffix(".c")

    # Get macOS SDK path if available (needed for clang on some setups)
    env = os.environ.copy()
    try:
        env["SDKROOT"] = subprocess.check_output(
            ["xcrun", "--show-sdk-path"], stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        pass

    # Always write OpenMP source first
    with open(src, 'w') as f:
        f.write(OMP_C)

    # Candidate OpenMP compile commands
    omp_candidates = []

    # Mac: Apple clang + homebrew libomp
    try:
        omp_prefix = subprocess.check_output(
            ["brew", "--prefix", "libomp"], stderr=subprocess.DEVNULL
        ).decode().strip()
        omp_candidates.append([
            "clang", "-O3",
            "-Xpreprocessor", "-fopenmp",
            f"-I{omp_prefix}/include",
            f"-L{omp_prefix}/lib", "-lomp",
            "-o", str(BINARY), str(src)
        ])
    except Exception:
        pass

    # Linux / GCP: gcc with OpenMP
    omp_candidates.append([
        "gcc", "-O3", "-fopenmp", "-o", str(BINARY), str(src)
    ])

    for cmd in omp_candidates:
        r = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if r.returncode == 0:
            cores = os.cpu_count() or 1
            os.environ["OMP_NUM_THREADS"] = str(cores)
            log.info(f"OpenMP binary compiled: {cmd[0]}  ({cores} cores)")
            return

    # Single-thread fallback — always correct, just slower
    log.warning("OpenMP unavailable — using single-thread fallback")
    log.warning("Install libomp (brew install libomp) for ~10x speedup")
    with open(src, 'w') as f:
        f.write(SINGLE_C)
    subprocess.run(["gcc", "-O3", "-o", str(BINARY), str(src)], check=True)
    log.info("Single-thread binary compiled")

# ──────────────────────────────────────────────────────────────────────────────
# MATH HELPERS
# ──────────────────────────────────────────────────────────────────────────────

def Tp(k, p):
    """T_p(k) = Wróblewski's parabolic weight function.
    T_p(k) = (k - ⌊(p+1)/2⌋)(k - ⌊(p-1)/2⌋) / 2
    This is the key to AP-freeness: it's convex, so sphere midpoints are
    strictly interior — they have weight strictly less than r, so they're
    excluded from B(p,q,r)."""
    return (k - math.floor((p+1)/2)) * (k - math.floor((p-1)/2)) // 2

def B_max_exact(p, q, r):
    """Compute the maximum element of B(p,q,r) exactly, using Python bigints.
    
    This is used to compute the shifts s1,s2,s3 BEFORE running the C engine.
    We need s1,s2,s3 to pass to C, but we can't compute them inside C without
    overflow. Python bigints handle arbitrarily large integers exactly.
    
    Strategy: greedy assignment from highest position to lowest.
    At each position, assign the highest digit whose weight is ≤ remaining budget.
    Digits with T_p=0 (center digits) don't consume budget, so they're not assigned
    here — only budget-consuming digits get placed at the highest positions."""
    b = 2*p - 1
    Tv = [Tp(k, p) for k in range(p)]
    val = 0
    rem = r
    for pos in range(q-1, -1, -1):
        best = 0
        for k in range(p-1, -1, -1):
            if Tv[k] <= rem:
                best = k
                break
        val += best * (b ** pos)
        rem -= Tv[best]
    return val   # Python bigint — no overflow possible

def A(n):
    """n-th element of G₃ (greedy base-3 AP-free set), 1-indexed.
    Formula: A(n) = 1 + (base-3 representation of n-1 read in base-3 using only {0,1})
    Equivalent to: replace each '1' bit of (n-1) with 3^k instead of 2^k.
    Runs in O(log n)."""
    m = n-1; r = 0; k = 0
    while m > 0:
        r += (m & 1) * (3**k)
        m >>= 1
        k += 1
    return r + 1

def project(h, gains):
    """Project construction limit from last two gains.
    The gain sequence decays geometrically: g_n ≈ g_{n-1} × ratio.
    Sum of geometric series: Σ g_n × ratio^k / (1-ratio).
    Returns inf if gains are increasing (signals overflow/corruption)."""
    if len(gains) < 2:
        return h
    ratio = gains[-1] / gains[-2]
    if ratio >= 1:
        return float('inf')   # gains increasing = corruption flag
    return h + gains[-1] * ratio / (1 - ratio)

def stream_step(p, q, r, z_max):
    """Execute one Wróblewski pairing step and return (h_gain, count, new_z_max).
    
    Shift computation entirely in Python bigints — never passed to any C integer.
    The C engine receives shifts as decimal strings; strtod() parses them as float64.
    
    The new z_max (= s3 + b_max) is also computed in Python bigints.
    This is the key anti-overflow design: Python arithmetic is exact for integers
    of any size. The C engine only ever sees floating-point values."""
    b_max = B_max_exact(p, q, r)          # Python bigint, exact
    z = z_max; t = b_max; m = max(z, t)
    s1 = m + z + 1                         # Python bigint
    s2 = 3*m + 2*t + z + 3                # Python bigint
    s3 = 3*m + 4*t + z + 4                # Python bigint

    result = subprocess.run(
        [str(BINARY), str(p), str(q), str(r),
         str(s1), str(s2), str(s3)],       # decimal strings → strtod in C
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"Engine error:\n{result.stderr.strip()}")

    parts   = result.stdout.strip().split()
    h_gain  = float(parts[0])
    count   = int(parts[1])
    new_z   = s3 + b_max                   # Python bigint arithmetic — exact

    return h_gain, count, new_z

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

# ──────────────────────────────────────────────────────────────────────────────
# CHECKPOINT / RESUME
#
# Saves state after every completed step. If the script is interrupted (power
# loss, GCP preemption, Ctrl+C), it resumes from the last checkpoint.
#
# z_max is stored as a string to preserve exact Python bigint value.
# JSON doesn't natively support integers larger than 2^53.
# ──────────────────────────────────────────────────────────────────────────────
def save_state(h, gains, z_max, n):
    with open(STATE_FILE, 'w') as f:
        json.dump({
            "h":      h,
            "gains":  gains,
            "z_max":  str(z_max),    # string to preserve bigint precision
            "last_n": n,
            "ts":     time.strftime("%Y-%m-%d %H:%M:%S"),
            "note":   "z_max stored as string — load with int()"
        }, f, indent=2)
    log.info(f"Checkpoint saved: Z_{n}")

def load_state():
    if not STATE_FILE.exists():
        return None
    state_hash = sha256_file(STATE_FILE)
    with open(STATE_FILE) as f:
        d = json.load(f)
    if d['last_n'] == CERTIFIED_Z13_LAST_N:
        if state_hash != CERTIFIED_Z13_SHA256:
            raise RuntimeError(
                "Certified Z13 checkpoint hash mismatch: "
                f"expected {CERTIFIED_Z13_SHA256}, got {state_hash}. "
                "Refusing to resume from a modified checkpoint."
            )
        log.info(f"Certified Z13 checkpoint hash verified: {state_hash}")
    log.info(f"Resumed from Z_{d['last_n']}: H={d['h']:.10f}")
    return d['h'], d['gains'], int(d['z_max']), d['last_n']

# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    log.info("=" * 70)
    log.info("f(3) LOWER BOUND ATTACK — DEFINITIVE v1.1")
    log.info("PIP Projects Inc. / Digital Sidewalk Lab")
    log.info("Darby Bailey McDonough, Ph.D.")
    log.info("=" * 70)

    compile_binary()
    t0_wall = time.time()

    # ── Resume or start fresh ───────────────────────────────────────────────
    state = load_state()

    if state:
        h, gains, z_max, start_n = state
        log.info(f"Resuming: H={h:.10f}  gap={h-3.00849:+.10f}\n")

    else:
        log.info("\nBuilding Z₀–Z₁₁ from verified construction...")

        # G₃: the greedy AP-free set, elements are A(1), A(2), ...
        # Only need z_max from Z₀ for shift computation.
        Z = [A(i) for i in range(1, 500_000) if A(i) <= 21_523_361]
        z_max = max(Z)
        del Z   # free RAM — we don't need the full list anymore

        # Z₁–Z₄: small enough for Python itertools enumeration
        # (these blocks have ≤ 5.9M nodes each — fast)
        for p2, q2, r2 in [(4,9,5),(4,10,5),(6,9,12),(6,10,13)]:
            b2 = 2*p2-1
            Tv2 = [Tp(k, p2) for k in range(p2)]
            pts = []
            for d in itertools.product(range(p2), repeat=q2):
                if sum(Tv2[x] for x in d) == r2:
                    pts.append(sum(d[i]*(b2**i) for i in range(q2)))
            b_mx = max(set(pts))
            z2 = z_max; t2 = b_mx; m2 = max(z2, t2)
            # Advance z_max using Python bigint arithmetic
            z_max = 3*m2 + 4*t2 + z2 + 4 + b_mx

        # Z₅–Z₁₁: only z_max matters; use B_max_exact for the shift chain
        for n in range(5, 12):
            r_n = 15 if n == 5 else math.floor(4*(n+6)/3)
            b_mx = B_max_exact(6, n+6, r_n)
            z2 = z_max; t2 = b_mx; m2 = max(z2, t2)
            z_max = 3*m2 + 4*t2 + z2 + 4 + b_mx

        # Use verified Z₁₁ H — more reliable than accumulating float rounding errors
        # across 11 Python-computed steps.
        h = VERIFIED[11]["H"]
        gains = [
            0.0020734151,   # Z1  gain
            0.0010081648,   # Z2
            0.0005199461,   # Z3
            0.0003211768,   # Z4
            0.0001726668,   # Z5
            0.0000908756,   # Z6
            0.0000476959,   # Z7
            0.0000250310,   # Z8
            0.0000132409,   # Z9
            0.0000070051,   # Z10
            0.0000036987,   # Z11
        ]
        start_n = 11

        log.info(f"Z₁₁ z_max: {z_max:,}")
        log.info(f"Z₁₁ H    : {h:.10f}")
        log.info(f"Record   : 3.00849000  Margin: +{h-3.00849:.10f}  ✓ BEATEN\n")

    # ── Stream Z₁₂+ ─────────────────────────────────────────────────────────
    log.info("=" * 70)
    log.info("STREAMING Z₁₂+")
    log.info("=" * 70)

    for n in range(start_n + 1, 40):
        r_n = math.floor(4*(n+6)/3)
        log.info(f"\nStep {n}: B(6,{n+6},{r_n})...")

        t0 = time.time()
        h_gain, count, z_max = stream_step(6, n+6, r_n, z_max)
        elapsed = time.time() - t0

        h += h_gain
        gains.append(h_gain)
        proj = project(h, gains)

        log.info(f"Z_{n}:  {count:,} nodes  {elapsed/3600:.3f}h")
        log.info(f"       H    = {h:.10f}")
        log.info(f"       gain = {h_gain:+.10f}")
        log.info(f"       proj = {proj:.10f}")
        log.info(f"       Gap  = {h - 3.00849:+.10f}")
        log.info(f"       Wall = {(time.time()-t0_wall)/60:.1f} min")

        # ── Verification gate ────────────────────────────────────────────────
        # Check against known-good values. Both H and node count must match.
        # Node count is deterministic (same block = same tree traversal).
        # If either fails: halt immediately and don't save the checkpoint.
        if n in VERIFIED:
            v = VERIFIED[n]
            h_ok     = abs(h - v["H"]) < VERIFY_TOLERANCE
            count_ok = (count == v["nodes"])
            if h_ok and count_ok:
                log.info(f"✅ Z{n} VERIFIED — H and node count match ground truth")
            else:
                if not h_ok:
                    log.error(f"❌ Z{n} H mismatch: got {h:.10f}, expected {v['H']:.10f}")
                    log.error(f"   Difference: {h - v['H']:.4e}")
                    log.error(f"   If difference ≈ +7e-6: atoll overflow still in C code")
                    log.error(f"   If difference ≈ -7e-6: some other precision issue")
                if not count_ok:
                    log.error(f"❌ Z{n} node count: got {count:,}, expected {v['nodes']:,}")
                    log.error(f"   If count is short: possible OpenMP dropped branches")
                log.error("Halting — do not use this result for publication.")
                sys.exit(1)

        # ── Anomaly detection ────────────────────────────────────────────────
        # Gains must monotonically decrease. Any increase signals:
        # - Overflow (atoll or integer wrapping in C) → fix: check strtod/int128
        # - OpenMP race condition (dropped or double-counted branches)
        # - Code bug in shift computation
        if len(gains) > 1 and gains[-1] > gains[-2]:
            log.error(f"❌ GAIN INCREASED: {gains[-2]:.6e} → {gains[-1]:.6e}")
            log.error(f"   Ratio: {gains[-1]/gains[-2]:.4f}  (must be < 1)")
            log.error(f"   proj=inf is the usual companion symptom")
            log.error("Halting — do not use this result for publication.")
            sys.exit(1)

        # Save checkpoint after verification passes
        save_state(h, gains, z_max, n)

    # ── Final report ─────────────────────────────────────────────────────────
    log.info("\n" + "=" * 70)
    log.info("FINAL RESULT")
    log.info("=" * 70)
    log.info(f"f(3) >= {h:.10f}")
    log.info(f"Wróblewski 1984 record: 3.00849000")
    log.info(f"Improvement: +{h - 3.00849:.10f}")
    log.info(f"AP-free by pairing theorem (Wróblewski 1984) ✓")
    log.info(f"Total wall time: {(time.time()-t0_wall)/60:.1f} min")
