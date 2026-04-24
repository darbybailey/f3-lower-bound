# DEVNOTES.md — Engineering Notes

Development history, overflow bugs found and fixed, and deployment notes.
The core math and construction are in the paper; this file documents
the implementation journey.

---

## Bug History

### BUG 1 — atoll() Integer Overflow (Critical)

**Symptom:** Z₁₂ gain appeared plausible (+0.0000026), but Z₁₃ gain exploded to
+0.000041 with proj=inf.

**Root cause:** Shifts s₂ ≈ 1.61×10¹⁹ and s₃ ≈ 2.17×10¹⁹ at Z₁₂ exceed
LLONG_MAX = 9.22×10¹⁸. C's `atoll()` silently wraps them to large negative numbers.
`1/(val + negative)` produces wrong-signed, inflated values.

**Insidious detail:** The Z₁₂ result from `atoll()` was +0.0000026404 — plausible
enough to fool us for a week. The wrapping happened to produce a number in the right
ballpark rather than obviously exploding. Only Z₁₃ made the corruption obvious.

**Fix:** Replace `atoll()` with `strtod()`. Shifts passed as decimal strings from
Python; `strtod()` reads them as float64 — sufficient for 1/(val+s) at these scales.

**Lesson:** Any C code parsing shifts as 64-bit integers will silently fail at Z₁₂+.
Always use `strtod()`. Always pass shifts as decimal strings from Python bigints.

---

### BUG 2 — long / double Overflow for Node Values (Critical)

**Symptom:** Results at Z₁₃ would be wrong even after fixing Bug 1.

**Root cause:** Node values at Z₁₃ reach ~3×10¹⁹. float64 represents integers exactly
only to 2⁵³ ≈ 9×10¹⁵. `long` overflows at 9.22×10¹⁸. Both fail silently.

**Fix:** `unsigned __int128` for all internal arithmetic (`bpow[]`, `val`). Cast to
`double` only at the final `1.0/(v+s)` step. `__int128` is exact to ~1.7×10³⁸.

**Lesson:** Any code using `long` or `double` for node values fails at Z₁₃+.
`unsigned __int128` is mandatory, not optional.

---

### BUG 3 — Misdiagnosed Race Condition (Developer Error)

**Symptom:** The first correct OpenMP Z₁₂ result (H=3.0084948386) was lower than the
"verified" single-core result (H=3.0084955124). We concluded a race condition was
dropping 25% of contributions.

**True cause:** The "verified" single-core result was itself corrupted by Bug 1.
The atoll overflow on s₂ and s₃ happened to inflate H into a plausible-looking range.

**Evidence that there was no race condition:** Both runs produced exactly
7,563,706,368,000 nodes. If a race condition dropped branches, node counts would differ.

**The `dummy=0.0` Kahan issue:** The concern was that using a fresh `kc=0.0` in the
OpenMP critical section loses inter-thread Kahan compensation. This is true, but the
error is ~1.6×10⁻²² — negligible. It is not a correctness issue.

**Lesson:** Always verify the reference value before diagnosing a deviation from it.
"Plausible-looking" results can still be corrupted.

---

### BUG 4 — Signed 64-bit Node Counter Overflow Horizon

**Symptom:** Not triggered in the certified Z11–Z13 runs, but the original `long long`
node counter becomes unsafe if the script is extended far enough.

**Root cause:** Node counts grow by about 5.8x per step. Starting from
44,239,186,034,688 nodes at Z13, a signed 64-bit counter crosses `2^63-1`
around Z20. Signed overflow in C is undefined behavior, so later counts would
eventually become unreliable even if the harmonic sum arithmetic remained correct.

**Fix:** Promote node counters in both C kernels to `unsigned __int128` and print
them with an explicit decimal formatter.

**Attribution:** Identified by ChatGPT Codex, April 2026.

**Lesson:** Value arithmetic and metadata arithmetic need independent overflow
budgets. `__int128` was already mandatory for node values; it is also the right
choice for long-horizon node counts.

---

## Ground Truth (Two Independent Correct Runs, April 2026)

| Step | H | gain | nodes |
|---|---|---|---|
| Z₁₁ | 3.0084928720 | +0.0000036987 | 1,293,635,551,232 |
| Z₁₂ | 3.0084948386 | +0.0000019666 | 7,563,706,368,000 |
| Z₁₃ | 3.0084958842 | +0.0000010455 | 44,239,186,034,688 |

Note: an earlier Z₁₂ result of +0.0000026404 was **corrupted by Bug 1** (atoll overflow
on s₂ and s₃). The correct Z₁₂ gain is +0.0000019666.

---

## Script History

| Script | Status | Notes |
|---|---|---|
| `claude_win_c.py` | Superseded | Original C streamer, atoll overflow at Z₁₂+ |
| `script_z_cstream.py` | Superseded | atoll overflow |
| `script_godmode.py` | Superseded | atoll overflow |
| `gemini_ftw.py` | Reference | First correct Z₁₂ (strtod + __int128 + OpenMP) |
| `z_push_v4.1.py` | Reference | Grok's version, all fixes, clang + collapse(4) |
| **`f3_lower_bound.py`** | **Current** | All lessons, clean header, Zenodo/GitHub ready |

---

## GCP Deployment Notes

The script auto-detects Linux and compiles with `gcc -fopenmp`. On a 176-core instance:

```bash
sudo apt-get install -y python3 gcc libomp-dev
python3 f3_lower_bound.py
```

`OMP_NUM_THREADS` is set automatically from `os.cpu_count()`. Scale cores by changing
the instance type — no code changes needed.

Approximate wall times scale as (Z₁₂ time) × 5.8^(n-12) per step.
Z₁₄ on 176 cores: ~15–20 hours.

Recommended: Spot/Preemptible instance. The `z_state.json` checkpoint means a preemption
costs at most one step's compute time.

---

## Overflow Reference Table

| Step | s₁ | s₂ | s₃ | atoll safe? |
|---|---|---|---|---|
| Z₁₁ | 4.55×10¹⁷ | 1.47×10¹⁸ | 1.97×10¹⁸ | ✓ all safe |
| Z₁₂ | 5.00×10¹⁸ | **1.61×10¹⁹** | **2.17×10¹⁹** | s₁ ✓, s₂/s₃ OVERFLOW |
| Z₁₃ | 5.50×10¹⁹ | 1.77×10²⁰ | 2.39×10²⁰ | ALL OVERFLOW |
