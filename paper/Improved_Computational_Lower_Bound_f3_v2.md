---
title: "An Improved Computational Lower Bound for $f(3)$ via Wróblewski's Recursive Pairing of Behrend Sphere Blocks"
author: "Darby Bailey McDonough"
orcid: "0000-0003-4072-6191"
affiliation: "Pip Projects, Inc. · Digital Sidewalk Lab · Independent Researcher"
doi: "10.5281/zenodo.24486265"
concept-doi: "10.5281/zenodo.19722201"
v1-doi: "10.5281/zenodo.19722202"
version: "2.0"
date: "August 2026"
license: "CC BY 4.0"
---

# An Improved Computational Lower Bound for $f(3)$ via Wróblewski's Recursive Pairing of Behrend Sphere Blocks

**August 2026 revision:** analytic Jensen extension of the pairing chain through $Z_{32}$

Darby Bailey McDonough\
Pip Projects, Inc. · Digital Sidewalk Lab · Independent Researcher\
ORCID: 0000-0003-4072-6191

Preprint — April 2026, revised August 2026 (v2)\
CC BY 4.0

Correspondence: ORCID 0000-0003-4072-6191\
This version DOI: 10.5281/zenodo.24486265\
Canonical v1 DOI: 10.5281/zenodo.19722202 (concept 10.5281/zenodo.19722201)\
GitHub: https://github.com/darbybailey/f3-lower-bound/

Keywords: additive combinatorics; arithmetic progressions; harmonic sums; certified computation; Behrend construction; Wróblewski pairing; Jensen inequality; exact rational certificate

---

## Abstract

We report a revised lower bound
$$f(3)\ge 3.008500858304852623,$$
obtained from a finite explicit $3$-AP-free subset of the positive integers. The April 2026 version of this paper established the certified computational bound
$$f(3)\ge 3.0084928720$$
at the pairing step $Z_{11}$, by a C+Python enumeration of Wróblewski's recursive pairing of Behrend sphere blocks, together with an IEEE 754 error certificate. Supporting streamed values at $Z_{12}$ and $Z_{13}$ were $H=3.0084948386$ and $H=3.0084958842$. The construction in that version used $p=6$ through $Z_{13}$.

The present August 2026 revision retains the $Z_{11}$ theorem and records two subsequent extensions of the same pairing construction. First, an interval-certified adaptive chain $Z_{14}$--$Z_{23}$ (certificate of 24 July 2026) yields
$$f(3)\ge 3.008500752707.$$
An outward-rounded moment-certificate checkpoint at $Z_{24}$, consumed as exact integers, gives the machine-interval lower
$$H\ge 3.008500791468429.$$
There is no $Z_{24}$ JSON file in the packet; $Z_{24}$ is not a Jensen step. Second, eight exact-rational Jensen steps $Z_{25}$--$Z_{32}$ on even-$p$ blocks produce the $18$-place downward truncation above. No floating-point value enters the arithmetic of that last theorem.

The number is a finite-construction certificate. It is not a proof that the block schedule is globally optimal. The block search used local/IMR-style continuation screening over scale-compatible $(p,q,r)$ families; that screening is a search method only and is not imported as a theorem.

---

## 1. The problem

For $k\ge 3$, define
$$
f(k)=\sup\Bigl\{\sum_{n\in A}\frac1n:A\subseteq\mathbb{N}^+,\ A\text{ is $k$-AP-free}\Bigr\}.
$$
A set $A\subseteq\mathbb{N}^+$ is $k$-AP-free if it contains no $k$-term arithmetic progression. The finiteness of $f(3)$ was proved by Bloom and Sisask (2020) as a consequence of their density bounds on $3$-AP-free sets; the explicit value of $f(3)$ remains unknown.

The best published numerical lower bound prior to the April 2026 version of this paper was $f(3)\ge 3.00849$ (Wróblewski 1984, *Mathematics of Computation* 43, 261–262). That paper used the recursive pairing of Behrend sphere blocks which is followed here.

The April 2026 theorem of the present work is $f(3)\ge 3.0084928720$, a strict improvement of Wróblewski's printed bound, certified at $Z_{11}$. The August 2026 revision does not withdraw that theorem. It extends the same pairing to a longer explicit chain and replaces later streamed floating-point increments by interval-certified and then exact-rational lower bounds.

---

## 2. Construction

We follow Wróblewski (1984). The two lemmas below are the same ingredients used for the $p=6$ chain in the April text; they are written here at the generality required by the later adaptive and Jensen steps. The pairing theorem itself is classical. The interval form and the even-$p$ mean identity are recorded because they are the statements used by the implementation.

### 2.1. Behrend sphere blocks

For integers $p\ge 2$ and $0\le k<p$, define the weight
$$
T_p(k)=\frac{\bigl(k-\lfloor(p+1)/2\rfloor\bigr)\bigl(k-\lfloor(p-1)/2\rfloor\bigr)}{2}.
$$
For integers $q\ge 1$ and $r\ge 0$, the Behrend block $B(p,q,r)$ is the set of integers
$$
B(p,q,r)=\Bigl\{\sum_{i=0}^{q-1}d_i(2p-1)^i:0\le d_i<p,\ \sum_{i=0}^{q-1}T_p(d_i)=r\Bigr\}.
$$
Write $b=2p-1$ for the underlying base. When $p$ is even and $h=p/2$, the weights
$$
w_p(k)=\frac{(k-h)(k-h+1)}{2},\qquad 0\le k<p,
$$
coincide with $T_p$. This identity is used by the $Z_{25}$--$Z_{32}$ tables, all of which have even $p$; it was checked directly for $p\in\{4,6,8,64,80\}$.

**Lemma 1 (Block AP-free).** For all integers $p\ge 2$, $q\ge 1$, and $r\ge 0$, the set $B(p,q,r)$ is $3$-AP-free.

*Proof.* Let $\alpha=\lfloor(p+1)/2\rfloor$ and $\beta=\lfloor(p-1)/2\rfloor$, so that $2T_p(k)=(k-\alpha)(k-\beta)$. Expanding gives
$$
2T_p(k)=k^2-(\alpha+\beta)k+\alpha\beta.
$$
Hence, for any integers $x,y,z$,
$$
2\bigl(T_p(x)+T_p(z)-2T_p(y)\bigr)=\bigl(x^2+z^2-2y^2\bigr)-(\alpha+\beta)(x+z-2y).
$$
If $x+z=2y$, the linear term vanishes and
$$
x^2+z^2-2y^2=\frac{(x-z)^2}{2},
$$
so
$$
T_p(x)+T_p(z)-2T_p(y)=\frac{(x-z)^2}{4}.
$$
The right-hand side is a non-negative integer whenever $x+z=2y$ (equivalently, $x-z$ is even). This identity was additionally checked as an integer equality on $20136$ triples with $x+z=2y$ and $0\le x,y,z<p$, including the odd values $p=5,7$ and the even values $p\in\{4,6,8,64,80,172\}$.

Now suppose $x,y,z\in B(p,q,r)$ satisfy $x+z=2y$ as integers, and write each in base $b=2p-1$ with digits in $\{0,\dots,p-1\}$. Digitwise sums satisfy $x_i+z_i\le 2(p-1)=b-1$, so the addition $x+z$ has no carries. Therefore $x+z=2y$ holds digitwise: $x_i+z_i=2y_i$ for every position $i$. Summing the second-difference identity over digits yields
$$
\sum_i\bigl(T_p(x_i)+T_p(z_i)-2T_p(y_i)\bigr)=\sum_i\frac{(x_i-z_i)^2}{4}.
$$
The left-hand side is $r+r-2r=0$ because $x,y,z$ lie on the same level. Each summand on the right is non-negative, so each vanishes, so $x_i=z_i$ for all $i$, so $x=z$, and then $x=y=z$. Thus $B(p,q,r)$ contains no $3$-term arithmetic progression. $\square$

The April text stated the same fact as a convexity proposition: $T_p$ is strictly convex, and for $a\neq c$ in $B$ the midpoint fails the level constraint. The argument above is that proposition written in coordinates, and it does not use evenness of $p$.

### 2.2. Pairing

For finite nonempty $3$-AP-free sets $Z,T\subset\mathbb{N}^+$ write $z=\max Z$, $t=\max T$, and $m=\max(z,t)$. Set
$$
s_1=m+z+1,\qquad s_2=3m+2t+z+3,\qquad s_3=3m+4t+z+4,
$$
and define
$$
\mathrm{pair}(Z,T)=Z\cup(T+s_1)\cup(T+s_2)\cup(T+s_3).
$$
The four components are contained in the intervals
\begin{align*}
I_0&=[1,z],\\
I_1&=[s_1+1,s_1+t]=[m+z+2,\,m+z+t+1],\\
I_2&=[s_2+1,s_2+t]=[3m+2t+z+4,\,3m+3t+z+3],\\
I_3&=[s_3+1,s_3+t]=[3m+4t+z+5,\,3m+5t+z+4].
\end{align*}
These are the containing intervals used by the implementation (and by the later normalized chain, in which each appended block is translated so that its minimum is $1$).

**Lemma 2 (Pairing).** If $Z$ and $T$ are $3$-AP-free, then $\mathrm{pair}(Z,T)$ is $3$-AP-free.

This is Wróblewski (1984, Theorem 1). Each translate $T+s_i$ is $3$-AP-free because $T$ is. The shifts $s_1<s_2<s_3$ are chosen so that a $3$-term progression with terms in two or more distinct components cannot close inside the union. The full case analysis is in Wróblewski (1984). The present paper does not replace that analysis by a computer search.

The interval hull above is recorded because it is the form used in code. If the four intervals admit no mixed $3$-term progression, neither do any subsets placed in those intervals. An interval-hull audit of the implemented shifts found no overlap, no configuration with two terms in one interval and the third in another, and no three-block mixed progression, on a fixed list of pairs including $(z,t)=(64,2512)$ and on $400$ random pairs. That audit is a supporting check of the implemented formulae. It is not a proof of Lemma 2.

**Consequence.** Starting from a $3$-AP-free seed and iteratively pairing with Behrend blocks, the entire chain is $3$-AP-free by induction on Lemmas 1 and 2. No computational search for progressions in the constructed set is required or performed.

### 2.3. Harmonic recurrence

If $Z_{n+1}=\mathrm{pair}(Z_n,B_{n+1})$ and $B_{n+1}$ has been translated into $\{1,\dots,t\}$, the harmonic sum obeys
$$
H(Z_{n+1})=H(Z_n)+\sum_{b\in B_{n+1}}\Bigl(\frac1{b+s_1}+\frac1{b+s_2}+\frac1{b+s_3}\Bigr).
$$
The April computation evaluated the increment by enumerating $B_{n+1}$. The July certificate lower-bounds the increment by an outward-rounded moment expansion. The August steps $Z_{25}$--$Z_{32}$ lower-bound it by Jensen's inequality at the exact mean, in $\mathbb{Q}$.

### 2.4. Even-$p$ Jensen increment

**Proposition 1 (even-$p$ Jensen).** Let $p$ be even, $b=2p-1$, and let $B=B(p,q,r)$ be nonempty. Complementing digits by $d\mapsto p-1-d$ preserves $T_p$ and sends a raw block value $v$ to $(b^q-1)/2-v$. Consequently the mean raw value on the level is
$$
\mathbb{E}[v]=\frac{b^q-1}{4}.
$$
This mean identity was checked by complete enumeration on the slice $p=8$, $q=4$, $r=6$ ($n=224$), where both sides equal $12656$. After translating $B$ so that its minimum is $1$, write $\ell$ for the raw minimum and $\mathrm{off}=\ell-1$, and set
$$
\mu=\frac{b^q-1}{4}-\mathrm{off}.
$$
The map $x\mapsto 1/x$ is convex on $(0,\infty)$. Jensen's inequality therefore gives, for each shift $S>0$,
$$
\sum_{v\in B}\frac1{S+v}\ge\frac{|B|}{S+\mathbb{E}[v]}=\frac{n}{S+\mu},
$$
where $n=|B|$ and the expectation is taken after the same translation. Applying the inequality to $S\in\{s_1,s_2,s_3\}$ produces the exact rational lower bound
$$
\Delta H\ge n\sum_{j=1}^{3}\frac1{s_j+\mu}\in\mathbb{Q}.
$$
No floating-point value is used in this arithmetic.

The hypothesis that $p$ is even is essential for the complement symmetry as stated. Odd $p$ is not used in $Z_{25}$--$Z_{32}$. Lemma 1 itself does not require even $p$.

---

## 3. Computational implementation, $Z_0$--$Z_{13}$

This section is the April 2026 computational record, condensed. It is the basis of Theorem 1 below and is not superseded by later interval or Jensen bounds.

### 3.1. Block sequence through $Z_{11}$

The seed $Z_0$ is the standard greedy $3$-AP-free set
$$
G_3=\{n\in\mathbb{N}^+:\text{the base-$3$ digits of $n-1$ lie in $\{0,1\}$}\},
$$
taking those elements that are at most $21{,}523{,}361$ ($65{,}536$ nodes). Subsequent blocks and node counts are as follows.

| Step | Block | Nodes |
|---|---|---|
| $Z_0$ | $G_3$, elements $\le 21{,}523{,}361$ | $65{,}536$ |
| $Z_1$ | $B(4,9,5)$ | $64{,}512$ |
| $Z_2$ | $B(4,10,5)$ | $258{,}048$ |
| $Z_3$ | $B(6,9,12)$ | $1{,}053{,}696$ |
| $Z_4$ | $B(6,10,13)$ | $5{,}959{,}680$ |
| $Z_5$ | $B(6,11,15)$ | $34{,}174{,}976$ |
| $Z_6$ | $B(6,12,16)$ | $197{,}345{,}280$ |
| $Z_7$ | $B(6,13,17)$ | $1{,}139{,}294{,}208$ |
| $Z_8$ | $B(6,14,18)$ | $6{,}578{,}044{,}928$ |
| $Z_9$ | $B(6,15,20)$ | $38{,}260{,}604{,}928$ |
| $Z_{10}$ | $B(6,16,21)$ | $222{,}697{,}619{,}456$ |
| $Z_{11}$ | $B(6,17,22)$ | $1{,}293{,}635{,}551{,}232$ |

For $n\ge 3$ the $p=6$ radius follows $r_n=\lfloor 4(n+6)/3\rfloor$ with $r_5=15$ (Wróblewski 1984).

### 3.2. Overflow hazards and the C kernel

Three numerical hazards arise at $Z_{12}$ and beyond and do not affect $Z_0$--$Z_{11}$. At $Z_{12}$ the shifts $s_2\approx 1.61\times 10^{19}$ and $s_3\approx 2.17\times 10^{19}$ exceed the $64$-bit signed maximum $9.22\times 10^{18}$. All shift and $z_{\max}$ arithmetic is performed in Python arbitrary-precision integers; shifts are passed to the C kernel as decimal strings and parsed with `strtod()`. Node values at $Z_{13}$ reach approximately $3\times 10^{19}$, beyond exact double representation ($2^{53}\approx 9\times 10^{15}$) and beyond $64$-bit `long`. Internal C arithmetic uses `unsigned __int128`, with a cast to `double` deferred to the final $1.0/(v+s)$ step. With $12$ OpenMP threads, each accumulates a thread-local Kahan sum; the final merge introduces at most approximately $1.6\times 10^{-22}$ absolute error, which is negligible at $10$ decimal places.

Core enumeration (simplified):

```
void go(int pos, int rem, unsigned __int128 val,
        double *H, double *kc, unsigned __int128 *cnt) {
    if (pos == Q) {
        if (rem == 0) {
            (*cnt)++;
            double v = (double)val;    /* sole cast point */
            double x = 1.0/(v+S1)+1.0/(v+S2)+1.0/(v+S3) - *kc;
            double t = *H + x;
            *kc = (t - *H) - x;       /* Kahan compensation */
            *H = t;
        }
        return;
    }
    if (rem > (Q-pos)*MAX_W) return;
    for (int d = 0; d < P; d++) {
        if (Tv[d] > rem) continue;
        go(pos+1, rem-Tv[d],
           val + (unsigned __int128)d * bpow[pos], H, kc, cnt);
    }
}
```

Shifts `S1`, `S2`, and `S3` are read via `strtod()`. `unsigned __int128` is available on GCC $\ge 4.6$ and Clang $\ge 3.1$ on standard $64$-bit platforms.

### 3.3. Floating-point error certificate at $Z_{11}$

Let $N=3\times 1{,}293{,}635{,}551{,}232\approx 3.88\times 10^{12}$ be the number of floating-point additions in the $Z_{11}$ step, let $\varepsilon=2^{-53}\approx 1.11\times 10^{-16}$ be IEEE 754 double-precision machine epsilon, and let $\max|t_i|\le 1/s_1^{\min}\approx 2.20\times 10^{-18}$ be the largest possible term magnitude. Denominator rounding contributes at most $\delta_{\mathrm{den}}\le\varepsilon N\max|t_i|$. Summation error, by standard backward analysis (Higham 2002, Theorem 4.1), contributes at most $\delta_{\mathrm{sum}}\le 2\varepsilon N\max|t_i|$. Combined,
$$
\delta H\le 3\varepsilon N\max|t_i|\approx 2.85\times 10^{-21}.
$$
The margin $\Delta=3.0084928720-3.00849=2.872\times 10^{-6}$ exceeds this conservative bound by a factor exceeding $10^{15}$. Therefore $H(Z_{11})>3.00849$ is numerically certified under IEEE 754 double-precision arithmetic.

### 3.4. Verification protocol and streamed results

At $Z_{12}$, the first step requiring overflow safeguards, two independent implementations (single-threaded and $12$-core OpenMP) were required to agree on the exact node count $7{,}563{,}706{,}368{,}000$ and on $H=3.0084948386\pm 10^{-8}$ before continuation. Node-count agreement rules out branch-dropping; agreement in $H$ rules out silent overflow.

| Step | $H$ | Gain | Nodes | Gain ratio |
|---|---|---|---|---|
| $Z_0$ | $3.0042100$ | — | $65{,}536$ | — |
| $Z_1$ | $3.0062834$ | $+0.0020734$ | $64{,}512$ | — |
| $Z_2$ | $3.0072915$ | $+0.0010082$ | $258{,}048$ | $0.4864$ |
| $Z_3$ | $3.0078115$ | $+0.0005199$ | $1{,}053{,}696$ | $0.5156$ |
| $Z_4$ | $3.0081327$ | $+0.0003212$ | $5{,}959{,}680$ | $0.6178$ |
| $Z_5$ | $3.0083053$ | $+0.0001727$ | $34{,}174{,}976$ | $0.5376$ |
| $Z_6$ | $3.0083962$ | $+0.0000909$ | $197{,}345{,}280$ | $0.5261$ |
| $Z_7$ | $3.0084439$ | $+0.0000477$ | $1{,}139{,}294{,}208$ | $0.5247$ |
| $Z_8$ | $3.0084689$ | $+0.0000250$ | $6{,}578{,}044{,}928$ | $0.5249$ |
| $Z_9$ | $3.0084822$ | $+0.0000132$ | $38{,}260{,}604{,}928$ | $0.5296$ |
| $Z_{10}$ | $3.0084892$ | $+0.0000070$ | $222{,}697{,}619{,}456$ | $0.5306$ |
| $Z_{11}$ | $3.0084928720$ | $+0.0000037$ | $1{,}293{,}635{,}551{,}232$ | $0.5312$ |
| $Z_{12}$ | $3.0084948386$ | $+0.0000020$ | $7{,}563{,}706{,}368{,}000$ | $0.5317$ |
| $Z_{13}$ | $3.0084958842$ | $+0.0000010$ | $44{,}239{,}186{,}034{,}688$ | $0.5317$ |

The bound $f(3)\ge 3.0084928720$ is established at $Z_{11}$. Steps $Z_{12}$ and $Z_{13}$ are supporting computations. The April text noted that the $p=6$ gain ratio had stabilized near $0.5317$ and recorded the geometric extrapolation $H_\infty\approx H(Z_{11})+g_{12}/(1-0.5317)\approx 3.0084970$ as a ceiling for that $p=6$ schedule, not as an upper bound on $f(3)$. The later adaptive chain uses larger even $p$ and is not constrained by that $p=6$ figure.

Wall times on an Apple MacBook Pro (Apple M3, $12$-core OpenMP): $Z_{11}$ $6.0$ h (single-threaded C, March 2026); $Z_{12}$ $2.619$ h; $Z_{13}$ $15.1$ h, completed $2026$-$04$-$05$ $02$:$49$:$56$ as recorded in the April text. The $Z_{13}$ checkpoint hash recorded in the April deposit is
`cb371461c83d2a07d6eff58d778c27846dda5f5bda6ac1d8060b80bb08d97da4`.
The streamed $Z_{12}$ gain was $+0.0000019666$; the streamed $Z_{13}$ gain was $+0.0000010455$.


---

## 4. Adaptive extension, $Z_{14}$--$Z_{24}$

Enumeration of $p=6$ blocks becomes impractical past $Z_{13}$. The 24 July 2026 certificate continues the pairing by an interval-certified adaptive schedule: each appended block is a Behrend sphere $B(p,q,r)$ with even $p$ allowed to grow, counts are dual-verified (digit dynamic programming against a generating-function convolution), and endpoints are printed by integer truncation — floor for lowers, ceiling for uppers — never by floating-point division. The chain through $Z_{13}$ in that certificate is a normalized pairing (each block translated so that its minimum is $1$). Those $12$-decimal cumulative lowers are a separate computation from the April IEEE table in §3 and are not a replacement of it. From $Z_{14}$ onward the certificate is the record used here.

### 4.1. Compact Jensen check from the published $Z_{11}$ value

The same certificate records a compact, hand-auditable increment from the published April value $3.0084928720$, using exact-rational Jensen bounds on three blocks:

| Block | Count | Gain lower bound |
|---|---|---|
| $B(6,18,24)$ | $7{,}563{,}706{,}368{,}000$ | $1.9426354533\times 10^{-6}$ |
| $B(6,19,25)$ | $44{,}239{,}186{,}034{,}688$ | $1.0329304755\times 10^{-6}$ |
| $B(28,12,372)$ | $905{,}700{,}370{,}477{,}056$ | $1.7427784961\times 10^{-6}$ |

This is Certificate Theorem 1: $f(3)\ge 3.0084975903$. It is an intermediate check, not the headline of the adaptive chain.

### 4.2. Interval-certified adaptive chain $Z_{14}$--$Z_{23}$

The following schedule is copied from the per-step integer certificate. Every listed count is dual-verified. Cumulative lowers are the certificate's $12$-decimal floor truncations.

| Step | Block | Count | $\mathrm{cum\_lower}$ ($12$ d.p.) |
|---|---|---|---|
| $Z_{14}$ | $B(28,12,372)$ | $905{,}700{,}370{,}477{,}056$ | $3.008497699557$ |
| $Z_{15}$ | $B(36,12,633)$ | $11{,}213{,}058{,}566{,}356{,}992$ | $3.008498838757$ |
| $Z_{16}$ | $B(44,12,930)$ | $83{,}324{,}249{,}853{,}149{,}184$ | $3.008499492524$ |
| $Z_{17}$ | $B(52,12,1320)$ | $442{,}756{,}477{,}229{,}973{,}504$ | $3.008499872235$ |
| $Z_{18}$ | $B(44,13,1026)$ | $3{,}525{,}553{,}663{,}296{,}258{,}048$ | $3.008500179465$ |
| $Z_{19}$ | $B(52,13,1431)$ | $22{,}112{,}064{,}204{,}648{,}570{,}880$ | $3.008500384418$ |
| $Z_{20}$ | $B(44,14,1104)$ | $149{,}432{,}474{,}571{,}740{,}086{,}272$ | $3.008500529887$ |
| $Z_{21}$ | $B(52,14,1531)$ | $1{,}108{,}124{,}557{,}322{,}905{,}321{,}472$ | $3.008500634222$ |
| $Z_{22}$ | $B(60,14,2052)$ | $6{,}171{,}384{,}832{,}138{,}164{,}387{,}840$ | $3.008500699504$ |
| $Z_{23}$ | $B(52,15,1654)$ | $55{,}717{,}204{,}865{,}228{,}008{,}488{,}960$ | $3.008500752707$ |

This is Certificate Theorem 2:
$$f(3)\ge 3.008500752707.$$
Assertions recorded as active at every step of that certificate: dual independent counts; stratified zeroth moments equal to the count; positive expansion centres; the even-$p$ symmetry lemma. The pairing shifts and terminal maxima are stored as exact integers in the certificate; they are omitted from the body.

Already at $Z_{14}$ the $12$-decimal lower $3.008497699557$ lies above the April $p=6$ geometric ceiling $3.0084970$. That is the expected behaviour once $p$ is allowed to grow; it does not affect the April $Z_{11}$ theorem.

### 4.3. The $Z_{24}$ checkpoint

The $29$ August 2026 analytic-extension packet continues from a $Z_{24}$ checkpoint, not from a $Z_{24}$ Jensen JSON file. There is no `Z24_jensen.json` in the packet, and none is claimed.

The checkpoint is an outward-rounded moment-certificate state consumed as exact integers by the verifier `verify_frontier.py`:

- block $B(64,15,2512)$;
- terminal maximum
  $$z=151034783145739756303329159115181;$$
- integer lower endpoint
  $$\ell_{24}=30085007914684292896021976970576819728782,$$
  used as the rational $\ell_{24}/10^{40}$;
- machine-interval lower, $15$ decimal places,
  $$H\ge 3.008500791468429.$$

The $15$-decimal figure is the packet's published truncation of that rational. The theorem in §6 uses the integer pair $(\ell_{24},z)$, not a reconstructed floating-point $H$. The step from $Z_{23}$ to $Z_{24}$ is part of the moment-certificate chain; it is not one of the eight exact-rational Jensen steps.

---

## 5. Analytic Jensen extension, $Z_{25}$--$Z_{32}$

### 5.1. The eight even-$p$ blocks

From the $Z_{24}$ integers the packet appends eight explicit even-$p$ sphere blocks, each paired by Lemma 2 and lower-bounded by Proposition 1:

| Step | Block |
|---|---|
| $Z_{25}$ | $B(80,15,3910)$ |
| $Z_{26}$ | $B(98,15,5878)$ |
| $Z_{27}$ | $B(122,15,9151)$ |
| $Z_{28}$ | $B(108,16,7633)$ |
| $Z_{29}$ | $B(132,16,11464)$ |
| $Z_{30}$ | $B(164,16,17779)$ |
| $Z_{31}$ | $B(144,17,14536)$ |
| $Z_{32}$ | $B(172,17,20803)$ |

All eight values of $p$ are even, so the complement-symmetry mean of Proposition 1 applies. For each step the verifier recomputes the exact count and extrema tables, forms
$$
\mu=\frac{b^q-1}{4}-(\ell-1),\qquad
\mathrm{gain}=n\sum_{j=1}^{3}\frac1{s_j+\mu}\in\mathbb{Q},
$$
and checks the stored numerator and denominator against the corresponding `Z*_jensen.json` record. The new terminal maximum is the pairing formula $z\leftarrow 3m+5t+z+4$. No floating-point value is used in this arithmetic.

Exact node counts stored in those JSON records (and rechecked by the verifier) are collected in the short appendix. The exact rational gains are large fractions; they are stored as `gain_num` / `gain_den` in the packet and are not restated here.

### 5.2. Truncated bound and increment

Starting from $H\leftarrow \ell_{24}/10^{40}$ and adding the eight exact rational gains, the verifier downward-truncates the resulting rational to $18$ decimal places. The packet, and an independent recomputation on $2026$-$08$-$29$, both print
```
PASS
f(3) >= 3.008500858304852623
```
That $18$-place string is the new truncated lower bound. The increment from the $Z_{24}$ $15$-decimal checkpoint to this truncation is
$$
3.008500858304852623-3.008500791468429=6.6836423623\times 10^{-8},
$$
as printed by `z32_increment.py`.

### 5.3. How the verifier is run

From the unpacked analytic-extension directory,
```
python3 verify_frontier.py
```
The script imports `tables` from `fastscan.py`, reads each `Z*_jensen.json` for $k=25,\dots,32$, and asserts equality of count, block minimum, raw maximum, gain numerator and denominator, and the new $z$. It then asserts that the $18$-place truncation equals `3.008500858304852623`. A run in the preparation of this revision reproduced `PASS` and that string.

The block schedule was obtained by local/IMR-style continuation screening over scale-compatible $(p,q,r)$ families, with limited lookahead. That sentence records a search method. The mathematical claim uses only the eight explicit blocks, Lemmas 1 and 2, Proposition 1, and the $Z_{24}$ integers.

### 5.4. Claim boundary

The packet README states the boundary, which this paper adopts without weakening or expanding:

> The number is a finite-construction certificate, not a proof that this schedule is globally optimal. The block search used local/IMR-style continuation screening and limited lookahead over scale-compatible $(p,q,r)$ families. The mathematical claim still rests on the same two construction lemmas as the $Z_{14}$--$Z_{23}$ chain: AP-freeness of $B(p,q,r)$, and the generalized Wróblewski pairing lemma for arbitrary finite $3$-AP-free blocks.

No result of that screening is imported beyond the eight explicit blocks listed above.

---

## 6. Theorems

**Theorem 1 (April 2026; retained).**
$$f(3)\ge 3.0084928720.$$

*Proof.* By Lemmas 1 and 2 (equivalently, by the convexity proposition and pairing theorem of the April text), $Z_{11}$ is $3$-AP-free. By §3.3, the computed harmonic sum $3.0084928720$ is certified to exceed $3.00849$: the margin $\Delta=2.872\times 10^{-6}$ exceeds the conservative total floating-point error bound $\delta H\le 2.85\times 10^{-21}$ by a factor exceeding $10^{15}$. Two independent implementations with matching node count $1{,}293{,}635{,}551{,}232$ and matching $H$ values provide a cross-check against computational error. Since $Z_{11}$ is a $3$-AP-free subset of $\mathbb{N}^+$ with that certified harmonic sum, $f(3)\ge 3.0084928720$. $\square$

**Theorem 2 (July 2026 certificate).**
Assume Lemmas 1 and 2, the dual-count and integer-truncation assertions of the 24 July 2026 certificate, and the adaptive schedule of §4.2. Then
$$f(3)\ge 3.008500752707.$$

*Proof.* Each block in the schedule is $3$-AP-free by Lemma 1. Each pairing step is $3$-AP-free by Lemma 2. The certificate's per-step integer lower bounds, truncated downward to $12$ decimal places, accumulate to $3.008500752707$ at $Z_{23}$. That number is therefore a harmonic-sum lower bound for an explicit finite $3$-AP-free set. $\square$

**Theorem 3 (August 2026; analytic Jensen extension).**
Assume:

1. Lemma 1: $B(p,q,r)$ is $3$-AP-free for all integers $p\ge 2$, $q\ge 1$, $r\ge 0$;
2. Lemma 2: $\mathrm{pair}(Z,T)$ is $3$-AP-free whenever $Z$ and $T$ are (Wróblewski 1984, Theorem 1);
3. the $Z_{24}$ outward-rounded moment-certificate checkpoint supplies the exact integers
   $$
   \ell_{24}=30085007914684292896021976970576819728782,
   \qquad
   z_{24}=151034783145739756303329159115181,
   $$
   so that the harmonic sum of the finite $3$-AP-free set at that step satisfies
   $$H_{24}\ge \frac{\ell_{24}}{10^{40}},$$
   after pairing with the block $B(64,15,2512)$;
4. the eight subsequent steps use the even-$p$ blocks of §5.1, and the gain of each step is the exact rational of Proposition 1, computed in $\mathbb{Q}$ and stored as the `gain_num` / `gain_den` pair of the corresponding `Z*_jensen.json`.

Then the resulting finite $3$-AP-free set $Z_{32}$ satisfies
$$H(Z_{32})\ge 3.008500858304852623,$$
the $18$-place downward truncation of the exact rational $H_{24}$ plus the eight Jensen gains. Hence
$$f(3)\ge 3.008500858304852623.$$

*Proof.* The set at the $Z_{24}$ checkpoint is $3$-AP-free by Lemmas 1 and 2 applied along the preceding pairing chain, with the checkpoint integers taken as given in hypothesis 3. Each of $Z_{25},\dots,Z_{32}$ is obtained by pairing with a block listed in hypothesis 4, hence is $3$-AP-free by Lemmas 1 and 2. Proposition 1 supplies, for each of those eight steps, an exact rational lower bound on the harmonic increment. Adding those eight rationals to $\ell_{24}/10^{40}$ produces an exact rational lower bound on $H(Z_{32})$. Downward truncation to $18$ decimal places yields $3.008500858304852623$, as asserted by `verify_frontier.py` and independently reprinted on $2026$-$08$-$29$. $\square$

Theorem 3 is conditional on the $Z_{24}$ integers. Those integers are not derived in this file from a $Z_{24}$ JSON artefact, because no such file is in the packet. The claim remains a finite explicit construction bound.

---

## 7. Reproducibility

**April 2026 computational spine.** The public repository https://github.com/darbybailey/f3-lower-bound/ contains `f3_lower_bound.py`, the verified $Z_{13}$ checkpoint `z_state.json`, `README.md`, `SHA256SUMS.txt`, `CITATION.cff`, and `LICENSE`. Requirements: Python $3.8$ or later (standard library only) and GCC $\ge 4.6$ or Clang $\ge 3.1$. OpenMP is optional; the two implementations agree at the stated $Z_{12}$ checkpoint. Running `python f3_lower_bound.py` compiles the C kernel, reconstructs the $Z_{11}$ shift chain, and resumes from `z_state.json` when that file is present. The canonical v1 deposit is Zenodo $10.5281/\mathrm{zenodo}.19722202$ (concept $10.5281/\mathrm{zenodo}.19722201$). This revision is reserved as version DOI $10.5281/\mathrm{zenodo}.24486265$.

**July 2026 certificate.** The file `F3_CERTIFICATE.txt` (24 July 2026) is the integer record of the normalized chain through $Z_{13}$ and the adaptive chain $Z_{14}$--$Z_{23}$. Its own reproduction line is `python3 f3_certify.py all`. The $12$-decimal lowers in §4.2 are taken from that file.

**August 2026 analytic-extension packet.** The packet `f3_lower_frontier_Z32_2026-08-29.zip` has SHA-256
```
120aa0b44fb90065364ec5c533265b1ad9c19c71af8f973396c788726f1d1a86
```
It contains `verify_frontier.py`, `fastscan.py`, `cert_jensen_one.py`, `frontier_summary.json`, `Z25_jensen.json` through `Z32_jensen.json`, and a README. It does not contain a $Z_{24}$ JSON file; `verify_frontier.py` hard-codes $\ell_{24}$ and $z_{24}$. To recheck Theorem 3:
```
python3 verify_frontier.py
```
Expected output:
```
PASS
f(3) >= 3.008500858304852623
```
A supporting increment check, `z32_increment.py`, prints $Z_{32}-Z_{24}=6.6836423623\times 10^{-8}$ from the published truncations. The block-identity and pairing-interval audits used in §2 live with the same date stamp and printed `PASS` on the recorded runs (`T_p` identity on $20136$ triples; pairing hull, zero fail flags).

This revision has reserved version DOI $10.5281/\mathrm{zenodo}.24486265$, issued as a new version of $10.5281/\mathrm{zenodo}.19722202$ under the concept DOI $10.5281/\mathrm{zenodo}.19722201$. The analytic-extension packet is attached to that record. The public GitHub default branch has not yet been updated to include the $Z_{14}$--$Z_{32}$ materials.

---

## 8. Prior work

Wróblewski (1984) established $f(3)\ge 3.00849$ using the identical pairing construction on an ODRA 1305 mainframe. The April 2026 version of the present paper extended the $p=6$ computation through $Z_{13}$, added the overflow safeguards required at $Z_{12}$ and beyond, and certified $f(3)\ge 3.0084928720$ at $Z_{11}$. That version is the canonical public record: DOI $10.5281/\mathrm{zenodo}.19722202$.

Walker (2022, rev. 2025) conducted an exhaustive branch-and-bound search over Kempner sets in bases $b\le 120$, finding no improvement to $f(3)$. Wróblewski's construction uses Behrend sphere blocks with level-set constraints on digit sums; Kempner sets restrict digit membership. The two approaches are structurally distinct.

On the upper-bound side, there has been substantial recent progress on $3$-term progression-free sets. Bloom and Sisask (2020) proved quantitative density bounds on $3$-AP-free sets that imply $f(3)<\infty$. Kelley and Meka (2023) obtained substantially stronger quantitative bounds. Bloom and Sisask (2023) gave an exposition and refinement of the Kelley–Meka method. Raghavan (2026) further improved the exponent to $1/6$ up to a $\log\log N$ factor. These works concern upper bounds for $r_3(N)$, not explicit numerical upper bounds for the harmonic sum $f(3)$. To the author's knowledge, no explicit numerical upper bound for $f(3)$ has yet been published.

The best explicit lower bound prior to the April 2026 computation remained Wróblewski's construction. The April theorem improved the printed digits to $3.0084928720$. The present revision records a longer explicit chain and a larger truncated lower bound, still of the same pairing type, and still a finite-construction statement.

Behrend (1946) is the classical source of the sphere-in-the-box construction of which the blocks $B(p,q,r)$ are a discrete form. Higham (2002) is the reference for the $Z_{11}$ summation-error bound.

---

## 9. AI collaboration and validation disclosure

The April 2026 text was developed with the assistance of AI systems from OpenAI, Anthropic, Google, and xAI during March–April 2026, used for editorial feedback, code review discussion, literature cross-checking, and reproducibility review. All long production computations, local reruns, smoke tests, checkpoint validation, and final judgment about the April claims were performed by the author in a local Apple MacBook Pro terminal environment.

The August 2026 analytic extension used Grok (xAI) and Gemini (Google) as advisory tools: for discussion of the even-$p$ Jensen increment, for interval-hull and $T_p$-identity audit scripts, and for editorial comments on claim boundaries. Those systems did not author the lemmas, did not supply unsourced numerical claims, and were not treated as a proof assistant of record. All AI-assisted suggestions in both periods were treated as advisory rather than authoritative.

The author takes full responsibility for the correctness, wording, and scope of all claims. The recursive pairing construction and the proof of $3$-AP-freeness of the pairing are due to Wróblewski (1984). The April computational implementation, overflow safeguards, floating-point error certificate, and the $Z_{11}$ theorem are the contribution of the first version of this paper. The written general forms of Lemmas 1 and 2, Proposition 1, the recording of the July adaptive certificate, the $Z_{24}$ integer checkpoint, and the exact-rational extension through $Z_{32}$ are the contribution of this revision.

Author note. This line of inquiry grew out of the author's prior work in audio signal processing, communications systems, and cybernetics, where questions of structure, resonance, and computational constraint motivated the present investigation.

---

## References

1. Wróblewski, J. (1984). A nonaveraging set of integers with a large sum of reciprocals. *Mathematics of Computation*, 43(167), 261–262. https://doi.org/10.2307/2007580
2. Bloom, T., & Sisask, O. (2020). Breaking the logarithmic barrier in Roth's theorem on arithmetic progressions. arXiv:2007.03528. https://arxiv.org/abs/2007.03528
3. Kelley, Z., & Meka, R. (2023). Strong bounds for 3-progressions. FOCS 2023 Best Paper. arXiv:2302.05537. https://arxiv.org/abs/2302.05537
4. Bloom, T. F., & Sisask, O. (2023). The Kelley–Meka bounds for sets free of three-term arithmetic progressions. *Essential Number Theory*, 2(1), 15–44. https://doi.org/10.2140/ent.2023.2.15 (arXiv:2302.07211)
5. Raghavan, R. (2026). Improved bounds for 3-progressions. arXiv:2603.27045. https://arxiv.org/abs/2603.27045
6. Walker, A. (2022, rev. 2025). Integer sets of large harmonic sum which avoid long arithmetic progressions. arXiv:2203.06045. https://arxiv.org/abs/2203.06045
7. Behrend, F. A. (1946). On sets of integers which contain no three terms in arithmetical progression. *Proceedings of the National Academy of Sciences*, 32(12), 331–332. https://doi.org/10.1073/pnas.32.12.331
8. Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms* (2nd ed.). SIAM. https://doi.org/10.1137/1.9780898718027

---

## Appendix. Stored $Z_{25}$--$Z_{32}$ counts

Exact node counts from the packet JSON files, rechecked by `verify_frontier.py`. These integers are not used as floating-point data.

| Step | Block | Count |
|---|---|---|
| $Z_{25}$ | $B(80,15,3910)$ | $15{,}062{,}258{,}981{,}484{,}173{,}320{,}028{,}160$ |
| $Z_{26}$ | $B(98,15,5878)$ | $210{,}685{,}685{,}278{,}136{,}881{,}913{,}856{,}000$ |
| $Z_{27}$ | $B(122,15,9151)$ | $3{,}633{,}095{,}646{,}662{,}660{,}213{,}063{,}024{,}640$ |
| $Z_{28}$ | $B(108,16,7633)$ | $77{,}945{,}839{,}849{,}627{,}820{,}343{,}268{,}737{,}024$ |
| $Z_{29}$ | $B(132,16,11464)$ | $1{,}293{,}410{,}754{,}478{,}786{,}391{,}940{,}445{,}634{,}560$ |
| $Z_{30}$ | $B(164,16,17779)$ | $26{,}990{,}593{,}877{,}393{,}764{,}352{,}384{,}963{,}706{,}880$ |
| $Z_{31}$ | $B(144,17,14536)$ | $610{,}989{,}486{,}686{,}452{,}172{,}948{,}206{,}096{,}744{,}448$ |
| $Z_{32}$ | $B(172,17,20803)$ | $8{,}775{,}434{,}006{,}765{,}867{,}312{,}757{,}618{,}212{,}864{,}000$ |

---

Correspondence: ORCID 0000-0003-4072-6191\
This version DOI: 10.5281/zenodo.24486265\
Canonical v1 DOI: 10.5281/zenodo.19722202 (concept 10.5281/zenodo.19722201)\
© 2026 Darby Bailey McDonough. CC BY 4.0.

*v2, August 2026.*
