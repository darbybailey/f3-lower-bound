# f(3) Lower Bound — Improved Computational Result

**f(3) ≥ 3.0084928720**

This repository contains the implementation for the paper:

> McDonough, D. B. (2026). *An Improved Computational Lower Bound for f(3)
> via Wróblewski's Recursive Pairing of Behrend Sphere Blocks.* Preprint.
> DOI: 10.5281/zenodo.19722202

This result improves the best published lower bound due to Wróblewski (1984).

---

## Repository Layout

```
README.md           This file
LICENSE             MIT License
CITATION.cff        Citation metadata
f3_lower_bound.py   Main script — construction, verification, and extension
z_state.json        Verified Z₁₃ checkpoint (resume without recomputing Z₁₂)
DEVNOTES.md         Engineering history, overflow bugs found and fixed, GCP notes
SHA256SUMS.txt      File checksums
```

---

## What is f(3)?

`f(3)` is the supremum of harmonic sums `Σ 1/n` over all 3-AP-free subsets of the
positive integers. Bloom and Sisask (2020) proved `f(3)` is finite; the explicit
value is unknown. The construction here is Wróblewski's (1984) recursive pairing of
Behrend sphere blocks — the AP-free property follows from the pairing theorem by
induction and requires no computational verification.

---

## Requirements

- Python 3.8+ (standard library only)
- GCC ≥ 4.6 **or** Clang ≥ 3.1 (required for `unsigned __int128`)
- OpenMP (optional, but recommended for multi-core speedup)

**macOS:**
```bash
brew install libomp
```

**Linux / GCP:**
```bash
sudo apt-get install -y gcc libomp-dev
```

---

## Running

```bash
# Fresh run — builds from Z₁₁, verifies Z₁₂, continues
python f3_lower_bound.py

# Resume from Z₁₃ checkpoint (z_state.json must be in working directory)
python f3_lower_bound.py
```

The script compiles the C kernel automatically on first run.

---

## Verification Checkpoint

At Z₁₂ the script checks two independent implementations agree before continuing:

```
Z_12:  7,563,706,368,000 nodes  2.619h   (12-core OpenMP, Apple M3)
       H    = 3.0084948386
       gain = +0.0000019666
✅ Z12 VERIFIED: H and node count match ground truth
```

If the node count or H value differs, stop. See `DEVNOTES.md` for the overflow bugs
that caused earlier incorrect results.

---

## Results

| Step | H | Gain | Nodes |
|---|---|---|---|
| Z₁₁ ★ | **3.0084928720** | +0.0000036987 | 1,293,635,551,232 |
| Z₁₂ | 3.0084948386 | +0.0000019666 | 7,563,706,368,000 |
| Z₁₃ | 3.0084958842 | +0.0000010455 | 44,239,186,034,688 |

The lower bound **f(3) ≥ 3.0084928720** is established at Z₁₁.
Z₁₂ and Z₁₃ are supporting data confirming gain ratio stability at 0.5317.

---

## Approximate Wall Times (Apple M3, 12-core OpenMP)

| Step | Nodes | Time |
|---|---|---|
| Z₁₁ | 1.29 × 10¹² | 6.0 h (single-thread) |
| Z₁₂ | 7.56 × 10¹² | 2.619 h |
| Z₁₃ | 4.42 × 10¹³ | 15.1 h |
| Z₁₄ | ~2.56 × 10¹⁴ | ~87 h |

Each step takes approximately 5.8× longer than the previous.

---

## Numerical Certification

The margin above Wróblewski's bound (Δ = 2.872×10⁻⁶) exceeds the conservative
IEEE 754 floating-point error bound (δ_H ≤ 2.85×10⁻²¹) by a factor exceeding 10¹⁵.
See Section 3.2 of the paper for the full derivation separating denominator rounding
and summation error (Higham, 2002, Thm. 4.1).

---

## Checksums

```text
392e31b4e0ccecc29969f33585519341ffb59c64fe8ac964dc9814eec66a858b  f3_lower_bound.py
cb371461c83d2a07d6eff58d778c27846dda5f5bda6ac1d8060b80bb08d97da4  z_state.json
4dd451d1138fc3ee45e35b5355a4f15fb0162a5723c56097b052dd1c8a396828  DEVNOTES.md
```

Verify with: `sha256sum f3_lower_bound.py z_state.json DEVNOTES.md`

---

## References

1.	Wróblewski, J. (1984). A nonaveraging set of integers with a large sum of reciprocals. Mathematics of Computation, 43(167), 261–262. https://doi.org/10.2307/2007580 
2.	Bloom, T., & Sisask, O. (2020). Breaking the logarithmic barrier in Roth's theorem on arithmetic progressions. arXiv:2007.03528. https://arxiv.org/abs/2007.03528 
3.	Kelley, Z., & Meka, R. (2023). Strong bounds for 3-progressions. FOCS 2023 Best Paper. arXiv:2302.05537. https://arxiv.org/abs/2302.05537 
4.	Bloom, T. F., & Sisask, O. (2023). The Kelley–Meka bounds for sets free of three-term arithmetic progressions. Essential Number Theory, 2(1), 15–44. https://doi.org/10.2140/ent.2023.2.15 (arXiv:2302.07211) 
5.	Raghavan, R. (2026). Improved bounds for 3-progressions. arXiv:2603.27045. https://arxiv.org/abs/2603.27045 
6.	Walker, A. (2022, rev. 2025). Integer sets of large harmonic sum which avoid long arithmetic progressions. arXiv:2203.06045. https://arxiv.org/abs/2203.06045 
7.	Behrend, F. A. (1946). On sets of integers which contain no three terms in arithmetical progression. Proceedings of the National Academy of Sciences, 32(12), 331–332. https://doi.org/10.1073/pnas.32.12.331 
8.	Higham, N. J. (2002). Accuracy and Stability of Numerical Algorithms (2nd ed.). SIAM. https://doi.org/10.1137/1.9780898718027

---

## License and Citation

Code: MIT License — see `LICENSE`.  
Paper: CC BY 4.0.

To cite: see `CITATION.cff` or use the Zenodo DOI: [10.5281/zenodo.19722202](https://doi.org/10.5281/zenodo.19722202)

*Pip Projects, Inc. · Digital Sidewalk Lab · Independent Researcher*  
ORCID: [0000-0003-4072-6191](https://orcid.org/0000-0003-4072-6191)
