# f(3) Lower Bound — Improved Computational Result

**Version 2.0 (published 2026-08-29): f(3) ≥ 3.008500858304852623**

This repository contains the implementation and certificates for:

> McDonough, D. B. (2026). *An Improved Computational Lower Bound for f(3)
> via Wróblewski's Recursive Pairing of Behrend Sphere Blocks.* Preprint,
> Version 2.0, published 2026-08-29.
> Version DOI: [10.5281/zenodo.24486265](https://doi.org/10.5281/zenodo.24486265).
> Concept DOI: [10.5281/zenodo.19722201](https://doi.org/10.5281/zenodo.19722201).

Version 2.0 retains the April 2026 computational theorem and records two later
extensions of the same Wróblewski pairing of Behrend sphere blocks:

| Version | Bound | Basis | DOI |
|---|---|---|---|
| **2.0** (2026-08-29) | **f(3) ≥ 3.008500858304852623** | Exact-rational Jensen steps Z₂₅–Z₃₂, conditional on the Z₂₄ checkpoint integers | [10.5281/zenodo.24486265](https://doi.org/10.5281/zenodo.24486265) |
| 2.0, Theorem 2 | f(3) ≥ 3.008500752707 | Interval-certified adaptive chain Z₁₄–Z₂₃ (24 July 2026 certificate) | same Version 2.0 record |
| **1.0** (2026-04-24) | **f(3) ≥ 3.0084928720** | Certified C+Python enumeration at Z₁₁, with IEEE 754 error certificate | [10.5281/zenodo.19722202](https://doi.org/10.5281/zenodo.19722202) |

The Version 1.0 computational spine in this repository is unchanged. The Version 2.0
paper, changelog, July certificate, and Z₃₂ packet are added alongside it.
The typeset Version 2.0 PDF and the Z₃₂ zip are on the Zenodo record; checksums
are in `SHA256SUMS.txt`. Those two binaries are not duplicated in this git tree.

This result improves the best published lower bound due to Wróblewski (1984).

---

## Repository Layout

```
README.md           This file
LICENSE             MIT License
CITATION.cff        Citation metadata (Version 2.0)
f3_lower_bound.py   v1 main script — construction, verification, and extension
z_state.json        v1 verified Z₁₃ checkpoint (resume without recomputing Z₁₂)
DEVNOTES.md         Engineering history, overflow bugs found and fixed, GCP notes
SHA256SUMS.txt      v1 spine checksums plus all Zenodo Version 2.0 file SHA-256s
paper/Improved_Computational_Lower_Bound_f3.pdf
                    v1 typeset preprint (unchanged)
paper/Improved_Computational_Lower_Bound_f3_v2.md
                    Version 2.0 manuscript source
f3_v2_CHANGELOG.md  Version 2.0 vs 1.0 changelog (byte-identical to the Zenodo file)
certificate/        24 July 2026 interval certificate (Theorem 2)
  F3_CERTIFICATE.txt
  f3_certificate.json
z32_packet/         August 2026 analytic-extension packet (from the Zenodo zip)
  verify_frontier.py
  fastscan.py
  cert_jensen_one.py
  Z25_jensen.json … Z32_jensen.json
  frontier_summary.json
  README.md
  SHA256SUMS.txt
```

The changelog file still describes the Version 2.0 DOI as reserved, because that
is the text deposited on Zenodo. The record is now published: Version DOI
[10.5281/zenodo.24486265](https://zenodo.org/records/24486265), concept DOI
[10.5281/zenodo.19722201](https://zenodo.org/doi/10.5281/zenodo.19722201).

---

## Version 2.0: analytic Jensen extension (Z₂₅–Z₃₂)

Starting from the outward-rounded Z₂₄ moment-certificate integers (block
B(64,15,2512); there is no Z₂₄ JSON file), eight explicit even-p sphere blocks
are appended. Complement symmetry gives the exact mean digit value; Jensen's
inequality gives an exact rational lower bound for all three translated copies.
No floating-point value is used in that theorem arithmetic. The downward-truncated
bound is

    f(3) >= 3.008500858304852623.

Verify with:

```bash
cd z32_packet
python3 verify_frontier.py
```

A successful run prints `PASS` and `f(3) >= 3.008500858304852623` (standard
library only). Packet file checksums are in `z32_packet/SHA256SUMS.txt`. The
parent zip on Zenodo has SHA-256

    120aa0b44fb90065364ec5c533265b1ad9c19c71af8f973396c788726f1d1a86  f3_lower_frontier_Z32_2026-08-29.zip

The number is a finite-construction certificate. It is not a proof that the
block schedule is globally optimal.

---

## What is f(3)?

`f(3)` is the supremum of harmonic sums `Σ 1/n` over all 3-AP-free subsets of the
positive integers. Bloom and Sisask (2020) proved `f(3)` is finite; the explicit
value is unknown. The construction here is Wróblewski's (1984) recursive pairing of
Behrend sphere blocks — the AP-free property follows from the pairing theorem by
induction and requires no computational verification.

---

## Requirements (Version 1.0 computational spine)

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

The Version 2.0 frontier verifier needs only Python 3 (standard library).

---

## Running the Version 1.0 enumerator

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

## Version 1.0 results (unchanged)

| Step | H | Gain | Nodes |
|---|---|---|---|
| Z₁₁ ★ | **3.0084928720** | +0.0000036987 | 1,293,635,551,232 |
| Z₁₂ | 3.0084948386 | +0.0000019666 | 7,563,706,368,000 |
| Z₁₃ | 3.0084958842 | +0.0000010455 | 44,239,186,034,688 |

The Version 1.0 lower bound **f(3) ≥ 3.0084928720** is established at Z₁₁.
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

The Version 1.0 margin above Wróblewski's bound (Δ = 2.872×10⁻⁶) exceeds the conservative
IEEE 754 floating-point error bound (δ_H ≤ 2.85×10⁻²¹) by a factor exceeding 10¹⁵.
See Section 3.2 of the Version 1.0 paper for the full derivation separating denominator rounding
and summation error (Higham, 2002, Thm. 4.1).

---

## Checksums

See `SHA256SUMS.txt` for the Version 1.0 spine files and for every file in the
published Zenodo Version 2.0 record, including
`f3_lower_frontier_Z32_2026-08-29.zip` and
`Improved_Computational_Lower_Bound_f3_v2.pdf`.

Verify the v1 spine with: `sha256sum f3_lower_bound.py z_state.json DEVNOTES.md`

Verify the unpacked Z₃₂ packet with: `sha256sum -c z32_packet/SHA256SUMS.txt`

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

To cite Version 2.0: see `CITATION.cff` or use the Version DOI
[10.5281/zenodo.24486265](https://doi.org/10.5281/zenodo.24486265).
The concept DOI [10.5281/zenodo.19722201](https://doi.org/10.5281/zenodo.19722201)
always resolves to the latest version. Version 1.0 remains citable at
[10.5281/zenodo.19722202](https://doi.org/10.5281/zenodo.19722202).

*Pip Projects, Inc. · Digital Sidewalk Lab · Independent Researcher*  
ORCID: [0000-0003-4072-6191](https://orcid.org/0000-0003-4072-6191)
