# CHANGELOG — f(3) preprint v2 vs v1

**Status.** Reserved unpublished new version of Zenodo 10.5281/zenodo.19722202, under the existing concept DOI 10.5281/zenodo.19722201. Version DOI: 10.5281/zenodo.24486265. Do not overwrite v1. Do not open a new concept DOI.

**v1 (April 2026; Zenodo 10.5281/zenodo.19722202; concept 10.5281/zenodo.19722201).**
Computational theorem $f(3)\ge 3.0084928720$ at $Z_{11}$, by Wróblewski pairing of Behrend blocks $B(6,q,r)$ (and the earlier $p=4$ steps) on the greedy seed $G_3$, with a C+Python enumeration kernel and an IEEE 754 error certificate. Supporting streamed values $H(Z_{12})=3.0084948386$ and $H(Z_{13})=3.0084958842$. Construction in v1 ran $p=6$ through $Z_{13}$.

**What this v2 revision adds (August 2026).**

1. **Expanded construction lemmas** (paper §2). A written Block AP-free lemma for general $B(p,q,r)$, with the carry-free base-$2p-1$ argument and the second-difference identity $T_p(x)+T_p(z)-2T_p(y)=(x-z)^2/4$. A written Pairing lemma that records the implemented interval form of Wróblewski (1984, Theorem 1) and treats the interval-hull audit as supporting verification of those shifts, not as a replacement for the classical case analysis. An even-$p$ Jensen proposition with exact rational arithmetic.

2. **Adaptive interval-certified chain $Z_{14}$–$Z_{23}$** (paper §4), from the 24 July 2026 certificate. Compact schedule (block, dual-verified count, $12$-decimal cumulative lower). Certificate Theorem 2: $f(3)\ge 3.008500752707$. Certificate Theorem 1 ($f(3)\ge 3.0084975903$, three Jensen steps from the published $Z_{11}$ value) is recorded as a compact intermediate, not as the headline.

3. **$Z_{24}$ checkpoint** (paper §4). Outward-rounded moment-certificate integers consumed by the verifier: block $B(64,15,2512)$, $z=151034783145739756303329159115181$, $\ell_{24}/10^{40}$ with $\ell_{24}=30085007914684292896021976970576819728782$, and the $15$-decimal machine-interval lower $H\ge 3.008500791468429$. **There is no $Z_{24}$ JSON file** and no `Z24_jensen.json`. The paper says so.

4. **Analytic Jensen extension $Z_{25}$–$Z_{32}$** (paper §5). Eight explicit even-$p$ blocks; exact rational gains; no float in the theorem arithmetic. New truncated bound $f(3)\ge 3.008500858304852623$. Independent `python3 verify_frontier.py` run: `PASS` and that $18$-place truncation. Packet zip SHA-256 `120aa0b44fb90065364ec5c533265b1ad9c19c71af8f973396c788726f1d1a86`. Increment $Z_{32}-Z_{24}=6.6836423623\times 10^{-8}$.

5. **Theorems restated** (paper §6). v1 theorem retained. New theorem for $3.008500858304852623$ with explicit hypotheses: the two lemmas, the $Z_{24}$ checkpoint integers, and the exact Jensen steps.

6. **Reproducibility** (paper §7). Points at the v1 GitHub tree and at the analytic-extension packet; records how to run the verifier; records the zip checksum; states that $Z_{24}$ is a hardcoded integer checkpoint.

7. **Claim boundary** (quoted from the packet README). Finite explicit construction. Not a proof that the schedule is globally optimal. Block search mentioned only as local/IMR-style continuation screening. No result of that screening is imported beyond the explicit blocks.

8. **AI disclosure updated** (paper §9) for the August 2026 analytic extension. Grok and Gemini used as advisory tools; the author remains responsible.

**What v2 does not do.**
It does not replace or withdraw the April $Z_{11}$ computational certificate. It does not claim a $Z_{24}$ JSON artefact that is not in the packet. It does not invent per-step floating-point $H$ values beyond those printed in v1, in `F3_CERTIFICATE.txt`, in `frontier_summary.json`, in the $Z_{25}$–$Z_{32}$ JSON files, or in the verifier output. It is a new version of the published v1 record, not a second concept DOI.

**Canonical v1 DOI (unchanged).** 10.5281/zenodo.19722202.

**This version DOI.** 10.5281/zenodo.24486265 (new version of 10.5281/zenodo.19722202; concept 10.5281/zenodo.19722201).

**Upload pack for reserved record 24486265.** Typeset PDF, this markdown, this CHANGELOG, and `f3_lower_frontier_Z32_2026-08-29.zip` (SHA-256 `120aa0b44fb90065364ec5c533265b1ad9c19c71af8f973396c788726f1d1a86`). Headline Theorem 1 remains 3.0084928720; the Z32 number is Theorem 3 only, with the Z24 hypotheses. Recommended additional attachments (not in the Z32 zip; needed to inspect Theorem 2's July certificate): `F3_CERTIFICATE.txt` and `f3_certificate.json`. A later public GitHub update to `darbybailey/f3-lower-bound` should add the v2 packet, checksums, and citation DOI, and should also attach the source, checkpoint, and checksums v1 already promised so the record is no longer PDF-only.
