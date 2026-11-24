# f(3) lower-bound frontier — 2026-08-29

This packet extends the normalized adaptive Wróblewski/Behrend-block chain.

## Result

Starting from the exact outward-rounded Z24 checkpoint

    H >= 3.008500791468429...

with terminal state

    z = 151034783145739756303329159115181,

we append eight explicit even-p sphere blocks using the generalized pairing construction:

    Z25 B(80,15,3910)
    Z26 B(98,15,5878)
    Z27 B(122,15,9151)
    Z28 B(108,16,7633)
    Z29 B(132,16,11464)
    Z30 B(164,16,17779)
    Z31 B(144,17,14536)
    Z32 B(172,17,20803)

For each appended block, complement symmetry gives the exact mean digit value and Jensen's inequality gives an exact rational lower bound for all three translated copies. No floating-point value is used in the theorem arithmetic.

The resulting finite explicit construction gives the downward-truncated bound

    f(3) >= 3.008500858304852623.

## Claim boundary

The number is a finite-construction certificate, not a proof that this schedule is globally optimal. The block search used local/IMR-style continuation screening and limited lookahead over scale-compatible (p,q,r) families. The mathematical claim still rests on the same two construction lemmas as the Z14-Z23 chain: AP-freeness of B(p,q,r), and the generalized Wróblewski pairing lemma for arbitrary finite 3-AP-free blocks. Those lemmas should be stated and proved at full generality in the paper before promotion as a standalone theorem.

`verify_frontier.py` recomputes the exact count/extrema tables and exact rational Jensen gains for Z25-Z32 and checks the final 18-place truncation.
