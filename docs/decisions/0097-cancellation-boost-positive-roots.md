# ADR-0097: input-encoded cancellation at ln/asin/lgamma/digamma — three fixes, one mechanism

- **Status**: accepted
- **Date**: 2026-06-11

## Context

The 2026-06-10 workspace review (epic pf-8iji; pf-smcb, pf-rylv, pf-wmv7)
confirmed the largest Theme-1 cluster: the Ziv driver's first-iteration error
model — relative half-width `|y|·2^-(w−guard)` — is violated by kernels whose
evaluation cancels against structure carried in the input, and the interval
test then certifies a wrong value. The RC2 fix (`cancellation_boosted`,
review 2026-05-29) covered the negative-axis reflection branches of
lgamma/digamma; these paths had no protection:

- **ln just below 1** (pf-smcb): `ln(1 − 2^-80)` at p100 → 53 returned a
  value 2^15 ulps wrong, certified INEXACT. The range reduction
  `ln(x) = ln(m) + e·ln2` at `e = −1` computes `ln(2x) − ln2`, a
  cancellation whose depth is the input's proximity to 1.
- **asin near ±1** (pf-rylv): the `x²` rounding error is amplified `2^d`
  through `1 − x²` (`d` = the proximity depth), then propagated through
  sqrt and atan. The review's named reproducer `1 − 2^-200` at p400 was a
  **misadjudication** — a single-bit delta yields the sparse
  `x² = 1 − 2^-199 + 2^-400`, exactly representable at the first working
  precision, and the recorded review output was limb-exact against mpmath
  RN400 (the reviewer compared Display-truncated digits; the same
  reviewer-oracle failure class as J_100(200) and the beta truth value).
  The mechanism is real with a **dense** delta: `x = 1 − RN150(π)·2^-202`
  (exact at p400) was ~2^18 ulps@400 wrong, verified by run pre-fix.
- **lgamma/digamma positive roots** (pf-wmv7): `lgamma(2 + 2^-100)` at
  p120 → 53 had relative error 2.5e-3; digamma at the p100 rounding of its
  positive root 1.46163… certified garbage. The Stirling-with-shift (and
  Spouge) compositions cancel near the roots at 1, 2, and 1.46163…, with
  the additional subtlety that `z_min` sized from the *target* precision
  caps the truncation accuracy regardless of working precision, so no
  amount of plain Ziv guard growth can recover.

## Decision

Three fixes matched to where the cancellation lives, most-structural first.

**ln: remove the cancellation (structural), plus a proximity boost.** For
`x ∈ [0.5, 1)` (`e = −1`) the atanh-series argument
`u = (x − 1)/(x + 1) ∈ [−1/3, 0)` is already inside the series' convergent
range, so the reduction is skipped entirely: `m = x`, `e = 0`. `x − 1` is
Sterbenz-exact at the working precision, so the direct path carries full
relative accuracy relative to `x_w` with no extra precision and no second
evaluation; it also drops the `ln2` computation. Only `e = −1` manufactures
cancellation (`e = 0` is already direct; every other `e` keeps the result
`≥ ln2·|e| − ln2` away from 0). The structural fix alone left the
input-rounding gap, caught by the slice's adversarial verification: when
the proximity depth exceeds `w`, `round_to_precision(x, w)` collapses
`x_w` to exactly 1, the series returns exact 0, and `half_width(0) = 0`
certifies +0 on the first iteration (`ln(1 − 2^-200 @p400) → 53` returned
+0 INEXACT; ADR-0050's cancellation-to-zero class). The asin pattern
closes it: `x − 1` computed exactly at the input precision, and every Ziv
working precision boosted past the depth when `exponent(x − 1) ≤ −2`, on
both sides of 1.

**asin: deterministic pre-boost.** The amplification depth is input-encoded
and computable exactly up front: `gap = 1 − |x|` (Sterbenz-exact at the
input precision). When `exponent(gap) ≤ −2`, every Ziv evaluation runs at
`w + boost` with `boost = (−exponent(gap) − 1) + 8 ≈ d + 7`, where
`d = −exponent(1 − x²)`. Derivation: the `x²` rounding error `~2^-(w+boost)`
becomes a relative error `2^(d − w − boost)` of `1 − x²`, halves through
sqrt, and passes through the atan composition with sensitivity ≤ 1, so the
result's relative error is `≤ 2^(d − w − boost + c) ≤ 2^-(w + 5)` for small
`c` — inside the half-width model with margin. No probing re-evaluation
(unlike `cancellation_boosted`): the boost is exact, bounded by the input
precision, and zero for inputs away from ±1.

**lgamma/digamma: conditional realised-cancellation boost.** Inside the
positive-root windows — `[3/4, 5/4] ∪ [7/4, 9/4]` for lgamma (roots 1, 2),
`[5/4, 7/4]` for digamma (root 1.46163…) — the positive branch routes
through `cancellation_boosted`, mirroring the negative-branch RC2 fix. The
positive-branch evaluation is extracted to return `(value, operand_scale)`
(Stirling-shift: max exponent of `lgamma(z)` and `ln(product)`; Spouge:
max exponent of the `(z+1/2)ln(z+a)`, `(z+a)`, `lnΓ(z+1)`, `ln z` chain;
direct Stirling: the value's own exponent — no cancellation). Three
load-bearing details: (1) `z_min` is re-derived from the *boosted*
precision inside the closure (`z_min_for_target(w)`), because the
truncation floor of the shifted asymptotic is set by `z_min`, not by the
working precision; (2) the windows are compared on the *original* input
with exact dyadic quarter constants, so the trigger is
precision-independent; (3) Spouge's reported scale must also charge the
sum `S(z, a)`'s internal alternating cancellation (its largest term
exceeds `S` by ~0.1·w bits, hidden behind the `ln`) — without it the
boost stops short and `lgamma(2 + 2^-500 @p520) → 400` certified a value
~2^21 ulps wrong (found by the slice's adversarial verification;
the same under-report erodes accuracy on every outside-window Spouge
call too, filed as pf-dzty). Outside the windows `|lgamma| ≥ ~2^-3.4` and
`|ψ| ≥ ~2^-2.1`, so the realised cancellation stays inside the 24-bit Ziv
guard for the differential-lane precision range; the un-boosted path pays
nothing, and the in-window cost is one extra evaluation (the
`cancellation_boosted` probe), confined to two half-unit windows.

## Consequences

- The regression lane pins all four kernels bit-exactly against mpmath at
  exact dyadic inputs: ln at depth 80 (plus the above-1 control), asin with
  the dense-delta reproducer (plus the sparse control documenting the
  misadjudication), lgamma at both roots, digamma at the p100-parse of its
  root (plus a shallow control). gamma/beta/zeta-FE compositions inherit
  the lgamma fix transitively.
- ln in [0.5, 1) gets cheaper (no ln2), not just correct. asin pays the
  boost only within `2^-2` of ±1. lgamma/digamma pay double evaluation
  only inside the root windows.
- The corrected pf-rylv reproducer is recorded here and in the lane: the
  review's named input was never broken (sparse-delta exactness), the
  dense-delta band was. Verify-the-verdict held again.
- Failure modes considered (inverted): (1) ln's direct path at the
  convergence edge `x = 0.5` gives `u = −1/3` exactly, the documented
  series bound; correctness there is guarded by the seam probes — and the
  first draft of this very fix HAD the input-collapse failure (proximity
  deeper than `w`), found by the adversarial verification and closed with
  the proximity boost; the lane pins both sides of 1 at depth 200.
  (2) asin's boost could under-cover if the amplification analysis missed
  a factor; the adversarial verification re-derived it (boost `d + 8`
  bounds both the input-rounding and the promoted `x²` error under
  `2^-(w+8)` for arbitrary depth) and probed deeper and denser deltas than
  the lane, all bit-exact. (3) The lgamma boost's operand scale could hide
  internal cancellation — it DID, in Spouge's sum, the first draft's
  second failure; charged now, with the deep-root reproducer pinned in
  the lane. (4) The lgamma/digamma windows could be too narrow at high
  Stirling targets: measured +2..+4-bit guard-model violations just
  outside the windows at w = 600 (values still correct; pre-existing,
  trace-identical at baseline), filed as pf-dzty with pf-2thy
  cross-linked rather than silently widened here.
- Adjacent finding (not this slice): acos near 1 has the same
  input-collapse family with a worse posture — `+0 with Status OK` on a
  transcendental result (the INEXACT force fires only on `Class::Normal`).
  Pre-existing, verified at baseline; filed as pf-9761
  (certified-wrong-with-OK family).
- In-window cost measured by the verifier: 2.2–3.2× per call against the
  un-boosted self (the `cancellation_boosted` probe evaluation plus the
  boosted-precision `z_min`), confined to the half-unit windows;
  outside-window and ln/asin costs are noise-level.

## Related

- Issues: pf-smcb, pf-rylv, pf-wmv7 (closed by this ADR), epic pf-8iji;
  pf-2thy (adjacent probe, the outside-window marginal band note).
- Review: `~/.claude/plans/pfloat-workspace-review-2026-06-10.md` Theme 1
  items 1–3; reproducer checks F1/F2(+F2b)/F3 in
  `~/.claude/plans/pfverify-harness/`.
- Other ADRs: ADR-0050 (tanh cancellation-to-zero), the RC2 review fix
  (`cancellation_boosted`), ADR-0095/0096 (this arc), ADR-0063/0064
  (INEXACT posture).
