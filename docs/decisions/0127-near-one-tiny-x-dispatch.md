# ADR-0127: the near-1 tiny-x family (exp and cosh short-circuit; cos and the reciprocals were already handled)

- **Status**: accepted
- **Date**: 2026-07-05

## Context

pf-767j continued pf-7nnw (ADR-0121, which fixed the near-0 `x ± c·x³`
family: `sin`, `tan`, `asin`). The remaining deep-tiny directed-mode
1-ulp defects were framed as a near-1 family where the base is 1 (not x)
and the correction is `x²` or linear (not `x³`): `cos(x) = 1 − x²/2`,
`cosh(x) = 1 + x²/2`, `exp(x) = 1 + x`, plus `sec = 1/cos`, with `csc`
and `cot` (`~1/x`) flagged to probe. The mechanism: for x tiny enough
that the result rounds to 1 (or its neighbour), the eval past the Ziv
guard cap collapses to exactly 1, and the nearest modes are right but the
directed modes certify 1 where the neighbour is due (the saturation-
analogue of `tanh`/`erf`, and the near-0 analogue of ADR-0121).

Reproducing RED first corrected the scope. Empirically, at `x = 2^-1200`
(past the cap) and target 53:

| kernel | pre-fix | why |
|---|---|---|
| `exp` (both signs) | WRONG | plain `ziv_round`, no depth mechanism |
| `cosh` | WRONG | plain `ziv_round`, no depth mechanism |
| `cos` | correct | `ziv_round_with_depth` + `reduction_depth_hint` |
| `sec`, `csc`, `cot` | correct | same depth-hint driver |

`cos` and the reciprocals dispatch through the trig-reduction path, whose
`ziv_round_with_depth` depth hint (ADR-0103) already grows the working
precision to the near-0 / near-grid proximity depth for any input, so
they correctly round `cos(2^-1200)` to `pred(1)` under `TowardZero`,
`csc(2^-1200)` to `succ(2^1200)` under `TowardPositive`, and so on. The
bead's "confirmed cos red" premise did not hold on the current tree; the
`ziv_round_with_depth` wiring predates it. So the fix is exactly `exp`
and `cosh`, which use plain `ziv_round`.

## Decision

Short-circuit the two plain-`ziv_round` kernels to
`round_with_infinitesimal(1, Positive, subtracts_magnitude, target,
mode)`, matching the `expm1`/`sinh`/`tanh` tiny-x pattern (ADR-0059).

- **`exp`: `e_x ≤ −(target+2)`.** The correction is the LINEAR term `x`
  (`exp(x) = 1 + x + …`), so the perturbation exponent is `e_x`, not
  `2·e_x`. The threshold puts `|x| ≤ 2^{−(target+1)}`, strictly below the
  half-ulp above 1 (`2^{−target}`) and below the half-ulp below 1
  (`2^{−target−1}`), so both signs round correctly. `subtracts_magnitude
  = x.is_sign_negative()`: `exp(x) < 1` iff `x < 0`.
- **`cosh`: `2·|e_x| ≥ target+2`.** The correction is `x²/2`
  (`cosh(x) = 1 + x²/2 + …`), so the deviation exponent is `2·e_x − 1`.
  The threshold puts the deviation upper bound `2^{2·e_x+1}` at
  `≤ 2^{−target−1}` (half the ulp BELOW 1, the tighter of the two
  asymmetric half-ulps around 1), so `cosh` rounds to 1 (nearest) or
  `succ(1)` (directed-up). `subtracts_magnitude = false` (grows).

The threshold derives the `x²`-vs-linear distinction from the series, not
from analogy to ADR-0121's `x³` threshold `e ≤ −(target+2)`, which is
wrong for a quadratic or linear correction. Placed after the special-case
and exponent-rim dispatch and disjoint from the `e_x ≥ 62` extreme-x rim.

## Consequences

- `exp` and `cosh` are correctly rounded under every mode for tiny x, at
  any depth (the short-circuit is O(1), input-independent). The 20-row
  sweep (`tests/regression_review_2026_07_05_r53.rs`) brackets each
  threshold under all five modes; the deep rows are RED pre-fix.
- **Inversion: the bead over-scoped; reproduce-RED-first caught it.** Four
  of the six named kernels (`cos`, `sec`, `csc`, `cot`) were already
  correct via the depth-hint driver. Fixing them would have added a
  redundant short-circuit and churned correct code. The audit is recorded
  as `cos`/`sec` control rows so a future reader sees they were checked,
  not missed.
- **Inversion: the short-circuit must not over-fire.** A threshold too
  loose would divert an x where the driver gives a DIFFERENT correct
  answer. The boundary rows (`exp` at `e = −53…−56`, `cosh` at
  `e = −26…−29`) sit on both sides of the threshold and pass under all
  modes both pre-fix (driver, within cap) and post-fix (short-circuit) —
  so the two agree in the overlap, and only the deep rows past the cap
  needed the fix.
- **Inversion: the half-ulp around 1 is asymmetric.** The ulp just below 1
  is `2^{−target}`; just above, `2^{−target+1}`. `cos`/`cosh` (deviation
  toward and away from 1) bind on the SMALLER (below-1) half-ulp, so the
  `x²` threshold is `target+2`, not `target+1`; a `+1` margin would have
  mis-rounded `cos` at the boundary (`cos` is handled elsewhere, but the
  same margin governs any future near-1-from-below kernel).

## References

- pf-767j (the bug), pf-7nnw / ADR-0121 (the near-0 sibling), ADR-0059
  (the tiny-x `round_with_infinitesimal` pattern), ADR-0103
  (`ziv_round_with_depth` depth hint that already covers cos and the
  reciprocals), ADR-0096/0107 (`round_with_infinitesimal` rim safety).
- `tests/regression_review_2026_07_05_r53.rs`; oracle generator
  `scratchpad/gen_r53.py`.
