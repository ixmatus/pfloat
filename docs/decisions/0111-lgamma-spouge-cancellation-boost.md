# ADR-0111: route lgamma's Spouge regime through cancellation_boosted (the pf-rlrb certified-wrong-value sibling of ADR-0110)

- **Status**: accepted
- **Date**: 2026-07-03

## Context

Fixing the digamma interior-zero defect (pf-0r1l, ADR-0110) surfaced a
pre-existing sibling in `lgamma` (pf-rlrb, P1, epic pf-8iji). The two
share a kernel and a mechanism.

`gamma_stirling::spouge_lgamma_scaled` returns `(value, operand_scale)`:
its `ln S` step (S the Spouge sum) is a **near-total internal
alternating cancellation** whose depth grows with the argument
(~0.1·`working` at `z ≈ 2.5`, ~0.4·`working` at `z ≈ 1e6`). The value
it returns is therefore accurate to only `working − (that
cancellation)` bits, and the returned scale reports the depth so a
caller can charge it. Its docstring states outright: *"The caller MUST
re-drive through `cancellation_boosted`."*

`lgamma`'s positive branch honoured that contract **only inside the
positive root windows** `[3/4, 5/4] ∪ [7/4, 9/4]`
(`lgamma_at_w`, ADR-0097). For `x` outside those windows the branch
fell through to `lgamma_positive_at_w(x, z_min, working_prec).0`, which
for `working_prec > STIRLING_REACH_THRESHOLD = 600` dispatches to
`spouge_lgamma_scaled` and **discarded `.1`**. The Ziv half-width model
`|y|·2^-(working − guard)` was then violated by the uncharged
cancellation and a **wrong value certified**.

Run-verified (mpmath 1.4.1 @6000 bits, exact input `5/2`): `lgamma(5/2)`
at target 1024 returned a value diverging from truth at bit ~857 where
1024 correctly-rounded bits were promised. The input `5/2` is exactly
representable, so it encodes no proximity depth — the cancellation is
**purely internal to the Spouge sum**, the mechanism the digamma
sibling already handles.

The defect was pre-existing and untested: `differential_gamma` and its
lane siblings cap at p256 (working ≤ 320 < 600), so no lane exercised
the Spouge path. `gamma` (`exp(lgamma)`) and `beta`
(`lgamma(a)+lgamma(b)−lgamma(a+b)`) call `lgamma_round` at their own Ziv
working precision, so both certified wrong values transitively at high
targets (verified red-before/green-after: `gamma(5/2)@1024`,
`beta(5/2,3/2)@1024` both failed with the fix stashed).

The digamma probe pf-2thy ("digamma lacks lgamma's >600-bit Spouge
dispatch; 920-bit ceiling + 2^28-iteration shift loops") was already
resolved by ADR-0110 — digamma's Spouge dispatch and
`cancellation_boosted` routing exist — and is closed by the
high-precision guard added here.

## Decision

Extend `lgamma_at_w`'s positive-branch routing condition to mirror the
digamma sibling exactly:

```rust
if in_positive_root_window(x) || working_prec > STIRLING_REACH_THRESHOLD {
    return super::ziv::cancellation_boosted(working_prec, |w| {
        lgamma_positive_at_w(x, z_min_for_target(w), w)
    });
}
lgamma_positive_at_w(x, z_min, working_prec).0
```

`cancellation_boosted` (geometric growth, ADR-0110) re-runs the closure
at increasing `w`, charging the reported scale until the value carries
`working` accurate bits. For an exactly-representable input like `5/2`
the probe sees the true input on the first iteration and converges in
~2 Spouge evaluations at ≤ 2·`working` — the same cost profile the
digamma sibling accepted, and input-proportional (the DoS-budget
posture). Below the threshold the shift-Stirling path has no sum
cancellation and runs directly, so the `differential_gamma` lane
(p ≤ 256) keeps its fast path unchanged.

No new kernel, no new dispatch, no API change: `spouge_lgamma_scaled`
already returned the scale, and `lgamma_positive_at_w` already
propagated it. The fix is a one-condition change that makes lgamma
honour the contract its scaled kernel documents.

Verification (`tests/regression_review_2026_06_10.rs`, mpmath 1.4.1
oracle): `lgamma(5/2)@1024` in NE/TZ/TP, a `gamma`+`beta` transitive
guard at 1024, a `digamma(5/2)@1024` guard closing pf-2thy, and a
Stirling-path control at p256 confirming the fast path is untouched.
The MPFR differential lanes (lgamma/gamma/beta/digamma/zeta, p ≤ 256)
confirm no fast-path regression.

## Consequences

- Positive: the last certified-wrong-VALUE defect in the gamma family
  is closed; `lgamma`, `gamma`, and `beta` are now correctly rounded at
  every target, not only where the differential lanes reached.
- Cost: high-target lgamma off the root windows now runs ~2 Spouge
  evaluations instead of one, at up to 2·working. Bounded and
  input-proportional; the differential lanes (p ≤ 256) are unaffected
  because they never crossed the threshold.
- The `spouge_lgamma_scaled` docstring's "caller MUST re-drive"
  contract is now satisfied at every call site (the crate-wide sweep
  found exactly one caller; the two remaining `_positive_at_w(..).0`
  scale-discards are both below-threshold Stirling fall-throughs with
  no sum cancellation, correct to discard).

### Inversion (failure paragraphs considered)

- *"The root-window path already boosted, so the bug is the window
  bounds, not the routing."* Refuted by reproduction: `5/2 = 2.5` is
  outside `[7/4, 9/4] = [1.75, 2.25]`, and widening the window would not
  cover arbitrary large `z` where the Spouge cancellation is largest.
  The Spouge cancellation is argument-driven and independent of the
  root windows; the routing, not the bounds, is the fix.
- *"Internal absorption (compute S at `working + working/8 + 96`) is
  simpler than re-driving."* Rejected for consistency and honesty: the
  sibling (digamma, ADR-0110) re-drives, the scale is already computed
  and propagated, and a fixed `working/8` margin under-covers the
  ~0.4·`working` cancellation at large `z` (the pf-gg96/lying-inner-
  kernel failure shape). The scale-driven iteration adapts to the
  realised depth; a fixed margin is a guess.
- *"gamma/beta have their own Spouge dispatch and need separate
  fixes."* Refuted by the crate-wide sweep: `gamma` composes
  `exp(lgamma)` and `beta` composes lgamma sums; neither has its own
  Spouge path. Both are fixed transitively, confirmed by the
  red-before/green-after guards.

## References

- pf-rlrb (this defect), pf-2thy (digamma probe closed here), epic
  pf-8iji.
- ADR-0110 (digamma sibling, the routing template), ADR-0097
  (`cancellation_boosted` and the operand-scale charge), ADR-0098
  (input-structure-aware conditioning).
- `src/math/lgamma.rs` (`lgamma_at_w`), `src/math/gamma_stirling.rs`
  (`spouge_lgamma_scaled`), `src/math/ziv.rs` (`cancellation_boosted`).
- `tests/regression_review_2026_06_10.rs`:
  `lgamma_high_precision_spouge_path_is_boosted`,
  `gamma_beta_high_precision_inherit_the_lgamma_boost`,
  `digamma_high_precision_spouge_path_is_boosted`.
