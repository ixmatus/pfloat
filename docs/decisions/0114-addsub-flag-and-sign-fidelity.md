# ADR-0114: add/sub result fidelity — directed-mode sign mirroring and exponent-saturation flag parity

- **Status**: accepted
- **Date**: 2026-07-03

## Context

Two independent fidelity gaps in `src/ops/addsub.rs`, both from the
2026-06 review (epic pf-8iji), both mechanical:

1. **pf-egxm — `0 − x` mis-rounds under directed modes.** `zero_plus_finite`
   handled `0 ± b` and `a ± 0` by rounding the operand at its *stored*
   sign and flipping the sign afterward. Directed rounding is relative to
   the true signed value, so when the effective (result) sign differs
   from the operand's stored sign — the `0 − b` case — the round went the
   wrong way: `sub_round(+0, 2^60+1 @64, 53)` under TowardPositive rounded
   the positive `2^60+1` toward +∞ (to `2^60+256`) then negated, yielding
   `−(2^60+256)` where `−2^60` is due (TowardPositive on a negative value
   rounds toward less-negative). TowardNegative mirrored the error.

2. **pf-kh3z — add/sub exponent saturation dropped its flag.** When the
   result exponent exceeds the i64 range it is clamped (pfloat has no
   `emax`; `i64::MAX`/`i64::MIN` is a saturated finite), but addsub
   clamped silently with Status OK, while `mul`/`div`/`fma`/`remainder`/
   `scale` all raise OVERFLOW/UNDERFLOW on the same clamp. `a + a` with
   `a` at exponent i64::MAX returned Status OK.

## Decision

1. `zero_plus_finite`: apply the effective sign to a clone of the operand
   *before* `round_to_precision`, so the directed mode sees the correct
   direction. The add case and the `a ± 0` cases are unaffected (the
   effective sign equals the stored sign there, so the clone-and-apply is
   an identity).

2. At the exponent clamp, mirror `mul` exactly: compute an
   `exp_saturation` status (OVERFLOW above `i64::MAX`, UNDERFLOW below
   `i64::MIN`, OK in range) and OR it into the final status after the
   rounding pipeline. Below the rims `i64::try_from` succeeds and the
   status is unchanged, so common-path results are untouched.

## Consequences

- `0 − x` is correctly rounded under every mode; the exponent-saturation
  flag now matches every sibling kernel, so downstream sticky-flag
  consumers see a consistent OVERFLOW/UNDERFLOW contract across all ops.
- No value change on any in-range result: the saturation flag is OK below
  the rims, and the sign-apply is an identity except for `0 − x`.
- Verified against MPFR: `differential_add` (5/5) and `differential_sub`
  (4/4) unchanged.

### Inversion (failure paragraphs considered)

- *"Applying the sign before rounding breaks the `a ± 0` / add cases."*
  Refuted: in those cases the effective sign equals the operand's stored
  sign, so the applied sign is identical and the round is unchanged
  (confirmed by the differential lanes).
- *"The saturated result is inf, so OVERFLOW is redundant."* Whether the
  clamp renders as a huge finite or an inf via the top-rim path, the
  transition past `i64::MAX` is an inexact overflow that the sticky flag
  must record; the sibling kernels set it, and a silent OK breaks flag
  parity for consumers that branch on OVERFLOW.

## References

- pf-egxm, pf-kh3z (epic pf-8iji); reproducers C1, D1 in the review
  harness `results.txt`.
- `src/ops/mul.rs` (the `exp_saturation` pattern mirrored here),
  `src/ops/addsub.rs` (`zero_plus_finite`, `add_finite_finite`).
- `tests/regression_review_2026_06_10.rs`:
  `zero_minus_finite_mirrors_directed_mode`,
  `addsub_exponent_saturation_flags_overflow`.
