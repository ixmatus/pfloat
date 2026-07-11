# 0129. `exp_round` returns `Result`, matching the `*_round` family

## Status

Accepted (2026-07-11).

## Context

Every mode-aware `*_round` kernel on `BigFloat` validates its
`target_precision` argument and returns `Result<(Self, Status),
BuildError>`, yielding `Err(BuildError::PrecisionZero)` when
`target_precision == 0`: `ln_round`, `log2_round`, `log10_round`,
`log1p_round`, `sin_round`, `cos_round`, `tan_round`, the reciprocal trig,
`asin_round`, `pow_round`, `sinh_round`, `cosh_round`, `tanh_round`,
`asinh_round`, and the rest.

`exp_round` was the lone exception. It returned `(Self, Status)` and guarded
precision only with `debug_assert!(target_precision >= 1)`. In a release
build the assert is compiled out, so `exp_round(0, mode)` fell through into
`exp_kernel`, which panics on the first `.expect("precision >= 1")` it
reaches. The features/no_std and review-tail work surfaced this as pf-291u
(review finding `elementary/type/4`): a public entry point that panics in
release on a well-typed input, inconsistent with every sibling.

The inconsistency is a real footgun: a caller that handles `Err` uniformly
across the elementary surface still gets a panic from `exp_round`, and the
`debug_assert` hides it from debug-built test suites.

## Decision

`exp_round` now returns `Result<(Self, Status), BuildError>` and returns
`Err(BuildError::PrecisionZero)` when `target_precision == 0`, exactly like
its siblings. The `exp` convenience method (`self.precision` is always `>= 1`
by the `BigFloat` invariant) unwraps with `.expect("self.precision >= 1 by
invariant")`, the same shape the other convenience wrappers use.

Every internal caller passes a precision that is nonzero by construction
(a validated `target_precision`, a working precision `>= target`, or a
literal), so each adapts with `.expect(...)` naming the invariant rather
than propagating the error:

- `exp10`, `expm1`, `cosh`, `sinh` (the elementary composers);
- `pfloat-complex` `cexp`/`clog`/`complex`;
- `pfloat-libm`: `exp` moves from the macro's `direct_sat` arm to
  `result_sat` (the arm that `.expect`s the kernel `Result`). `exp` was the
  last user of the tuple-returning `direct`/`direct_sat` arms, so both are
  removed and the macro doc updated; every kernel the macro emits now takes
  the `Result` path.

Alternatives considered and rejected:

- **Non-breaking sentinel** (`(qNaN, INVALID)` for precision 0): keeps the
  signature but is semantically wrong (precision 0 is not an IEEE invalid
  operation) and leaves `exp_round` the odd one out against the `Err`
  siblings.
- **Document the panic as contract**: cheapest, but enshrines the one
  panicking outlier and the cross-surface inconsistency the finding is
  about.

## Consequences

- The elementary `*_round` surface is now uniform: one precision-validation
  contract, no release panic on a well-typed input. Illegal input is
  surfaced in the type, per the total-function preference.
- This is a **breaking change** to `BigFloat::exp_round`'s public signature.
  pfloat is not yet published to crates.io, so no released contract breaks
  today; the next publish that carries this change requires a major version
  bump (or bundling with the other pre-publish breaking changes, e.g. the
  pf-9fq3 `Status` deserialize hardening). Recorded here so the bump is not
  forgotten at release time.
- The `pfloat-libm` `unary!` macro loses its `direct` and `direct_sat` arms.
  A future tuple-returning kernel would re-add a variant; none exists.

## References

- pf-291u (review finding `elementary/type/4`), epic pf-8iji.
- Sibling contract: `src/math/ln.rs` `ln_round`, `src/big.rs`
  `BuildError::PrecisionZero`.
- `tests/regression_nt21_specials.rs` (adjacent precision-0 contract for
  digamma/lgamma); the `exp_round` precision-0 guard is pinned in
  `tests/regression_review_2026_07_11_exp_round_result.rs`.
