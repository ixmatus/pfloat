# ADR-0072: Exact power-of-two scaling (`scale_by_pow2`)

- **Status**: accepted
- **Date**: 2026-06-06

## Context

Phase 4 builds rigorous enclosure arithmetic (`pfloat-ball`) on top of
the correctly-rounded scalar base. A ball op computes a midpoint at
working precision and then bounds a radius; almost every radius
manipulation, and the exact `lower`/`upper` inf-sup conversion, needs to
multiply a value by a power of two exactly and in `O(1)` time. The same
operation is IEEE 754-2019 §5.3.3 `scaleB` restricted to an exact `2^k`
scale, and it is useful well beyond the ball (decompositions of the form
`m × 2^e`, fast normalization, `ldexp`-style scaling).

Multiplying a binary float by `2^k` only shifts the unbiased exponent;
the mantissa, sign, and precision are untouched. pfloat already performs
this shift internally (the arithmetic kernels add to the `Normal`
exponent), but no public surface exposes it. The deliberate absence of a
raw-parts constructor (`big.rs`: "pfloat does not expose a converse
constructor from raw parts") means a caller outside the crate cannot
synthesize an exact scaling from `parts()`; it has to go through a full
multiplication by a constructed `2^k`, which allocates a mantissa and
routes through the rounding pipeline for an operation that is always
exact. The v1.0 API (frozen under ADR-0054) has no exact scaling
primitive.

A design wrinkle specific to pfloat: there is no `emax` or `emin`. The
exponent is an `i64` that **saturates** rather than producing `±∞` or a
subnormal. `ops::mul` already fixed this contract (the
`mul_extreme_exponent_saturates_not_panics` regression, fuzz-found via
Airy `bi_prime`): an exponent computed past the `i64` range clamps to
`i64::MAX`/`i64::MIN` as a finite value and raises `OVERFLOW`/`UNDERFLOW`,
never `±∞`. A scaling primitive must inherit exactly that contract so the
two operations agree.

## Decision

Add `scale_by_pow2(&self, k: i64) -> (Self, Status)` to `BigFloat`
(implemented in `src/ops/scale.rs`) and the same signature to
`FixedFloat<PREC>` (delegating through `to_big`). A `_with_flags` sibling
mirrors the house pattern on `BigFloat`. Purely additive, so the v1.0
freeze is preserved.

Semantics:

- **Finite non-zero**: the exponent shifts by `k`, computed in `i128`
  (the sum of two `i64`s cannot overflow `i128`) and clamped to the
  `i64` range. A shift past `i64::MAX` returns `i64::MAX` with
  `Status::OVERFLOW`; a shift past `i64::MIN` returns `i64::MIN` with
  `Status::UNDERFLOW`. Every in-range shift is exact, `Status::OK`, with
  the mantissa cloned byte-for-byte and the precision unchanged. This is
  the same saturating, no-`emax` contract as `ops::mul`.
- **`±0` / `±∞`**: returned unchanged (`±0 × 2^k = ±0`,
  `±∞ × 2^k = ±∞`), sign preserved, `Status::OK`.
- **NaN**: a quiet NaN propagates (sign and payload preserved); a
  signaling NaN raises `Status::INVALID` and is quieted, because
  `scaleB` is a general-computational operation and signals on sNaN like
  the other kernels (`sqrt`, `mul`).

`scale_by_pow2` takes no `target_precision` and returns no `BuildError`:
the operation never rounds and never changes precision, so there is no
fallible precision argument.

## Consequences

- The `pfloat-ball` radius operations and the exact `lower`/`upper`
  inf-sup conversion get an `O(1)`, allocation-light, exact primitive
  instead of a full multiplication. The clone of `self` is the only
  unavoidable cost; the exponent shift itself is constant-time.
- The validated-construction invariant is preserved: no raw-parts
  constructor is introduced, so the top-bit-set normalization, the
  `limbs_for(precision)` storage shape, and the precision-bound payload
  length stay type-checked invariants. The internal exponent mutation is
  surfaced safely, not by exposing the field.
- Saturation makes the operation total: no input panics, no input
  produces `±∞` from a finite operand (matching the no-`emax` model),
  and the `OVERFLOW`/`UNDERFLOW` flags let a caller detect the
  measure-zero saturation cases. For the ball this matters: a radius
  scaling that saturates upward (`OVERFLOW`) stays a sound upper bound,
  and one that saturates downward (`UNDERFLOW`) is flagged so the
  caller never silently treats a clamped radius as exact.
- The agreement with `mul` is testable and tested
  (`scaling_agrees_with_multiplication`), so the two paths cannot drift.

## Related

- Plan: `~/.claude/plans/melodic-knitting-ripple.md` (slice 1); design
  reference `~/.claude/plans/pfloat-phase4-scoping.md`.
- Beads: `pf-icgj.1` (under epic `pf-icgj`).
- Other ADRs: the saturating no-`emax` exponent contract originates in
  `ops::mul`; the deliberate absence of a raw-parts constructor is the
  `parts()` accessor decision (ADR-0016). `next_up`/`next_down`/`ulp`
  (ADR-0073) is the sibling slice-2 primitive.
