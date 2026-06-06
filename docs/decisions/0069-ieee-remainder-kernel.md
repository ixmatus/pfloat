# ADR-0069: IEEE 754-2019 remainder kernel

- **Status**: accepted
- **Date**: 2026-06-06

## Context

pfloat's arithmetic surface (add, sub, mul, div, sqrt, fma) shipped
without a remainder operation. IEEE 754-2019 §5.3.1 defines
`remainder(x, y) = x - n·y`, where `n` is `x / y` rounded to the
nearest integer with ties to even; the result is exact and satisfies
`|remainder| <= |y|/2`. Its absence was a genuine gap in the spec
surface.

The gap became load bearing for the num-traits adoption work
(pf-a4jh): `num_traits::Num` requires `core::ops::Rem`, and a float
`%` needs a remainder kernel. So the kernel is both a spec-completeness
item and the prerequisite that unblocks the `Num` and `Signed` trait
impls.

The naive implementation `x - trunc(x / y)·y` is a denial-of-service
vector at arbitrary precision: the quotient's integer part can be as
large as the `i64` exponent range, so forming it is infeasible.

## Decision

Implement IEEE 754-2019 §5.3.1 `remainder` as `BigFloat::remainder`
and `FixedFloat::remainder`, returning `(value, Status)`. The result is
exact (it never rounds), at a precision of `max(px, py)`. Wire
`core::ops::Rem` (`%`) and `RemAssign` to it for both types.

pfloat's `%` is therefore the IEEE remainder, not C `fmod`. The two
differ in how the quotient rounds (nearest-even versus toward zero).
The choice is documented at the operator impl rather than aliasing `%`
to a non-IEEE operation; an IEEE library's `%` being the IEEE remainder
is the least surprising reading for this crate, and `num_traits::Num`
does not mandate `fmod` semantics.

The algorithm reduces the mantissa integers modulo each other rather
than forming the quotient. With `|x| = Mx·2^ax` and `|y| = My·2^ay`,
the truncated remainder is `(X mod Y)·2^s` for `s = min(ax, ay)`:

- When `x` dominates (`ax >= ay`), the scaling factor `2^(ax-ay)` is
  reduced by modular exponentiation, `2^(ax-ay) mod My`, in `O(log)`
  multiplications. It is never materialized.
- When `y` dominates, the early exit `2|x| < |y| → x` (the nearest
  integer quotient is zero) bounds the opposite shift to at most
  `px - py + 1` bits, so that divisor is materialized safely.

The round-to-nearest-even adjustment compares `2R` with the scaled
`|y|`; on an exact tie it consults the truncated quotient's parity,
recovered with a second modular exponentiation modulo `2·My`, paid only
on ties. The signed-zero result of an exact multiple carries the
dividend's sign, per IEEE.

## Consequences

The arithmetic surface now covers IEEE §5.3.1, and the num-traits work
(pf-a4jh) can implement `Num` and `Signed`. The kernel is verified
bit-for-bit against MPFR's `mpfr_remainder` in
`tests/differential_remainder.rs` across random integer pairs,
power-of-two exponent gaps (driving the modular-exponentiation path and
the early exit), round-to-nearest-even ties, and fractional operands;
unit tests pin the special cases (`remainder(x, ±0)` and
`remainder(±∞, y)` invalid, `remainder(x, ±∞) = x`,
`remainder(±0, y) = ±0`, NaN propagation, sign rules).

The DoS-safety is structural: the reduction is `O(log)` in the exponent
gap, so a `10^9` gap costs the same handful of multiplications as a
small one.

`%` resolving to the IEEE remainder rather than `fmod` is a deliberate
divergence from `f64`'s operator; it is documented at the impl. A
`fmod`-style truncated remainder can be added later as a named method
if a caller needs the `f64`-consistent operation.

The kernel adds two limb helpers (`modpow2`, `sub_limbs`) local to the
module. It carries no Kani harnesses yet; the differential and unit
coverage stand in, and totality proofs are a possible follow-up.

## Related

- Plan: `~/.claude/plans/melodic-knitting-ripple.md` (Phase 3, slice C3c)
- Issues: `pf-2138` (discovered from `pf-a4jh`)
- Other ADRs: builds on ADR-0001 and ADR-0002 (mantissa layout and
  bit-level precision), ADR-0006 (`i64` exponent); unblocks the
  num-traits decision (ADR-0070)
