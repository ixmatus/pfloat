# ADR-0090: complex multiply and divide, and the rounding-mode and Status API

- **Status**: accepted
- **Date**: 2026-06-08

## Context

Slice C3 of pfloat-complex lands `Complex::mul` and `Complex::div`. Both
must be componentwise correctly rounded (the real and imaginary parts each
correctly rounded under the target real rounding mode, MPC's model). Two
problems sit on top of the arithmetic: the divide is genuinely harder than
the multiply, and the crate must settle two API shapes that recur across
every operation (the rounding mode and the `Status`).

## Decision

### Multiply: componentwise correct rounding for free

`(a + bi)(c + di) = (ac − bd) + (ad + bc)i`. Each component is one fused
two-product (the C1 primitive, ADR-0088): `re = a.mul_sub_mul(c, b, d)`,
`im = a.mul_add_mul(d, b, c)`. Because `mul_sub_mul` / `mul_add_mul` are
correctly rounded with a single rounding and no Ziv loop, the product is
componentwise correctly rounded by composition. Catastrophic cancellation in
`ac − bd` is exact (ADR-0088 Proof 1), so `z·conj(z) = |z|²` lands with an
exactly-zero imaginary part and no spurious `INEXACT`.

### Divide: a directed-pair enclosure Ziv loop

`(a + bi)/(c + di) = [(ac + bd) + (bc − ad)i] / (c² + d²)`. The quotient is
**not** correctly rounded by dividing a separately-rounded numerator by a
separately-rounded denominator (two roundings). Forming the exact numerator
then dividing once is correct in principle but **infeasible**: when one
product dominates the other (`ac ≫ bd`), the exact sum's bit-length tracks
the exponent gap and can exceed any representable precision.

The divide therefore runs a directed-pair enclosure Ziv loop, entirely on
pfloat's existing public API (no new core primitive needed):

1. At working precision `w = p + guard`, bracket the numerator and the
   denominator with their directed fused-two-product pairs:
   `N_lo = mul_add_mul(.., TowardNegative)`, `N_hi = (.., TowardPositive)`,
   and likewise `D_lo`, `D_hi` for `c² + d²` (which is `≥ 0`).
2. Form the quotient enclosure with directed division, then round both ends
   to the output precision `p` under the target mode.
3. Accept when the two ends agree (equal value **and** equal sign); else
   grow the guard (`64, 128, 256, 512, 1024`, capped at five iterations, the
   MPFR measure-zero caveat).

The working precision stays bounded for any exponent gap, because the
directed pair *brackets* the true value rather than representing it exactly.
Because `FixedFloat<PREC>` cannot hold `w > PREC`, the kernel runs in
`BigFloat` and the generic divide bridges through `RealScalar::to_big` /
`from_big`.

Two corrections to the adversarially-reviewed design, both load-bearing,
both found by re-deriving rather than trusting the review's "sound" verdict
(the project's verify-the-verdict discipline):

- **Sign-aware quotient enclosure.** The reviewed `[N_lo/D_hi, N_hi/D_lo]`
  holds only for a non-negative numerator. The imaginary numerator
  `bc − ad` is routinely negative or zero-straddling, where the correct
  enclosure (with `D ≥ 0`) is
  `lo = N_lo / (N_lo<0 ? D_lo : D_hi)`,
  `hi = N_hi / (N_hi<0 ? D_hi : D_lo)`.
  The naive formula would silently mis-round negative imaginary parts.
- **Exact-zero cancellation.** When the numerator cancels to exactly zero
  (the imaginary part of `z/z`, `bc − ad = 0`), the directed pair brackets
  `[−0, +0]` (a cancelling difference is `−0` under `TowardNegative`, `+0`
  otherwise), and the sign-aware convergence test never agrees, looping to
  the cap and over-reporting `INEXACT`. This is the cancellation-to-zero
  class (the tanh-defect family). The kernel special-cases an exact-zero
  numerator over a positive denominator as an exact signed zero, `+0` except
  in `TowardNegative`. A zero numerator over a zero denominator stays `0/0 =
  NaN + INVALID`.

`±0` convergence uses value equality plus a sign check, because IEEE
comparison treats `+0 == −0`; the Ziv test must not converge on the
wrong-signed zero.

### Componentwise division by zero is `NaN`, for now

For a zero divisor `c + di = 0`, each numerator (`ac + bd`, `bc − ad`) is
also zero (it is built from the divisor's components), so each component is
`0/0 = NaN + INVALID`. The C99 Annex G refinement (`z/0` for nonzero `z`
yielding a complex infinity) is a later slice (C4, the branch-cut and
special-value work); C3 is pure componentwise division.

### The rounding-mode and `Status` API: single, merged, by default

- **Single `RoundingMode`.** Every operation takes one
  `mode: RoundingMode`, applied to both components. A
  `ComplexRoundingMode { re, im }` pair (rounding the parts under different
  modes) is the spec-complete form but has no demonstrated consumer; it is
  deferred behind an additive variant, to be added when the complex ball or
  a 1788 face proves it load-bearing.
- **Merged `Status`.** Operations return one `Status`, the OR-monoid
  combination of the two component statuses (`INEXACT` if either part
  rounded, `INVALID` / `DIV_BY_ZERO` if either raised them). A per-component
  `ComplexStatus { re, im }` is deferred behind an accessor, for the same
  reason. The merged form is the honest default: a caller checking
  `INEXACT` wants "was anything rounded," and `|=` composition is exactly
  the IEEE sticky-flag discipline pfloat already uses.

## Consequences

- `Complex::mul` is componentwise correctly rounded by composition over the
  C1 primitive, with no Ziv loop; `Complex::div` is componentwise correctly
  rounded by the directed-pair Ziv loop, bounded for any exponent gap.
- The single-mode, merged-`Status` shape is the crate-wide default; the pair
  forms are additive, deferred, and recorded here so adopting one later is a
  conscious 1.x addition, not a reshape.
- The two design corrections (sign-aware enclosure, exact-zero cancellation)
  are the difference between a correct divide and a silently wrong one; both
  are pinned by unit tests (the all-five-modes real-axis cross-check against
  scalar division, `z/z` exact, `1/i = −i`).

## Related

- ADR-0088: the fused two-product primitive multiply composes from and
  divide brackets with.
- ADR-0089: the sealed `RealScalar` trait the bridge (`to_big`/`from_big`)
  extends.
- ADR-0091 (forthcoming): magnitude, phase, and the elementary functions
  with C99 Annex G branch cuts, where the complex-infinity division
  refinement also lands.
- Plan: `~/.claude/plans/plan-tower-expansion-scope-goofy-raven.md` (slice
  C3) and the scoping doc's decision 10 (the rounding-mode and Status calls).
