# pfloat-libm kernel list

This document is the per function record of what pfloat-libm v0.1 ships, where
each kernel is implemented, and how strongly it is verified. It is the libm
phase's first deliverable, and its load bearing purpose is to cite pfloat's
ADR-0032 against the six reciprocal and root kernels so the trivial alias
reflex never reintroduces a composition the status table would have to label as
merely faithful.

## How pfloat-libm computes

pfloat-libm is a thin shell over pfloat. Each function widens its hardware float
argument to a `BigFloat` exactly, evaluates pfloat's correctly rounded kernel at
a working precision, and rounds back to the hardware width under an outer Ziv
loop: it forms an enclosure of the true value from the kernel result and its
error bound, and commits a hardware float only when both ends of the enclosure
round there. The second rounding is therefore never blind, which is what closes
the `BigFloat` to float double rounding gap that a naive widen, compute, round
shell leaves open.

The arbitrary precision mathematics lives in pfloat, not here. The four net new
unary kernels and the two net new multi argument kernels below are implemented
as direct primary kernels in pfloat 1.1 (the conversion API and the six kernels
are the pfloat side of Phase 2); pfloat-libm wraps them exactly as it wraps the
elementary surface pfloat already shipped at 1.0.

## Verification tiers

- **Exhaustive f32.** Every one of the 2^32 `binary32` inputs is evaluated and
  compared against an independent oracle (MPFR via `rug`); the gate is zero
  mismatches under NearestEven, and the directed modes wherever the enclosure
  determines them. This is the claim no competing pure Rust libm can make.
- **Differential f64 plus worst case.** The 2^64 `binary64` space cannot be
  enumerated, so f64 rests on a large structured differential sample against
  MPFR plus the Lefevre and Muller hard to round vectors as adversarial seeds.
- **Differential multi argument.** Binary functions cannot be exhausted over
  either width and stay on differential plus worst case vectors for both.

## Elementary surface (inherited from pfloat 1.0, correctly rounded)

These already carry a correctly rounded status row in pfloat's own f32 oracle
sweep. pfloat-libm exposes them through the shell at both widths.

| Function | f32 tier | f64 tier |
| --- | --- | --- |
| `exp`, `exp2`, `exp10`, `expm1` | exhaustive | differential |
| `ln`, `log2`, `log10`, `log1p` | exhaustive | differential |
| `sqrt` | exhaustive | differential |
| `sin`, `cos`, `tan` | exhaustive | differential |
| `asin`, `acos`, `atan` | exhaustive | differential |
| `sinh`, `cosh`, `tanh` | exhaustive | differential |
| `asinh`, `acosh`, `atanh` | exhaustive | differential |

## Net new kernels (pfloat 1.1, direct per ADR-0032)

ADR-0032 (pfloat) is the gate against adding these as a few line composition
over the existing surface. The decision is recorded once there and cited here
against each function: `direct kernel required, not aliased`. The composition is
wrong for hard to round inputs, not merely slower.

### `cot`, `sec`, `csc` (trig reciprocals)

- **Why not an alias.** A correctly rounded `tan(x)` followed by a correctly
  rounded reciprocal is not a correctly rounded `cot(x)`: the two roundings
  compose into as much as one ULP of error in the reciprocal direction for hard
  to round inputs.
- **Direct kernel.** Compute `sin` and `cos` at an inflated working precision,
  reciprocate, and apply a single Ziv round at the target. The forbidden step is
  composing two target precision correctly rounded operations, not the inflated
  then rounded form. Near pole inputs (where `sin` or `cos` approaches zero) use
  the cancellation aware Ziv path.
- **Source.** DLMF section 4.14 for the reciprocal definitions and the range
  reduction shape; CRlibm and Sun fdlibm as behaviour oracles.
- **Tier.** Unary: exhaustive f32, differential f64.

### `cbrt` (cube root)

- **Why not an alias.** `cbrt(x) = pow(x, 1/3)` cannot be correctly rounded over
  the reals because `1/3` has no exact float representation, so the exponent fed
  to `pow` is already wrong.
- **Direct kernel.** Sign handling on the real cube root of negatives, then
  `|x|^(1/3)` by Newton or Halley iteration at working precision, with a Ziv
  round at the target.
- **Source.** Standard real cube root iteration; fdlibm as behaviour oracle.
- **Tier.** Unary: exhaustive f32, differential f64.

### `hypot` (Euclidean magnitude)

- **Why not an alias.** The naive `sqrt(x*x + y*y)` overflows or underflows on
  inputs where the true magnitude is finite and representable, and the squaring
  discards half the input precision before the root can use it.
- **Direct kernel.** Scale by `max(|x|, |y|)`, evaluate `max * sqrt(1 +
  (min/max)^2)` at working precision, Ziv round at the target. IEEE 754-2019
  section 9.2.1 fixes the special cases, including `hypot(infinity, NaN) =
  infinity`.
- **Source.** IEEE 754-2019 section 9.2.1.
- **Tier.** Binary: differential plus worst case at both widths.

### `rootn` (integer nth root)

- **Why not an alias.** Same inexact exponent trap as `cbrt`, with the added
  domain question of even `n` over negative `x`.
- **Direct kernel.** Integer root Newton iteration with the IEEE sign rules: odd
  `n` of a negative argument is the real root, even `n` of a negative argument is
  NaN with INVALID raised. IEEE 754-2019 section 9.2.
- **Source.** IEEE 754-2019 section 9.2.
- **Tier.** Binary (argument plus integer order): differential plus worst case.

## Out of scope for v0.1

- Multi argument elementary functions beyond `hypot` and `rootn` (`pow`,
  `atan2`): already in pfloat, exposed in a later pfloat-libm minor; differential
  rigor only.
- Special functions (gamma, erf and erfc, the integrals, Airy, Bessel, zeta):
  already in pfloat and cheap to expose through the same shell, deferred to
  pfloat-libm v0.2 because their f32 verification needs the Arb oracle stack back.
  v0.1 is deliberately MPFR only.

## References

- pfloat ADR-0032 (libm reciprocal and root kernels stay direct, not aliased).
- pfloat ROADMAP, Phase 2 (libm spinoff).
- DLMF, IEEE 754-2019 sections 9.2 and 9.2.1, as cited per kernel above.
