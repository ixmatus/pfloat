# ADR-0056: The six direct libm kernels (cot, sec, csc, cbrt, hypot, rootn)

Status: Accepted (2026-06-01)

## Context

ADR-0032 ruled that `cot`, `sec`, `csc`, `cbrt`, `hypot`, and `rootn`
ship as direct primary kernels, never as aliases over an
already-rounded sibling, and deferred them from the v1.0 surface to the
libm phase. Slice pf-lm1b implements all six in pfloat 1.1 so the
`pfloat-libm` spinoff can wrap them. ADR-0032 fixed the policy (direct,
not composed) but left the per-function algorithm, special-case table,
feature gate, and verification backend open. This ADR records those
choices, the non-obvious ones in particular, and serves as the
per-function kernel-list document ADR-0032 asked the libm phase to
produce.

## Decision

### 1. Reciprocal trig: inflated sin/cos, reciprocate, one Ziv round

`cot`, `sec`, and `csc` follow the pattern `tan` already proves
(`src/math/tan.rs`, ADR-0038): inside one `ziv_round` eval closure,
reduce the argument once via the shared Payne-Hanek `reduce`, evaluate
`sin_taylor` and `cos_taylor` at the Ziv working precision, form the
reciprocal or ratio at that precision, and let the driver round once to
the target. This is a direct kernel, not the forbidden composition: the
reciprocal is taken at the inflated working precision, so only one
rounding reaches the target and the double-rounding hazard ADR-0032
names does not arise. The per-quadrant identities give `cot = cos/sin`,
`sec = 1/cos`, `csc = 1/sin` with the quadrant signs. `reduce` already
scales its internal precision with the argument's binary exponent, so
large arguments inherit correct phase precision (the RC1 fix). No
`cancellation_boosted` path is needed: `cot` at its zeros reduces to a
small accurately computed quotient, and near the poles the reduced
reciprocal is a large but finite value, exactly as `tan` behaves near
its poles. The pole special cases differ from `tan`: `cot(±0)` and
`csc(±0)` are signed infinities raising `DIV_BY_ZERO` (both functions
are odd), while `sec(±0)` is exactly `1` (`sec` is even); all three at
`±∞` are `qNaN` with `INVALID`.

### 2. cbrt: exact integer cube root, no Ziv

`cbrt` is an exact-integer root in the mould of `sqrt`
(`src/ops/sqrt.rs`), not a Ziv-driven transcendental and not
`pow(x, 1/3)`. It extracts the mantissa as an integer, shifts so the
integer cube root carries `target + guard` bits with the result scale a
multiple of three (the cube-root generalization of sqrt's even-parity
shift), takes the floor cube root via the new general
`iroot_limbs(n, k)` primitive with `k = 3`, and rounds with
`pre_sticky = (remainder != 0)`. A perfect cube is exact; every other
cube root is irrational and cannot land on a half-way tie, so the
`(intermediate_precision, pre_sticky)` interface rounds all five modes
correctly with no Ziv loop and no calibration entry. `cbrt` is the real
cube root, so a negative operand is in domain (`cbrt(-8) = -2`) and the
sign is carried through; this is the central divergence from `sqrt`.

`iroot_limbs` is integer Newton from a bit-length over-estimate followed
by a `±1` correction loop that makes the floor exact by construction;
its cost is polynomial in the operand bit length and `log k`, never
linear in `k`.

### 3. rootn: bounded-complexity Newton, no exp/ln

`rootn(x, n)` reduces to the positive `|n|`-th root, computed by Newton's
iteration at the Ziv working precision, with the sign (odd roots
preserve it) and, for negative `n`, the reciprocal applied inside the
closure before one outer round. The Newton step raises `y^(m-1)` by
square-and-multiply, so the per-step cost is `O(log|n|)`, never `O(|n|)`:
a caller supplied `i32` order must not be able to drive a linear loop.
This is the same denial-of-service discipline as the parse exponent cap
(ADR-0031) and the `pow` integer-exponent cap. `rootn` does not use
`exp`/`ln`: that composition is forbidden by ADR-0032 and sits behind a
feature `rootn` does not require. The full IEEE 754-2019 §9.2.1
special-case table is implemented directly, including the order-zero
domain error, the signed-zero parity rules, the negative-order poles
that raise `DIV_BY_ZERO`, and the even-root-of-negative domain errors.

### 4. hypot: Ziv sqrt of the sum of squares, no Moler scaling

`hypot(x, y)` evaluates `sqrt(x² + y²)` at the Ziv working precision in
one closure and rounds once. The naive target-precision composition
loses half the input precision and is not correctly rounded; the direct
kernel computes at inflated precision. No scaling by `max(|x|, |y|)` is
used: that guards a fixed-width exponent field against `x²` overflowing,
but `BigFloat`'s exponent is an `i64` with saturating arithmetic, so the
field never overflows for a finite operand, and the only overflow that
can occur is the final result saturating, which the rounding pipeline
handles. The §9.2.1 special cases check infinity before NaN
(`hypot(±∞, NaN) = +∞`); a signaling NaN still raises `INVALID`, since
the infinity override fixes the value, not the §7.2 signal.

### 5. Feature gating: cbrt under `big`, the Ziv kernels under `exp-log`

`cbrt` is exact-integer with no Ziv path, so it lives under `big` and
ships in the default build alongside `sqrt`. `hypot` and `rootn` use the
`ziv_round` driver, which is gated `exp-log`; they therefore live under
`exp-log`, not bare `big`. `cot`, `sec`, `csc` live under `trig` with
their forward siblings. No new feature flags are introduced.

### 6. Verification backends

`cbrt` is MPFR-primary: the oracle sweep cross-checks it against rug's
`cbrt_ref`, and a differential lane sweeps signed integers (the negative
real-root branch). `cot`, `sec`, `csc` are Arb-primary: MPFR has no such
primitive, so they verify against Arb's native `arb_cot`/`arb_sec`/
`arb_csc` through the python-flint worker, an oracle independent of
pfloat's reduce-and-reciprocate path. The Arb worker gains a `-0` input
guard because Arb collapses signed zero and the generic path would lose
the cot/csc pole sign. `hypot` and `rootn` are two-argument, so they get
MPFR differential lanes rather than the unary sweep, matching the
atan2/pow precedent that multi-argument operations are not in the unary
status surface. rug 1.30.0 (MPFR 4.2.2) exposes the IEEE 754-2019
`rootn` directly via `root_ref` (`mpfr_rootn_ui`) and `root_i_ref`
(`mpfr_rootn_si`, signed order), so the rootn lane needs no reciprocal
synthesis and no dependency bump.

Every per-kernel error guard is `DEFAULT_ERROR_GUARD` (24 bits): the
reciprocal trio matches `tan`'s op count, `hypot` is four non-cancelling
operations, and `rootn`'s Newton iteration is self-correcting so the
error is the final step's handful of operations.

### 7. pf-gwum: cfg-gate specials-only calibration data

The slice folds in pf-gwum: the specials-family `ziv_calibration`
constants and the `agm_constants` `2/√π` and `ln(2π)` helpers are gated
by their consumer's feature, so a trig-only build (the profile
`pfloat-libm` uses) is free of dead code under `clippy -D warnings`.

## Consequences

- pfloat releases 1.1.0. The surface gap ADR-0032 recorded as "absent,
  deferred to the libm phase" closes; the per-function status table
  gains four correctly-rounded rows (cbrt, cot, sec, csc), and the
  tracked-function count moves from 63 to 67.
- `pfloat-libm` can now wrap the six kernels behind the `trig` feature it
  already depends on. The conversion API (ADR-0055) plus these kernels
  are the pfloat side of Phase 2.
- The reciprocal-trig pattern is now blessed by name with `tan` as the
  cited precedent, so a future reciprocal or composed kernel inherits it
  without reverse-engineering `tan`.
- `iroot_limbs` is a reusable integer k-th-root primitive; a future
  small-order `rootn` fast path can share it.
- The status table records cot/sec/csc as `oracle = "Arb"`,
  `oracle_independence = "independent"`. The correctness claim rests on
  Arb's wholly separate implementation; the cross-check baseline (the
  pf-hcz4-style sweep) for the new kernels is future work, so their rows
  carry no `[cross_check]` block yet.

## References

- ADR-0032: the direct-kernel policy these six implement.
- ADR-0055: the public f32/f64 conversion API; the other half of the
  pfloat side of Phase 2.
- ADR-0022, ADR-0038: the Ziv interval-test driver and its five-mode
  extension, reused by hypot, rootn, and the reciprocal trio.
- ADR-0031: the parse exponent cap, the precedent for bounding
  caller-supplied magnitudes (the rootn complexity discipline).
- DLMF §4.14 (reciprocal circular functions); IEEE 754-2019 §9.2.1
  (`hypot`, `rootn`, `cbrt` specifications).
- Plan: `~/.claude/plans/we-re-continuing-phase-2-modular-barto.md`.
