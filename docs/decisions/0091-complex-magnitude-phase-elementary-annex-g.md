# ADR-0091: complex magnitude, phase, and the elementary core with C99 Annex G branch cuts

- **Status**: accepted
- **Date**: 2026-06-08

## Context

Slice C4 of pfloat-complex lands the magnitude and phase surface (`abs`,
`norm_sqr`, `arg`, `to_polar`), the elementary core (`sqrt`, `exp`, `log`),
and the C99/C11 Annex G §G.5.1 complex-infinity refinement for `div` and `mul`
that C3 deferred. The hard, error-prone part is the branch-cut correctness:
Annex G fixes a semantic convention layered above rounding (which branch, which
signed zero), and a wrong row silently returns the wrong half of a branch cut.
Two failure classes recur. First, the C3 divide review found a real
sign-formula bug in a directed quotient enclosure (the naive endpoint selection
held for only one sign of the numerator); the same class threatens every
component enclosure here. Second, a result zero whose sign is fixed by an input
sign must carry that sign, and a zero formed by arithmetic (or stamped from the
rounding mode) loses it.

The branch-cut tables and the §G.5.1 recovery were derived from the Annex G
text and cross-checked, then each kernel was implemented and adversarially
verified by an independent re-derivation against the standard (the project's
verify-the-verdict discipline), not trusted from the first draft.

## Decision

### Rounding model: componentwise correct rounding

Complex numbers carry no total order, so a single complex directed rounding has
no meaning. The only coherent strong claim is componentwise correct rounding:
each of the real and imaginary parts is correctly rounded under its own real
rounding mode (MPC's model). Branch selection and signed-zero discrimination
are a documented Annex G convention on top of that, not a rounding guarantee.

### Magnitude and phase delegate to one correctly-rounded kernel each (class A)

`abs = hypot(re, im)`, `norm_sqr = mul_add_mul(re, re, im, im) = re² + im²` (one
fused rounding, ADR-0088, distinct from `abs²`), `arg = atan2(im, re)`. Each is
exactly one already-verified correctly-rounded scalar kernel, so it is
correctly rounded for free with no loop. The phase's entire negative-real-axis
branch cut and its `±0 → ±π/±0` signed-zero discrimination ride `atan2`'s IEEE
754-2019 §9.2.1 table in all five modes. The sealed `RealScalar` trait grows
(non-breaking) by exactly `hypot` and `atan2` for these; the other kernels do
not cross the trait boundary (see below).

### The elementary core uses a directed-pair enclosure, not a calibrated guard

`csqrt`, `cexp`, and `clog` have at least one component that is a composition of
transcendentals, not a single correctly-rounded kernel, so it needs a Ziv loop.
Every such component (class B) uses the C3 divide's directed-pair enclosure
(shared in `src/enclosure.rs`): bracket the component's true value with
`TowardNegative` / `TowardPositive` sub-evaluations at a growing working
precision (`GUARDS = [64, 128, 256, 512, 1024]`, cap five), round both ends to
the output precision under the target mode, and accept when they agree in value
**and** sign (the sign clause separates `±0`). `INEXACT` is computed from
whether the bracket collapsed (`bracket_is_exact`), never forced, so an exact
algebraic output reports `OK`.

The alternative, a per-function calibrated `error_guard` (pfloat's scalar-kernel
idiom), is rejected. A calibration constant is a recalled-tier fact whose
correctness rests on the sweep having hit the worst case (the directed-mode
saturation defect class, ADR-0080); the directed bracket carries no such number
and *is* the certificate. The enclosure also keeps one verified pattern across
the crate (`div`, `csqrt`, `cexp`, `clog`), and it handles the cancellation
regimes by guard-growth coupling rather than a half-width estimate. Because the
working precision exceeds any `FixedFloat<PREC>`, the kernels run in `BigFloat`
and the generic methods bridge through `RealScalar::to_big` / `from_big`, as
`div` does. The trait therefore needs no `sqrt` / `div` / `copysign` /
sign-predicate methods: those are inherent `BigFloat` methods the enclosures
call directly.

### Per-function construction

- **csqrt** (§G.6.4.2, principal branch `Re ≥ 0`, cut on the negative real
  axis, continuous from above). Kahan's robust form
  `w = sqrt((|x| + hypot(x, y))/2)` (which only ever adds, so it loses no bits
  to `|z| − x` cancellation), then `u = w`, `v = y/(2w)` for `x ≥ 0` and
  `v = copysign(w, y)`, `u = |y|/(2w)` for `x < 0`. The `y/(2w)` division is
  enclosed with sign-aware endpoint selection for both signs of `y` (the C3 bug
  class, handled explicitly per sign). The real-axis zeros (`y = ±0`) are
  stamped with `copysign(0, y)` directly, never routed through a division.
- **cexp** (§G.6.3.1, entire). Classified on the input real-part class. The
  finite interior encloses the two products `e^x cos y`, `e^x sin y` with a
  sign-aware product bracket; because `e^x > 0` the product sign follows the
  trig factor, so the endpoints are picked by the trig bracket's sign alone. A
  trig bracket straddling zero (`y` near `kπ/2`) does not converge and grows the
  guard (transient by Niven). The infinite-real-part rows with finite `y` take
  their component signs from `sign(cos y)`, `sign(sin y)`.
- **clog** (§G.6.3.2, cut on the negative real axis). `im = atan2(y, x)`
  (class A). `re = ln(hypot(x, y))` enclosed directly: `hypot` bracketed, then
  `ln` (increasing) applied to each end. The four poles, `∞ → +∞`,
  `NaN → NaN`, and `clog(1 + 0i) = +0` fall out of this composition (`ln(1)` is
  exact). Near `|z| = 1` the bracket straddles `ln(1) = 0` and the guard grows;
  the more aggressive `½·log1p((|x|−1)(|x|+1) + y²)` reformulation, which would
  push the convergence boundary much closer to the unit circle, is a deferred
  refinement, because the direct enclosure is correct and converges for `|z|`
  outside a band of width roughly `2^-1024` around 1 (the documented
  measure-zero cap caveat, shared with `div` and MPFR).

### The cross-function signed-zero rule

`div.rs::resolve`'s exact-zero short-circuit stamps the `±0` sign from the
rounding mode, which is correct for `z/z` cancellation but wrong for any zero
whose sign is fixed by an input (csqrt's real-axis imaginary part, cexp's
imaginary part at `y = ±0`, every Annex G directional zero). Such zeros are
stamped with explicit `copysign(0, input)` on the value path, never routed
through `resolve`. This is the one warning every kernel observes; the csqrt,
cexp, and clog adversarial verifications independently converged on it.

### §G.5.1 complex-infinity recovery for div and mul (full)

The naive componentwise formulas do not preserve a complex infinity (§G.3): a
zero divisor, an infinite dividend, or `(1 + 0i)·(∞ + ∞i)` collapses to `NaN`
where Annex G mandates an infinity. The full §G.5.1 recovery lands as exact
pre-dispatch in `src/specials.rs` (helper decomposition this crate's own,
derived from the algorithm, not copied), so the directed-pair divide and the
fused multiply only ever see finite operands:

- **D1** (any dividend / complex-zero divisor): directed infinity with the sign
  from the divisor real part `c` alone; a zero dividend part yields `∞·0 = NaN`
  for that part, and `0/0` falls out as `(NaN, NaN)`.
- **D2** (complex-infinite dividend / finite divisor) and **D3** (finite
  dividend / complex-infinite divisor): box the infinite operand
  (`∞ → ±1`, finite/NaN `→ ±0`) and scale the finite numerators by `∞` or `0`.
- **M1** (complex-infinity multiply recovery, applied only when the fused
  product already collapsed to `(NaN, NaN)`): the dual boxing.
- **M-OVF** (the finite-operand cross-product overflow branch) is named and
  intentionally absent: `BigFloat`'s i64 saturating exponent never overflows a
  finite product before the round (the same exemption `hypot` relies on,
  ADR-0032).

Per §G.5.1p5 and footnote 377 the values are normative and pinned; the raised
flags are not, so the code specifies values exactly and lets flags fall out
best-effort (the D1 `∞·0` path raises `INVALID`, not a mandated `DIV_BY_ZERO`).

### Feature gating, documented choices, and a flag nuance

The elementary surface is opt-in behind the `exp-log` and `trig` pass-through
features (mirroring pfloat and pfloat-ball): `hypot`/`exp`/`ln`/`log1p` and the
csqrt path are `exp-log`; `atan2`/`sin`/`cos`, cexp, and clog are `trig` (which
pulls `exp-log`). A bare `Complex<FixedFloat<PREC>>` keeps add/sub/mul/div
without the transcendental surface.

Where the standard leaves a sign genuinely unspecified (csqrt `−∞ + NaN·i`
imaginary sign; the `±∞` real parts of cexp's `x = +∞, y = ∞/NaN` rows; the
`±0` parts of cexp's `x = −∞, y = ∞/NaN` rows), the crate picks a deterministic
representative (`+∞` / `+0`) and the verification asserts magnitude, NaN-ness,
and the chosen sign with a comment that the standard permits either.

One flag nuance, surfaced by the cexp verification: cexp's `OVERFLOW` for a
finite `x` whose `e^x` saturates is essentially unreachable inside `BigFloat`
(i64 exponents saturate only near `x ≈ 10^19`), so the flag is material mainly
on the `from_big` bridge to a fixed-width scalar; the plumbing carries it
correctly through the exp directed pair if it ever fires.

## Consequences

- The magnitude, phase, and the imaginary part of `log` cost one verified
  kernel call each. The sealed trait stays minimal (two new methods), and the
  seal's "every component is a correctly-rounded pfloat scalar" story holds.
- One enclosure pattern (`src/enclosure.rs`) backs `div`, `csqrt`, `cexp`, and
  `clog`. The cost is recomputing each kernel twice per iteration and the
  measure-zero cap caveat on hard-to-round inputs, shared with `div`.
- The named failure modes for the disclosure: a wrong-branch result when a
  caller feeds an unsigned zero where the sign of zero was the only
  distinguishing information; and catastrophic cancellation in `ac − bd` /
  `ad + bc` and in `clog` near `|z| = 1`, on inputs no random-point sweep
  lands on.
- Deferred with a path to v1.x: `sin` / `cos` / `tan`, the hyperbolics, and
  inverse trig with their Annex G cuts; `pow` / `cis` / `from_polar`; and the
  `clog` `log1p` reformulation that tightens the `|z| ≈ 1` convergence band.
- Provenance: the csqrt / cexp / clog special-value tables were cross-checked
  against multiple independent transcriptions (cppreference, POSIX, the
  compiler-rt §G.5.1 reference for div/mul) because the ISO N1570 PDF was not
  machine-extractable in the deriving sessions. The C5 verification pins them
  against a paper N1570 §G.6 and an MPC/`acb` differential before the 1.0 cut.

## Related

- Plan: `plans/magical-skipping-lagoon.md` (C4-C6)
- Commits: `1c66c3b` (C4.1 magnitude/phase), `6ad990b` (C4.2 §G.5.1 recovery),
  `ec41069` (C4.3 csqrt), `5b4e26c` (C4.4 cexp), C4.5 clog (this slice)
- Other ADRs: builds on ADR-0088 (the fused two-product primitive),
  ADR-0089 (the sealed `RealScalar` trait), ADR-0090 (mul/div and the
  directed-pair enclosure); the verification posture is ADR-0092 (C5)
