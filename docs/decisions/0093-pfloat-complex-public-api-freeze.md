# ADR-0093: pfloat-complex 1.0 public API freeze

- **Status**: accepted
- **Date**: 2026-06-08

## Context

Slices C1 through C5 take pfloat-complex from 0.1.0 to 1.0.0. The surface is
the componentwise correctly-rounded complex cut: the fused two-product primitive
(ADR-0088), the sealed `RealScalar` trait (ADR-0089), `mul`/`div` (ADR-0090),
the magnitude/phase and elementary core (`csqrt`/`cexp`/`clog`) with C99 Annex G
branch cuts and the §G.5.1 infinity recovery (ADR-0091), and the five-lane
verification (ADR-0092). No new construction is part of the 1.0 cut. This ADR is
the API stability review that precedes the tag, the complex analog of ADR-0054
(scalar) and ADR-0086 (ball). Once 1.0 ships, the public surface is under semver:
additive changes only until 2.0.

The review found the surface freeze-ready and resolved the two API questions
ADR-0090 deferred (the merged-`Status` and single-`RoundingMode` calls). The
changes the ceremony still makes are the version bump `0.1.0` -> `1.0.0`, the
crate-root and README prose flip from the 0.x disclaimer to the 1.0 commitments,
the "How pfloat-complex is developed" disclosure, and the componentwise
rounding-posture doc (`docs/complex-rounding-status.md`).

## Decision

Freeze the public surface as it stands at this commit. The version bump to 1.0.0
lands in the same ceremony.

### The frozen surface, by feature

- **Always:** nothing standalone. Unlike the ball (whose `Mag` is a no_std
  alloc-free primitive), pfloat-complex has no public item below `big`: a
  `Complex<T>` is a pair of pfloat scalars, and the scalar engine needs `alloc`.
  A bare `--no-default-features` build compiles to an empty crate.
- **`big` (default):** `Complex<T>` and the sealed `RealScalar` trait, with the
  impl for `BigFloat`. The arithmetic surface (`new`, `re`, `im`, `is_nan`,
  `neg`, `conj`, `norm_sqr`, `add`, `sub`, `mul`, `div`) is available here,
  including the §G.5.1 complex-infinity recovery on `mul`/`div`.
- **`fixed`:** the `RealScalar` impl for `FixedFloat<PREC>`, so
  `Complex<FixedFloat<PREC>>`. `fixed` implies `big` (pfloat's `fixed`
  delegates through `big` today; the constraint is stated in the crate docs).
- **`exp-log`:** `Complex::abs` (= `hypot`) and `Complex::sqrt` (the Annex G
  §G.6.4.2 `csqrt`), plus the `RealScalar::hypot` method.
- **`trig`** (implies `exp-log`): `Complex::arg` (= `atan2`), `to_polar`,
  `Complex::exp` (§G.6.3.1 `cexp`), and `Complex::log` (§G.6.3.2 `clog`), plus
  the `RealScalar::atan2` method.

The arithmetic and elementary surface is exposed as inherent methods on
`Complex<T>`; the `complex` and `scalar` modules are public only as the home of
those impls and the `RealScalar` trait. `Complex`'s `re`/`im` fields are public.

### The componentwise rounding model is the frozen contract

Each operation rounds the real and imaginary parts each correctly under their own
real rounding mode (the MPC model; the only coherent strong rounding claim for a
type with no total order). Branch selection and signed-zero discrimination are a
documented C99 Annex G convention layered on top, not a rounding guarantee. This
is recorded per operation in `docs/complex-rounding-status.md`.

### Resolved: a single merged `Status`, not a per-component `ComplexStatus`

Every flag-producing operation returns `(value, Status)` where `Status` is the
single IEEE 754-2019 sticky-flag set, the OR-monoid merge of the two component
statuses (ADR-0007, the scalar/ball convention). The review keeps this for 1.0
rather than a per-component `ComplexStatus { re: Status, im: Status }`. A user who
needs to know *which* part was inexact can recompute each component through the
scalar API; the common need (did anything round, did anything signal) is the
merge. A per-component status is a strictly ADDITIVE v1.x extension (a new
accessor or a `*_components` variant), so deferring it costs no future
compatibility. `neg` and `conj` are exact sign-bit flips and carry no `Status`,
by design.

### Resolved: a single `RoundingMode`, not a per-component `ComplexRoundingMode`

Every rounding operation takes one `RoundingMode`, applied to both components.
This is the diagonal of MPC's `mpc_rnd_t` (which encodes a mode per part). The
review keeps the single-mode API for 1.0: it is the clean, complete, common
case, and per-component modes are a strictly ADDITIVE v1.x extension (a
`ComplexRoundingMode(RoundingMode, RoundingMode)` with `From<RoundingMode>`, or
`*_with_modes` method variants) that does not break the frozen signatures.
Recorded so the single mode is a deliberate scope choice, not an oversight.

### Public fields, no validity invariant (locked)

`Complex { re, im }` has public fields. A complex number carries no validity
invariant: any pair of real components, including NaN and infinity parts,
denotes a valid value (an infinite part is a complex infinity by Annex G §G.3),
so there is nothing for accessors to protect. The `re()`/`im()`/`is_nan()`
accessors are conveniences, not encapsulation. There is no `ComplexError` and no
fallible constructor: `new` is total.

### The sealed `RealScalar` trait can grow additively (locked)

`RealScalar` is sealed (ADR-0089): implemented only for `BigFloat` and
`FixedFloat<PREC>`, so no external crate can add an impl. Adding methods to it is
therefore NON-breaking (there are no external impls to break), which is how it
grew from the C1 arithmetic subset to the C4 `hypot`/`atan2` additions. The seal
is the v1.0 guarantee that "every `Complex` component is a verified,
correctly-rounded pfloat scalar" cannot be violated through this crate's surface.
The seal is scoped, not universal: a third party can still build
`num_complex::Complex<FixedFloat<P>>` outside this crate (pfloat ADR-0070), which
`RealScalar` does not and cannot prevent.

### Precision is `u32`, mirroring the scalar and ball crates

Component precision is `u32` (`RealScalar::precision() -> u32`), with no
`NonZeroU32`, for the same reason ADR-0086 records for the ball: a `Complex` over
a `u32`-precision scalar cannot enforce a stronger invariant than the type it
wraps. The `NonZeroU32` question stays a single, deferred, cross-crate 2.0
decision.

### Outside the 1.0 stability guarantee

The `differential-acb` feature (the per-release acb componentwise
certified-rounding lane, reached through the python-flint subprocess and never in
the shipped link graph) and the `cfg(kani)` harnesses (the `Status`-merge proofs)
are verification hooks, off in production builds, and are explicitly **not**
covered by the semver guarantee. The nightly toolchain pin (ADR-0011) applies:
the MSRV is the pinned nightly in `rust-toolchain.toml`, because the crate builds
pfloat, which needs `generic_const_exprs`.

### Deferred to v1.x, all strictly additive

`sin`/`cos`/`tan`, the hyperbolics, and inverse trig with their Annex G cuts;
`pow`/`cis`/`from_polar`; the `clog` `log1p` reformulation that tightens the
`|z| ≈ 1` convergence band; a per-component `Status`; per-component rounding
modes; and the `ComplexBall = Complex<Ball>` join. None of these breaks the
frozen surface.

## Consequences

- The surface is committed. From 1.0, changes to the items above are
  semver-significant: additive until 2.0. A user moving among the scalar, ball,
  and complex crates sees one method vocabulary (`(value, Status)`, `u32`
  precision, a single `RoundingMode` argument).
- The componentwise correct-rounding claim is backed by the enumerated Annex G
  tables, the exhaustive dispatch-totality enumeration, the algebraic identities,
  and the independent acb componentwise certified-rounding differential
  (ADR-0092). The branch-cut convention is the documented Annex G semantics, with
  the named failure mode (a wrong-branch result on an unsigned zero) recorded in
  the disclosure.
- The two ADR-0090 deferrals (merged `Status`, single `RoundingMode`) are closed
  for 1.0 as the simple form, with the richer per-component forms recorded as
  additive v1.x extensions.
- `differential-acb` and `cfg(kani)` consumers take no stability promise, by
  design. crates.io publication is deferred by choice (like pfloat `pf-4fi` and
  ball B0.5); the signed `pfloat-complex-v1.0.0` tag is the milestone.

## Related

- ADR-0054: the scalar crate's v1.0 public API freeze (the template).
- ADR-0086: the ball's v1.0 API freeze (the sibling this mirrors).
- ADR-0007: the `(value, Status)` return convention.
- ADR-0088 / ADR-0089 / ADR-0090 / ADR-0091: the surface this freezes (fused
  two-product, sealed `RealScalar`, mul/div, magnitude/phase/elementary core).
- ADR-0092: the verification posture backing the frozen claims.
- ADR-0011: the nightly toolchain pin.
- Plan: `plans/magical-skipping-lagoon.md` (C6).
