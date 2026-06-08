# ADR-0089: pfloat-complex and its independent sealed RealScalar trait

- **Status**: accepted
- **Date**: 2026-06-08

## Context

Slice C2 creates the `pfloat-complex` crate and its first surface:
`Complex<T>` and componentwise additive arithmetic. `Complex<T>` must be
generic over a real scalar, and the question is what constrains `T`. Two
decisions need recording: how the scalar bound is shaped (and why it is
sealed), and where it lives, given that `pfloat-ball` already defines a
`RealScalar` trait of the same name (ADR-0075).

## Decision

### A sealed `RealScalar`, defined in pfloat-complex

`Complex<T>` is bounded by a sealed `RealScalar` trait, implemented only for
`pfloat::BigFloat` and `pfloat::FixedFloat<PREC>`. Sealing (a private
supertrait an external crate cannot name) makes "each component is a
correctly-rounded pfloat scalar" a fact this crate's surface cannot be made
to break: a `Complex` built through pfloat-complex can never be instantiated
over an unverified or wrongly-rounded scalar type. This is
illegal-states-unrepresentable applied to the type parameter, the failure
mode `num-complex` leaves open (its `Complex<T>` admits any `T`, so a
`Complex<UnverifiedFloat>` whose magnitude is not correctly rounded
typechecks).

The seal is **scoped, not universal**, and the docs say so rather than
overclaiming. Because Phase 3 shipped `num_traits::Num` for
`FixedFloat<PREC>` (pfloat ADR-0070), a third party can still build a
`num_complex::Complex<FixedFloat<P>>` outside this crate. `RealScalar`
closes *pfloat-complex's own* inhabitant set, not the universe of generic
numeric code. Stating this honestly is the disclosure discipline applied to
a type-level guarantee.

### Defined in pfloat-complex, independent of the ball's `RealScalar`

The trait is defined in pfloat-complex, not imported from pfloat-ball, even
though the two are near-identical in shape. The reason is the DAG:
`pfloat-complex` and `pfloat-ball` are sibling Phase 4 stars with no
dependency edge between them (both descend from `pfloat`; only the eventual
`ComplexBall = Complex<Ball>` join consumes both). Reusing the ball's trait
would couple complex to the ball, inventing an edge the roadmap does not
draw. The duplication (two small sealed traits, one per consumer) is the
deliberate cost of the roadmap's rule: shared traits are extracted from
concrete crates *after* they exist, never pre-designed (the abstract-algebra
trait-graveyard warning). If a shared `RealScalar` is ever worth extracting,
it happens once there are three concrete consumers to validate its shape,
not now.

The trait presents only the subset C2 needs (`add`, `sub`, `negated`, and
the predicates the tests use), each delegating to the inherent pfloat kernel
of the same name via UFCS (`BigFloat::add(self, ..)`, which resolves to the
inherent method, not the trait method, so it does not recurse; a test guards
this). Because the trait is sealed, later slices extend it (the
`mul_add_mul` / `mul_sub_mul` forms for complex multiply and divide,
`hypot`/`atan2` for magnitude and phase, the elementary kernels) with no
breaking change.

### Componentwise rounding; public fields

Arithmetic is **componentwise correctly rounded**: `add`/`sub` round the
real and imaginary parts each under their own real rounding mode and OR the
two component statuses. This is MPC's model and the only coherent strong
rounding claim for complex numbers, which carry no total order, so a single
complex directed rounding has no meaning (ADR-0091 will record the same for
the elementary functions and their branch cuts). `neg` and `conj` are exact
sign-bit flips, so they take no mode and return no `Status`.

`Complex<T>` has **public `re` / `im` fields**. A complex number carries no
validity invariant: any pair of real components, NaN and infinity parts
included, denotes a valid value, so there is nothing for private fields and
accessors to protect. This differs from `Ball` (private fields, because a
ball has the radius-non-negativity and finite-midpoint invariants) and
matches `num-complex`.

## Consequences

- `pfloat-complex` depends only on `pfloat`; it does not depend on
  `pfloat-ball`, keeping the two stars independent as the DAG intends.
- The correctly-rounded-component guarantee is a type fact for code that
  builds `Complex` through this crate, with the `num-complex` escape hatch
  documented rather than papered over.
- A small, deliberate duplication (two `RealScalar` traits) is accepted now;
  extraction waits for the third concrete consumer.
- The crate ships a CI lane from its first commit (the pf-xyaq lesson,
  ADR-0053, applied up front), so the member never compiles in no per-push
  job the way the ball briefly did before its 1.0 ceremony.

## Related

- ADR-0075: pfloat-ball's `RealScalar`, the sibling trait this deliberately
  does not reuse.
- ADR-0070: `num_traits::Num` for `FixedFloat`, the escape hatch that scopes
  the seal.
- ADR-0088: the fused two-product primitive the complex multiply/divide
  slice (C3) will reach through an extended `RealScalar`.
- ADR-0053: the per-push CI coverage discipline the new member's lane
  applies.
- Plan: `~/.claude/plans/plan-tower-expansion-scope-goofy-raven.md` (slice
  C2).
