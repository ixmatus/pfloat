# ADR-0075: the sealed `RealScalar` trait

- **Status**: accepted
- **Date**: 2026-06-06

## Context

`Ball<T>` is generic over its midpoint scalar `T`. The whole soundness
story of a ball rests on one fact: the midpoint is computed by a
correctly-rounded pfloat kernel and the radius bounds that kernel's
rounding error. If `T` could be any numeric type, a third party could
instantiate `Ball<SomeWrongScalar>` whose midpoint is not correctly
rounded, and every enclosure built on it would be quietly invalid — the
same failure mode num-complex exposes when `Complex<T>` is allowed over
a `T: num_traits::Num` whose magnitude is not correctly rounded.

So `T` must be constrained to exactly the verified pfloat scalar types,
and the constraint must be one a downstream crate cannot widen.

## Decision

Introduce `RealScalar` in `pfloat-ball` (not in pfloat), **sealed** by a
private supertrait, and implement it only for `pfloat::BigFloat` and
`pfloat::FixedFloat<PREC>`. A third party cannot name the private
`Sealed` supertrait, so it cannot add an impl: the inhabitant set of
`RealScalar` is closed to the two correctly-rounded pfloat scalars.
`Ball<T: RealScalar>` therefore can never be instantiated over an
unverified midpoint type. This is illegal-states-unrepresentable applied
to a type parameter.

The trait surface is the minimum a ball kernel needs: precision and
classification accessors, sign operations, IEEE partial comparison, the
five directed arithmetic operations (`add`/`sub`/`mul`/`div`/`sqrt`,
each correctly rounded under a supplied mode), the slice-1/2 primitives
(`scale_by_pow2`, `next_up`/`next_down`/`ulp`), and two `Mag` bridges
(`magnitude_to_mag` narrows `|self|` up to a radius; `radius_to_scalar`
widens a radius back to a scalar, rounded up). Each impl delegates to the
type's inherent method; inherent-method resolution wins over the
same-named trait method, so the delegation does not recurse (a generic
round-trip test guards this). Every pfloat method that is fallible only
on a zero precision is presented in its always-valid form, because a
`RealScalar` value always has precision `≥ 1`.

Choosing `RealScalar` over a `num_traits::Num` bound is the
correctness-as-a-type-fact call: `Num` would reproduce num-complex's
compromise. The seal is honestly **scoped, not universal**: Phase 3
already shipped `num_traits::Num` for `FixedFloat<PREC>` (ADR-0070), so a
third party can still write `Complex<FixedFloat<P>>` through num-complex
outside this crate. `RealScalar` closes *pfloat-ball's own* surface, and
the doc comment says so rather than over-claiming a universal
prohibition.

Implementing a generic `impl ... for FixedFloat<PREC>` from outside
pfloat requires spelling the `where [(); limbs_for(PREC)]:` bound that
`FixedFloat<PREC>` carries (ADR-0011). `limbs_for` lived in a private
module, so pfloat now re-exports it (`pub use mantissa::limbs_for`) — a
purely additive change in service of downstream trait impls.

## Consequences

- A `Ball` is structurally incapable of wrapping a non-pfloat or
  wrongly-rounded scalar. The central promise (the midpoint is correctly
  rounded) is enforced by the compiler, not by documentation.
- The trait can grow as later slices need it (it is sealed, so adding a
  method is not a breaking change for any external implementor — there
  are none). Slice 8 adds the feature-gated elementary-function methods
  this way.
- `pfloat::limbs_for` is now public. It is a small, stable, total
  `const fn`; exposing it lets any downstream crate write generic
  `FixedFloat<PREC>` impls, not just pfloat-ball.
- The delegation-without-recursion pattern (inherent resolution beats the
  trait method) is load-bearing and tested; a future rename that broke it
  would stack-overflow, which the round-trip test catches.

## Related

- Plan: `~/.claude/plans/melodic-knitting-ripple.md` (slice 4); design
  reference `~/.claude/plans/pfloat-phase4-scoping.md` (decision 5).
- Beads: `pf-icgj.4` (under epic `pf-icgj`).
- Other ADRs: builds on `Mag` (ADR-0074) and the slice-1/2 primitives
  (ADR-0072, ADR-0073); the `num_traits::Num for FixedFloat` it declines
  to reuse is ADR-0070; the `FixedFloat<PREC>` const-generic bound is
  ADR-0011. `Ball<T>` (ADR-0076) consumes this trait.
