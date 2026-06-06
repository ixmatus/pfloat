# ADR-0070: num-traits impls scoped to FixedFloat

- **Status**: accepted
- **Date**: 2026-06-06

## Context

Roadmap Phase 3 calls for `num-traits` interop behind a feature, so
downstream generic numeric code (`fn f<T: Num>(...)`) can operate on
pfloat values.

Two facts shape the design. First, `num-traits`' constructors carry no
precision argument: `Zero::zero`, `One::one`, `Num::from_str_radix`,
`FromPrimitive::from_i64`. For `FixedFloat<PREC>` the precision is the
type parameter, so these are well defined. For `BigFloat` the precision
is a runtime field, so any impl would have to invent a hidden default
precision, exactly the kind of silent, surprising behavior this crate
avoids. Second, `num_traits::Num` requires `core::ops::Rem`, which only
became available with the remainder kernel (ADR-0069).

## Decision

Implement the num-traits stack for `FixedFloat<PREC>` only:
`Zero`, `One`, `Num`, `Signed`, `FromPrimitive`, `ToPrimitive`,
`NumCast`, and `Inv`. `BigFloat` is deliberately excluded; its dynamic
precision does not fit the no-precision-argument model, and dynamic
precision code uses `BigFloat`'s explicit precision-carrying methods
instead.

Every operation runs at the type's precision `PREC`, so there is no
hidden default anywhere. Specifics:

- `Zero::is_zero` is overridden to use the predicate rather than the
  default `*self == Self::zero()`, which would miss negative zero under
  the derived structural equality. `One::is_one` keeps its default: at
  a fixed precision `1.0` has a single canonical representation, so the
  derived equality is value equality for it.
- `Num::from_str_radix` supports radix 10 (parsing at `PREC` through
  `FixedFloat::parse_str`); other radixes return
  `RadixParseError::UnsupportedRadix` rather than being mishandled.
- The arithmetic the `Num` bound needs reuses the `core::ops`
  overloads from the `ops` feature (round to `PREC` under
  `NearestEven`, `Status` discarded). `Inv` is `one() / self`.
- `FromPrimitive` is lossless where the precision allows: `from_i64`
  and `from_u64` build the integer at `PREC` (the `u64 > i64::MAX` case
  splits as `(n - 2^63) + 2^63`); `from_f64` and `from_f32` widen the
  primitive exactly then round to `PREC`.
- `ToPrimitive` and `NumCast` route through `f64`, so a value with more
  than 53 significant bits loses precision converting to an integer or
  float. That is the standard lossy contract of those traits and is
  documented at the impl.

`Float` and `Real` are not implemented: pfloat carries no fixed
associated constants for them, and `PartialOrd` is intentionally absent
because comparison returns `(Option<Ordering>, Status)` per IEEE
754-2019 §5.11.

The feature is `num-traits = ["dep:num-traits", "fixed", "ops"]` and is
added to the CI test-matrix and clippy unions.

## Consequences

Generic numeric code over `T: Num`, `T: Signed`, and the conversion
traits works for `FixedFloat<PREC>` at any precision, verified in
`tests/num_traits_integration.rs` (generic sum and polynomial, zero and
one including negative zero, the signed predicates, primitive
round-trips, `NumCast`, `Inv`, `from_str_radix` plus its radix error,
and a second precision). The default build keeps its zero runtime
dependencies; `num-traits` enters the graph only when the feature is
on, and `default-features = false` keeps it no_std.

Excluding `BigFloat` means a caller holding a `BigFloat` cannot pass it
to generic `T: Num` code; that is the honest consequence of the
dynamic-precision mismatch, and the alternative (a hidden default
precision for `zero`/`one`/`from_str_radix`) was declined.

`ToPrimitive`'s `f64` routing is lossy beyond 53 bits; this is inherent
to converting an arbitrary-precision value to a 64-bit target and is
documented. A new public type `RadixParseError` is added for the
`from_str_radix` error. The new optional dependency `num-traits` is the
canonical ecosystem crate, justified by the interop the feature
provides.

## Related

- Plan: `~/.claude/plans/melodic-knitting-ripple.md` (Phase 3, slice C3c)
- Issue: `pf-a4jh`
- Other ADRs: depends on ADR-0069 (the remainder kernel that unblocks
  `core::ops::Rem` and thus `Num`); builds on ADR-0054 (v1.0 API
  freeze) and the `ops` feature's operator overloads
