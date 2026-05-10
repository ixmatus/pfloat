# ADR-0003: Dual API, `BigFloat` (dynamic) and `FixedFloat<const PREC: u32>` (const-generic)

- **Status**: accepted
- **Date**: 2026-05-10

## Context

Users of arbitrary-precision floats split into two profiles.

The dynamic-precision profile (MPFR's natural fit, `rug`'s shape)
chooses precision at runtime. Convergence sweeps walk through
increasing precision until a residual converges; libraries that take
user-specified precision through a configuration field need to hold
the value at whatever the user asked for. This profile demands a
precision field at runtime and a heap-allocated mantissa whose length
follows it.

The fixed-precision profile knows the precision at compile time.
Embedded users hard-code the precision their hardware budget allows;
hot loops at IEEE 754 binary64 (53 bits) want the optimizer to see
the precision as a constant; numerical code at binary128 (113 bits)
benefits from the same. This profile wants stack allocation, no
runtime branches on precision, and the ability to run on `no_std`
targets without an allocator.

A single API serving both profiles is possible. MPFR ships only the
dynamic shape; users wanting fixed precision pay for runtime
flexibility they do not use. Pure-Rust attempts so far (astro-float)
have followed MPFR's shape. Const generics make the fixed-precision
shape representable without touching the dynamic one.

## Decision

Ship two types in pfloat with the same operational surface.

```rust
pub struct BigFloat {
    class: Class,
    precision: u32, // bits
}
// mantissa stored as Vec<u64> inside Class::Normal

pub struct FixedFloat<const PREC: u32> {
    class: ClassFixed<PREC>,
}
// mantissa stored as [u64; ceil(PREC / 64)] inside ClassFixed::Normal
```

Both implement the same arithmetic, rounding, conversion, and
transcendental traits. The kernels share their core logic via a
private `Mantissa` trait that abstracts over `&[u64]` and
`&[u64; N]`. Users see only the concrete types.

Conversions are explicit:

- `BigFloat::from(FixedFloat<PREC>)` is exact and infallible.
- `FixedFloat::<PREC>::try_from(BigFloat)` rounds under a chosen
  mode, may set inexact / overflow / underflow, and threads the
  rounding mode through the call.

`big` and `fixed` are independent feature flags; either can ship
alone.

## Consequences

**Wins:**

- The fixed profile drops the runtime precision field and the heap
  allocation. Embedded code at known precision pays no per-value
  overhead for runtime flexibility.
- `FixedFloat<53>` enables a path to a correctly-rounded `f64`
  replacement at every rounding mode for callers that need IEEE
  rounding-mode control without hardware FPU support.
- The const-generic instantiation gives the optimizer compile-time
  visibility of mantissa length. Loops unroll cleanly; bounds checks
  fold; branch predictors stay sharp.
- `no_std`-without-`alloc` is a real, supported configuration:
  `--no-default-features --features=fixed`.
- `BigFloat` continues to serve the MPFR-shaped use case at full
  generality.

**Costs:**

- Test surface roughly doubles. Every operation is exercised against
  both types. Conformance corpora and differential lanes run twice.
- API documentation has to introduce both types and the conversion
  contract early. The dual surface is mildly harder to read than a
  single type would be.
- Const-generic precision arithmetic relies on Rust features that
  occasionally evolve. The MSRV (1.84) is the floor; if a newer
  feature is needed, the MSRV moves only when there is a concrete
  win to justify the bump.
- A small risk: callers who learn the API on `BigFloat` and then
  reach for `FixedFloat` may try to mutate precision at runtime.
  The type system rules that out at compile time; the failure is
  loud, not silent.

## Related

- ADR-0001 (limb representation, shared by both types)
- ADR-0002 (bit-level precision, shared by both types)
- ADR-0004 (storage strategy that distinguishes the two)
- ADR-0005 (`Class` and `ClassFixed`)
- DESIGN.md, "Type architecture" section.

## Update (2026-05-10)

The MSRV-on-stable stance in *Costs* is **superseded by [ADR-0011](0011-msrv-nightly-for-generic-const-exprs.md)**. pfloat now requires a date-pinned nightly toolchain so `FixedFloat<const PREC: u32>`'s mantissa storage can use `feature(generic_const_exprs)` directly. The dual-API design recorded above is otherwise unchanged: `BigFloat` and `FixedFloat<PREC>` ship together in 1.0 with the same operational surface, the same conversions, and the shared private `Mantissa` trait abstracting over their storage.
