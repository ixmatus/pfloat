# ADR-0119: make FixedFloat<0> uninstantiable at compile time

- **Status**: accepted
- **Date**: 2026-07-03

## Context

`FixedFloat<const PREC: u32>` (ADR-0003, ADR-0004) is the const-generic
counterpart to `BigFloat`, with a stack-allocated
`[u64; limbs_for(PREC)]` mantissa. `limbs_for(0)` returns `0`, so
`[u64; 0]` is a well-formed storage array and the type's only
where-clause bound, `[(); limbs_for(PREC)]:`, is satisfied for
`PREC == 0`. Nothing else rejected the value, so `FixedFloat<0>` was
instantiable (defect pf-sn3n, epic pf-8iji review remediation R4).

`FixedFloat<0>` is a nonsensical illegal state. A fixed-precision float
with zero significand bits can represent no finite non-zero value.
Worse, three concrete harms followed:

1. `FixedFloat::<0>::zero()` and the other five `const` special-value
   constructors built a `Self { class }` directly and handed back a
   usable value.
2. `to_big()` on such a value produced a `BigFloat` at precision 0,
   violating `BigFloat`'s own `precision >= 1` invariant (that type
   returns `BuildError::PrecisionZero` from every constructor).
3. `try_from_i64_round` and `try_from_big_round` reached a
   `.expect("PREC >= 1 by const-generic bound")` that named a bound
   which did not exist, so the panic message was a lie: the map claimed
   an invariant the territory did not hold.

Parnell's standing preference is illegal states unrepresentable, with
compile-time and type-level guarantees favored over runtime checks.
The constraint that shaped the fix: the public surface is frozen at
v1.0 (ADR-0054), the `const` constructors cannot be made fallible
(they return `Self`, not `Result`), and the crate leans on
`feature(generic_const_exprs)` (ADR-0011), which does not propagate
implied bounds. A bound added to the `FixedFloat` struct itself must be
restated at every one of the roughly 70 generic `impl` sites across the
arithmetic, transcendental, and special-function kernels, or the
library does not compile.

## Decision

Gate the value-birthing constructors with a **method-level** const-generic
bound rather than gating the struct. A new internal guard function

```rust
pub const fn require_nonzero_precision(precision: u32) -> usize {
    assert!(precision >= 1, "FixedFloat<0> is not a valid type: ...");
    0
}
```

returns `0` for every `precision >= 1` (so it sits in a zero-length
array bound without enlarging storage) and fails const-evaluation for
`precision == 0`. Because the bound is only ever evaluated inside a
const-generic where-clause, the failure surfaces as a compile error,
not a runtime panic. Placement follows the construction boundary:

- **Compile error at `PREC == 0`.** The six direct constructors that
  build storage without routing through `BigFloat` (`zero`, `neg_zero`,
  `infinity`, `neg_infinity`, `nan`, `signaling_nan`), plus the two
  rounding constructors (`try_from_i64_round`, `try_from_big_round`,
  which otherwise panicked), each carry
  `where [(); require_nonzero_precision(PREC)]:`. The `num-traits`
  construction traits that call these (`Zero`, `One`, `Num`, `Signed`,
  `FromPrimitive`, `NumCast`, `Inv`) inherit the bound at the `impl`
  level, so `FixedFloat<0>` implements none of them.
- **Typed error at `PREC == 0`, no gate.** `try_from_i64_exact`,
  `try_from_big_exact`, `parse_str`, and `serde` deserialization stay
  ungated. `try_from_big_exact` is the back-conversion primitive that
  all 55 arithmetic and transcendental kernels delegate their result
  through; gating it would force the struct-level cascade. None of
  these can produce a `FixedFloat<0>` regardless, because `BigFloat`
  rejects precision 0, so they already return a typed error (`BuildError`
  / `ParseError` / a serde error) and never a value or a panic.

The method-level bound gates only the method and its callers; the
ungated operations (`to_big`, `add`, `partial_cmp`, and the rest) keep
compiling for every `PREC`, so no cascade reaches the kernels. The two
now-truthful `.expect` messages cite this ADR.

Alternatives rejected:

- **Struct-level bound (strongest: makes the type unnameable).**
  Rejected because `generic_const_exprs` has no implied bounds, so the
  bound must be restated at roughly 70 `impl` where-clauses plus their
  imports across every kernel file. That churn is disproportionate to a
  narrow fix and inflates the review surface across unrelated code.
- **Pure runtime story.** Rejected because it cannot satisfy the hard
  requirement that a usable `FixedFloat<0>` be gone: the `const`
  constructors return `Self` and cannot be made fallible, so a runtime
  approach would still hand back a `FixedFloat<0>` value from `zero()`.

## Consequences

- No value of `FixedFloat<0>` can be born in safe code. The value-side
  illegal state is gone: direct construction is a compile error, the
  `num-traits` construction surface is a compile error, and every
  `BigFloat`-routed path returns a typed error.
- The diagnostic is honest and legible. In a consumer that enables
  `generic_const_exprs`, `FixedFloat::<0>::zero()` reports
  `error[E0080]: evaluation panicked: FixedFloat<0> is not a valid
  type: a fixed-precision float needs at least one significand bit
  (PREC >= 1)`, pointing at the constructor's bound.
- The type name stays inert. `Option<FixedFloat<0>>` is a well-formed
  type and may appear in a signature; it is simply uninhabited through
  the public API. This is weaker than the struct-level bound (which
  would make the name itself ill-formed) but achieves the operative
  goal (no usable value) at a fraction of the churn.
- Zero public-API growth. `require_nonzero_precision` is `pub` inside
  the private `mantissa` module (matching `limbs_for`, to satisfy the
  `private_interfaces` lint) but is not re-exported, so it adds nothing
  to the frozen v1.0 surface.
- `ClassFixed<0>` remains a nameable, constructable but inert residual.
  Gating it is not feasible: its bound feeds `FixedFloat`'s `class`
  field, which would re-trigger the struct-level cascade. A standalone
  `ClassFixed<0>` cannot reach `FixedFloat` (the `class` field is
  `pub(crate)`) and has no methods, so it is harmless.
- The change is confined to `src/fixed.rs`, `src/mantissa.rs`, and the
  seven `num-traits` impls in `src/num_traits_impls.rs`; the kernels
  are untouched.

### Inversion

The three most plausible ways this decision is wrong or harmful:

1. **A construction path was missed, and a `FixedFloat<0>` still leaks.**
   The audit traced every constructor: the six `const` specials and the
   two rounding constructors are gated; `try_from_i64_exact`,
   `try_from_big_exact`, `parse_str`, and `serde` deserialization all
   route through `BigFloat`, which returns `BuildError::PrecisionZero`,
   so they cannot build a value; the private `from_big_at_same_precision`
   is only reachable with a precision-0 `BigFloat`, which cannot exist.
   The `num-traits` closure was walked until the compiler reported no
   further `unconstrained generic constant`. Residual risk: an `unsafe`
   `mem::zeroed::<FixedFloat<0>>()` sidesteps the safe API, which is out
   of scope for a safe-code invariant and is the accepted boundary.
2. **The consumer-side `generic_const_exprs` requirement reads as a
   regression.** A downstream crate that names `FixedFloat<PREC>` with a
   concrete `PREC` must itself enable `#![feature(generic_const_exprs)]`,
   or it hits an unrelated `E0275` well-formedness overflow. This is
   pre-existing (it predates this ADR and is inherent to the const-generic
   storage of ADR-0011); it is documented on the regression test, not
   introduced here.
3. **Method-level bounds spread further than expected on the next
   change.** Any future method that constructs `±1` or a special value
   through a gated constructor will inherit the bound and, if it is a
   trait impl, restrict that trait for `FixedFloat<0>`. This is the
   intended shape (the construction boundary carries the invariant), but
   a contributor who adds an ungated sibling that calls a gated one will
   see the bound propagate; the fix is to add the same one-line bound,
   which the pattern in `src/fixed.rs` and `src/num_traits_impls.rs`
   demonstrates.

## Related

- Defect: pf-sn3n (epic pf-8iji, review remediation R4).
- Verification: `cargo test -p pfloat --features=fixed` (regression
  suite `tests/regression_pf8iji_fixed.rs`, 4 tests); full-union clippy
  at `-D warnings`; `cargo build --no-default-features --features=fixed`;
  `cargo fmt --check`. The negative case (`FixedFloat::<0>::zero()`
  yields `E0080`) was verified manually in a consumer crate.
- Other ADRs: ADR-0011 (`generic_const_exprs` for const-generic
  storage), ADR-0003 / ADR-0004 (`FixedFloat` design), ADR-0054 (v1.0
  API freeze).
