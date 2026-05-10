# ADR-0007: Rounding mode and exception flags

- **Status**: proposed
- **Date**: 2026-05-10

## Context

IEEE 754-2019 §4.3 defines five rounding modes: round to nearest with
ties to even (RNE, the default), round to nearest with ties away
from zero (RNA), round toward zero (RZ), round toward positive
infinity (RP), and round toward negative infinity (RM). It also
defines five sticky exception flags: inexact, overflow, underflow,
divide-by-zero, and invalid (§7).

A library has to pick how the rounding mode reaches each operation
and how the flags accumulate.

For rounding mode:

- **Function-call enum**, like MPFR. Pass a `RoundingMode` enum to
  every operation. Ergonomic for one-off calls; tedious for long
  expressions where every operation rounds the same way.
- **Typestate**, encoding the rounding mode in the type. The
  compiler refuses to mix rounding modes across an expression
  without an explicit conversion. Strong correctness guarantee;
  awkward at API boundaries.
- **Hardware-style global mode**, with `set_rounding_mode` /
  `get_rounding_mode`. C's `fenv.h` shape. Wrong for a software
  library that wants thread safety and reasoning about side
  effects.

For flags:

- **Thread-local sticky state**, like `fenv.h`. Free functions
  test, clear, and set. Ergonomic; works only with a thread-locality
  primitive (`std`).
- **Per-operation `(value, Status)` return**, like ferrodec. Status
  threads through every call explicitly. No global state; verbose at
  call sites.
- **Optional output reference**, like Rust's checked-arithmetic
  helpers. Caller passes `Option<&mut Status>`. Compromise between
  ergonomics and explicitness.

## Decision

For rounding mode: ship both a function-call enum and a typestate
wrapper.

```rust
pub enum RoundingMode {
    NearestEven,    // RNE (IEEE default)
    NearestAway,    // RNA
    TowardZero,     // RZ
    TowardPositive, // RP
    TowardNegative, // RM
}

pub trait RoundMode {
    const MODE: RoundingMode;
}
pub struct Rne; pub struct Rna; pub struct Rz; pub struct Rp; pub struct Rm;
// each with `impl RoundMode for ...`

pub struct Rounded<M: RoundMode, T>(pub T, PhantomData<M>);
```

The enum is the primary interface. The typestate wrapper is opt-in
for callers who want compile-time discipline across long expressions.
The wrapper threads the mode through arithmetic via `impl Add` /
`impl Mul` etc., matching ferrodec's `ops` feature shape.

For flags: gate by feature flag.

- `std`: thread-local `Cell<Status>`. Free functions
  `pfloat::flags::test`, `pfloat::flags::clear`,
  `pfloat::flags::set`. Each operation updates the thread-local
  before returning.
- `no_std` (without `std`): no global state. Operations take an
  explicit `Status` parameter:

  ```rust
  impl BigFloat {
      pub fn add_with_flags(&self, other: &Self, mode: RoundingMode, flags: &mut Status) -> Self;
      pub fn add(&self, other: &Self, mode: RoundingMode) -> (Self, Status);
  }
  ```

  The `&mut Status` form lets a caller accumulate flags across a
  long expression without per-operation `Status` returns. The
  `(Self, Status)` form is the simpler shape.

The two forms coexist; users on `std` who want explicit flags use
the `_with_flags` variants too.

## Consequences

**Wins:**

- The `std` user gets the IEEE-conventional sticky-flag mental model
  with no API noise.
- The `no_std` user pays no hidden cost: there is no thread-local
  primitive in `core`, so the explicit-status form is the honest
  shape for that target.
- The typestate wrapper is available to callers who want
  compile-time guarantees about rounding-mode propagation; they
  pay only for the wrapper they construct.
- Match between user mental model and IEEE 754 spec language is
  direct: "I called `add` under RNE; the inexact flag is now set
  because the result was not exactly representable."

**Costs:**

- Three calling conventions per operation (`add`, `add_with_flags`,
  `Rounded::add`). The combinatorics on the public API surface are
  managed by code generation through a macro that produces the
  three forms from one canonical kernel.
- Differential testing has to compare flag behavior as well as
  values. MPFR's flag conventions are close to IEEE 754 but not
  identical (the `mpfr_underflow` flag is set under different
  conditions than IEEE underflow). Translation lives in the
  `differential-mpfr` test module.
- Documentation has to make the "how do I get a flag" answer obvious.
  The README walkthroughs cover both shapes.

## Related

- ADR-0008 (differential testing oracle, where flag-translation
  lives)
- DESIGN.md, "Rounding modes and exception flags" subsection.
- ferrodec ADR-0002 (per-op status), the model for the no_std
  shape.
