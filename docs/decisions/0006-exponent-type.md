# ADR-0006: `i64` exponent

- **Status**: proposed
- **Date**: 2026-05-10

## Context

The exponent in a floating-point representation has to cover the
range of the most extreme values the library will ever produce.
Practical workloads stay within `±10^4932` (binary128's range).
Algorithmic intermediates (especially in special functions) can
push higher: factorial-of-large-argument and Bessel-function
recurrences regularly produce exponents in the thousands or tens
of thousands before reflection or reciprocation pulls the result
back.

The candidates:

- `i32`: range `±2^31`, comfortably ≈ `±10^646,000,000`. Halves the
  exponent storage compared to `i64`. Sufficient for any realistic
  use.
- `i64`: range `±2^63`, vastly exceeding any realistic use.
  Standard in MPFR and astro-float.
- A bignum exponent (arbitrary-precision integer for the exponent
  itself). Theoretically clean. In practice every operation pays
  per-step overhead for a feature no caller exercises.

## Decision

Use `i64` for the exponent in both `BigFloat` and `FixedFloat<PREC>`.
Match MPFR's choice.

The exponent describes the position of the most-significant bit of
the mantissa relative to the binary point: a normalized non-zero
value `v` has

```
v = sign × mantissa × 2^(exponent - precision + 1)
```

where `mantissa` is interpreted as an unsigned integer of
`precision` bits with the top bit set.

## Consequences

**Wins:**

- Every operation on the exponent (addition during multiplication,
  subtraction during division, comparison during alignment) is a
  single `i64` operation with no overflow drama at realistic
  precisions.
- Differential testing against MPFR is direct. No translation cost
  at the exponent boundary.
- The exponent fits in a register on every 64-bit target, and on
  32-bit targets the cost is two registers, the same as MPFR.

**Costs:**

- For `FixedFloat<53>` (single-limb mantissa), the exponent is
  larger than the mantissa it describes. The proportion is silly
  but harmless; `FixedFloat<53>` is not optimizing for storage
  density relative to `f64` (use `f64` if you want that), it is
  optimizing for correctly-rounded behavior at every IEEE rounding
  mode. The exponent overhead is paid once per value and disappears
  in any realistic working set.
- Saturating-arithmetic checks on the exponent (during overflow /
  underflow detection) operate on `i64`. The implementation has to
  use `checked_add` and `checked_sub` rather than relying on
  hardware overflow flags. Standard practice; not novel.

## Related

- ADR-0001 (limb representation)
- ADR-0005 (`Class` enum, where the exponent rides)
- DESIGN.md, "Exponent" subsection.
