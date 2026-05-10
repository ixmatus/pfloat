# ADR-0001: `u64` limb representation, sign-magnitude, top-bit-set normalization

- **Status**: proposed
- **Date**: 2026-05-10

## Context

A multi-precision floating-point library has to pick a layout for the
mantissa: limb width, limb order, normalization rule, and how the
sign is carried. The choice runs through every kernel and is awkward
to revisit later.

The candidates:

- **Limb width**: `u32`, `u64`, `u128`. `u64` is the default in the
  literature (Brent and Zimmermann, *Modern Computer Arithmetic*) and
  in MPFR. `u128` doubles the per-limb cost on most targets without
  doubling the per-bit throughput. `u32` doubles the limb count for
  the same precision.
- **Limb order**: little-endian (limb 0 is least significant) or
  big-endian. MPFR uses little-endian; so does astro-float; so does
  the bulk of the literature.
- **Normalization**: top bit of most-significant limb is `1` (no
  implicit bit), or top bit is `0` (one bit of headroom for additions
  before the carry pushes through).
- **Sign**: separate field, or fold into the mantissa as two's-complement.

astro-float uses signed `i64` limbs (two's-complement mantissa) to
sidestep the sign-tracking dance during subtraction. The cost is one
bit of precision per limb and a non-standard layout that diverges
from MPFR-shaped algorithms in the literature.

## Decision

- Limb type: `u64`, unsigned.
- Limb order: little-endian. `mantissa[0]` is the least significant
  64 bits; `mantissa[len - 1]` is the most significant.
- Normalization: the most significant bit of the most-significant
  limb is `1` for every normalized non-zero value. No implicit bit.
- Sign: separate field, lifted into the `Class` enum (see ADR-0005).
  Sign-magnitude representation throughout.

## Consequences

**Wins:**

- The `u64 × u64 → u128` width-doubling identity is a single
  primitive, available on every target via `wrapping_mul` and
  `widening_mul` (or the manual `(a as u128) * (b as u128)` form).
  Multiplication kernels read directly from MCA pseudocode.
- The `u128 / u64 → (u64, u64)` divrem primitive is the standard
  long-division building block; same story.
- Precision in bits maps cleanly to limb count via
  `ceil(precision / 64)`. ADR-0002 records the bit-level precision
  decision.
- Sign-magnitude matches IEEE 754 mental model and the spec language.
  No code path has to reason about whether `−0` is a different value
  than `+0` at the bit level; it is, by construction.

**Costs:**

- Subtraction has to detect which operand is larger and conditionally
  swap, instead of trusting the two's-complement addition unit. The
  cost is a few branches at the top of the kernel, paid once per
  operation.
- Targets where `u128` is software-emulated (some Cortex-M0+
  toolchains) pay a small constant for the `u64 × u64` step. The
  alternative (`u32` limbs) doubles the limb count and the
  per-operation overhead in the common case, so this remains the
  better tradeoff even there.
- The "no implicit bit" rule means the most-significant limb's top
  bit is always `1`, which is a redundant fact for the canonical
  form. The redundancy is the price of matching MPFR's algorithmic
  literature; an implicit-bit form would cost more in the kernel
  code than it saves in storage.

## Related

- ADR-0002 (bit-level precision)
- ADR-0004 (mantissa storage)
- ADR-0005 (`Class` enum, where the sign lives)
- DESIGN.md, "Numeric representation" section.
