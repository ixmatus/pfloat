# ADR-0002: Bit-level precision granularity

- **Status**: proposed
- **Date**: 2026-05-10

## Context

Precision can be tracked at the bit level (1-bit minimum, MPFR's
choice) or at the limb level (64-bit-aligned, astro-float's choice).
The choice flows into every API that takes a precision argument and
into the rounding pipeline.

Bit-level granularity lets callers ask for IEEE-defined precisions
directly: 24 bits for binary32, 53 for binary64, 113 for binary128,
and any custom width. Limb-level granularity rounds up internally;
the caller asking for 53 bits gets 64 in storage and has to
post-process if exact-width semantics matter.

The cost of bit-level precision is a single bit-mask on the
most-significant limb's high bits and a corresponding mask in the
rounding pipeline. Once paid, the rest of the arithmetic runs at
limb granularity regardless.

## Decision

Precision is tracked at the bit level. The minimum precision is one
bit. The maximum is `u32::MAX` (about 4 billion bits, vastly
exceeding any practical workload).

The mantissa storage rounds up to whole limbs:
`mantissa.len() == ceil(precision / 64)`. The high bits of the
most-significant limb that exceed the working precision are zero in
canonical form; the rounding pipeline maintains that invariant.

The precision field is `u32`. ADR-0003 records how the field carries
into `BigFloat` (runtime field) versus `FixedFloat<const PREC: u32>`
(const generic).

## Consequences

**Wins:**

- IEEE-defined precisions are first-class. `BigFloat::with_precision(53)`
  matches `f64`'s rounding behavior under any rounding mode;
  `FixedFloat<53>` does the same at compile time.
- Differential testing against MPFR is direct: pass the same
  precision to both, get the same rounding under matching modes.
- Conformance corpora from Lefèvre–Muller (which are tabulated at
  specific bit widths) plug in unmodified.

**Costs:**

- The most-significant-limb mask is a small per-operation overhead
  (one `&` and one bit-count for the rounding-position calculation).
  Negligible against the limb-level multiplication cost.
- Some bit-level precisions land awkwardly within a limb: precision
  65 needs a two-limb mantissa with 63 unused bits in the top limb.
  The waste is a property of the user's choice, not the
  representation; the rounding remains correct.
- The const-generic `FixedFloat<PREC>` synthesizes its mantissa size
  via `[u64; ((PREC + 63) / 64) as usize]`. Stable Rust supports this
  through `generic_const_exprs`'s subset; if a workaround is needed,
  the layout is a `const fn` over `PREC`.

## Related

- ADR-0001 (limb representation)
- ADR-0003 (dual API, where precision lives)
- DESIGN.md, "Precision granularity" subsection.
