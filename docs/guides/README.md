# Guides

Task oriented and learning oriented documentation for the pfloat family. The
reference material lives elsewhere: `docs/references.md` for the standards and
sources, `docs/algorithms.md` for the algorithm reading guide, and
`docs/decisions/` for the architecture decision records. These guides teach you
to use the crates.

Each guide names its type in the first line. A tutorial builds understanding and
holds your hand through a worked example; a how to gets one task done and
assumes you know why you want it. Every code block is mirrored by a compiled,
self checking program under the relevant crate's `examples/`, so the code in a
guide is code that runs.

## Available

- [01: your first correctly rounded value to N digits](01-first-value.md)
  (tutorial). Compute the square root of two to fifty correct digits; learn the
  precision in bits versus digits on output model and the `(value, Status)`
  shape. Companion: `examples/first_value.rs`.

## Planned

The set below is proposed; order is by value. The first three establish the
concepts the rest lean on.

- 02: read and act on the status flags (tutorial). `INEXACT`, `OVERFLOW`,
  `UNDERFLOW`, `INVALID`, `DIV_BY_ZERO`, and the thread local accumulator.
- 03: choosing a rounding mode, and why directed modes matter (tutorial). All
  five modes; a one ULP split under the two directed nearest modes.
- 04: rigorous error bounds with pfloat-ball (tutorial). Build a `Ball`, read
  the certified enclosure, see the true result guaranteed inside.
- 05: replace a C libm call with a correctly rounded pure Rust one (how to). Swap
  a `std` or `libm` call for `pfloat-libm`; select a mode.
- 06: reach a target accuracy without guessing a precision (how to).
  `refine_to_accuracy` grows the working precision to a certified accuracy.
- 07: test your own numerical code against a trusted oracle (how to). Use a
  `Ball` or a directed pair as a reference enclosure for an `f64` routine.
- 08: run pfloat in a no_std or embedded build (how to). `FixedFloat<PREC>`,
  feature flags, and why transcendentals still need `alloc`.
- 09: why correct rounding and interval soundness are hard (tutorial, an
  explanation). The table maker's dilemma, why a midpoint is not a bound, why a
  radius rounds outward.
