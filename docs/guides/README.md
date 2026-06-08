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

Tutorials, in reading order.

- [01: your first correctly rounded value to N digits](01-first-value.md). The
  precision in bits versus digits on output model and the `(value, Status)`
  shape. Companion: `examples/first_value.rs`.
- [02: read and act on the status flags](02-status-flags.md). `INEXACT`,
  `OVERFLOW`, `UNDERFLOW`, `INVALID`, `DIV_BY_ZERO`, and the thread local
  accumulator. Companion: `examples/status_flags.rs`.
- [03: choosing a rounding mode, and why directed modes matter](03-rounding-modes.md).
  All five modes and the one ULP directed bracket. Companion:
  `examples/rounding_modes.rs`.
- [04: rigorous error bounds with pfloat-ball](04-error-bounds.md). Build a
  `Ball`, read the certified enclosure, see the true result guaranteed inside.
  Companion: `pfloat-ball/examples/error_bounds.rs`.
- [09: why correct rounding and interval soundness are hard](09-why-hard.md). The
  table maker's dilemma, why a midpoint is not a bound, why a radius rounds
  outward. Companion: `examples/why_hard.rs`.

How-tos, by task.

- [05: replace a C libm call with a correctly rounded pure Rust one](05-replace-a-c-libm-call.md).
  Companion: `pfloat-libm/examples/correctly_rounded_exp.rs`.
- [06: reach a target accuracy without guessing a precision](06-target-accuracy.md).
  `refine_to_accuracy` grows the working precision to a certified accuracy.
  Companion: `pfloat-ball/examples/target_accuracy.rs`.
- [07: test your own numerical code against a trusted oracle](07-oracle-test.md).
  A `Ball` as a proven pass or fail bracket for an `f64` routine. Companion:
  `pfloat-ball/examples/oracle_test.rs`.
- [08: run pfloat in a no_std or embedded build](08-no-std-embedded.md).
  `FixedFloat<PREC>`, the feature flags, and the `alloc` caveat. Companion:
  `examples/fixed_precision.rs`.
