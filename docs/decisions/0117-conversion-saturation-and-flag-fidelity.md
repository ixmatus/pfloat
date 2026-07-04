# ADR-0117: conversion & parse fidelity — mode-aware parse-cap saturation, tininess-before-rounding underflow, and the sticky-flag lane

- **Status**: accepted
- **Date**: 2026-07-03

## Context

Three conversion/parse fidelity gaps from the 2026-06 review
(epic pf-8iji):

1. **pf-mw6u — parse cost-cap saturation is mode-blind.**
   `finite_to_bigfloat` (src/parse.rs) caps the decimal exponent at
   `MAX_DECIMAL_EXPONENT = 10^6` (ADR-0031, a `pow5` storage budget) and
   saturated past it to `±∞`/`±0` regardless of mode. But that cap fires
   far inside pfloat's exponent range — `10^(10^6) ≈ 2^3.3e6`, versus the
   `i64::MAX ≈ 9.2e18` binary-exponent rim — so the saturated value is a
   **directional approximation of a representable value pfloat declines to
   build**, not a true overflow. Mode-blind saturation broke directed
   rounding: `1e2000000` under TowardZero returned `+∞` (must stay
   finite), `1e-2000000` under TowardPositive returned `+0` (must round up
   to a nonzero).

2. **pf-k08c — `to_f32/f64` missed underflow on round-up-to-normal.**
   `to_ieee` (src/convert.rs) flagged UNDERFLOW from the DELIVERED
   exponent (`er < emin`), which misses the case where a subnormal input
   rounds UP to the smallest normal: `to_f32(2^-126 − 3·2^-152)` → the
   smallest normal `0x00800000` with INEXACT only, no UNDERFLOW, though
   the value is tiny (`< 2^-126`) before rounding.

3. **pf-sjgh — the conversions never fed the sticky-flag lane.**
   `to_f32_round`/`to_f64_round` returned the per-call Status but never
   called `auto_raise`, so `pfloat::flags::test()` stayed 0 for an
   sNaN→f32 INVALID — contradicting their own docstring ("the IEEE
   754-2019 sticky flags the conversion raised").

## Decision

1. pf-mw6u: make the parse-cap saturation mode-aware, mirroring the
   `overflow`/`tiny` conventions in `convert.rs`. Overflow side: round to
   `∞` only under nearest or toward `±∞` in the value's own direction;
   otherwise the largest finite. Underflow side: round to the smallest
   nonzero only toward the value's own sign; otherwise `±0`. OVERFLOW /
   UNDERFLOW | INEXACT as before. The saturated magnitude remains a
   cost-cap approximation (the true value is representable but unbuilt);
   the fix corrects only the *direction*, keeping the result on the
   correct side of zero/∞ under every mode.

2. pf-k08c: adopt IEEE 754 §7.5 **tininess detected before rounding**
   (the recommended default): flag UNDERFLOW when the value's own
   exponent is below the format's `emin` (`|value| < 2^emin`) and the
   result is inexact. This catches the round-up-to-normal case the
   delivered-exponent test missed; both IEEE detection methods call this
   input tiny, and before-rounding is adopted uniformly.

3. pf-sjgh: `auto_raise(status)` in `to_f32_round`/`to_f64_round` so the
   conversion's INVALID/INEXACT/OVERFLOW/UNDERFLOW reach the thread-local
   sticky-flag lane, as every arithmetic kernel already does and as the
   docstring promises.

## Consequences

- Parse and conversion now honour the rounding mode and the IEEE
  underflow contract, and conversions participate in the sticky-flag
  lane. No value change for in-range parses or non-tiny conversions.
- Verified against MPFR: `differential_add`/`differential_sub` unchanged
  (the addsub siblings); the parse/convert paths are covered by the three
  red-before/green-after reproducers.

### Deferred: pf-nyfz / pf-f7vg (fmt renders a huge finite as `inf`)

`to_decimal_string` and `Display` saturate a finite value past the format
cap to `"inf"`/`"0"` (ADR-0051), which is misleading for a finite (e.g.
`2^(i64::MAX-10)` prints `inf`), notably in `pfloat-ball`'s
`to_decimal_interval`. Making this mode-aware to match the now-mode-aware
parse runs into a hard constraint: there is **no bounded-cost finite
decimal rendering** for such a value under the `big`-only feature set
(`log10`/`exp10` are `exp-log`-gated), and changing the saturation token
is a change to a **public output format** (a one-way door). This warrants
an explicit display-convention decision rather than a unilateral change;
it is left open under pf-nyfz (P2) and pf-f7vg (P3) with this ADR
recording why. The `parse`↔`fmt` saturation symmetry ADR-0051 noted is
now intentionally one-sided pending that decision.

### Inversion (failure paragraphs considered)

- *"pf-mw6u should raise the cap and compute the value."* Rejected: the
  cap is ADR-0031's deliberate storage budget; recomputing an arbitrarily
  large `pow5` is the unbounded cost the budget exists to prevent. The
  mode-aware saturation is the in-budget fix.
- *"Detect tininess after rounding (using `er`)."* Rejected: it is a
  valid IEEE option but is exactly what missed the round-up-to-normal
  case; before-rounding is the recommended default and flags it.
- *"auto_raise inside `to_ieee` instead of the public methods."* Either
  works; the public methods are the single public entry points, so
  raising there covers every path (tiny/overflow/normal/nan) once.

## References

- pf-mw6u, pf-k08c, pf-sjgh (epic pf-8iji); reproducers Q2a/Q2b, L1, Q9
  in the review harness. pf-nyfz, pf-f7vg deferred (above).
- ADR-0031 (parse `pow5` budget), ADR-0051 (fmt saturation), the
  `convert.rs` `overflow`/`tiny` mode conventions mirrored by parse.
- `src/parse.rs` (`finite_to_bigfloat`), `src/convert.rs` (`to_ieee`,
  `to_f32_round`, `to_f64_round`).
- `tests/regression_review_2026_06_10.rs`:
  `parse_cap_saturation_is_mode_aware`,
  `to_f32_tiny_before_rounding_flags_underflow`,
  `to_f32_feeds_the_sticky_flag_lane`.
