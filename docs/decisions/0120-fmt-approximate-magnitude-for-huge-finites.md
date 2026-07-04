# ADR-0120: render a past-cap finite as an approximate magnitude `±1e{D}`, not `inf`/`0`

- **Status**: accepted
- **Date**: 2026-07-03
- **Supersedes**: the `inf`/`0` saturation token of ADR-0051 (the cap itself stands)

## Context

`to_decimal_string` and `Display` cap exact decimal rendering at
`MAX_FORMAT_DECIMAL_EXPONENT = 10^6` (ADR-0051): rendering `m·2^E`
exactly needs a `5^|shift|` bridge linear in `E`, and `log10`/`exp10`
(which would give a bounded scientific shortcut) are `exp-log`-gated out
of this `big`-only module. Past the cap the formatter emitted `inf` (too
large) or `0` (too small).

That token is misleading: a **finite** value is not infinite. A finite
`2^(i64::MAX-10)` printed `inf`, and — flagged by the review (pf-nyfz,
pf-f7vg) — `pfloat-ball`'s `to_decimal_interval` could print a finite
interval endpoint as `inf`, making a bounded interval read as unbounded.

Parnell chose (2026-07-03) to render an **approximate finite magnitude**
rather than keep the `inf`/`0` token or defer.

## Decision

Past the cap, render `±1e{D}` where `D` is the decimal exponent estimate
the cap check already computes (`approximate_log10_floor`, the rational
`log10(2) ≈ 30103/100000`). The value's sign is the leading `-`; `D`
carries its own sign, so a tiny value renders `±1e-{D}`.

The mantissa is **not resolved** — a documented saturation. Without
`exp10` the leading digits `10^frac` cannot be formed at bounded cost, so
the token uses a leading `1`. It is a valid, finite decimal literal that
conveys the order of magnitude; its low-order exponent digits are
approximate at astronomical magnitudes (the `30103/100000` residual), and
it is **mode-independent** (not a directed interval bound). This is the
explicitly accepted trade for a finite, honest-about-being-large rendering
over the misleading `inf`.

Both saturating sites (`format_normal` for `to_decimal_string`/`Display`,
and `format_shortest` for `to_shortest_decimal_string`) route through the
one `approximate_magnitude_string` helper. The `Saturation` enum and the
`inf`/`0` `saturated_string` are removed; genuine `Class::Infinity`/`Zero`
still render `inf`/`0` via `special_string` (unchanged).

## Consequences

- A finite past the cap now reads as a finite magnitude (`1e331…`,
  `-1e-331…`), never `inf`/`0`; ball printed intervals no longer show a
  finite endpoint as unbounded. The two fmt unit tests that encoded the
  old `inf`/`0` token are updated to assert the approximate rendering.
- The parse↔fmt saturation symmetry ADR-0051 noted is now deliberately
  asymmetric: parse still saturates to a value (mode-aware ±inf/±max/±0/
  ±MinPos, ADR-0117), fmt renders an approximate magnitude string. They
  serve different roles (a value vs a human-readable label); the round
  number is no longer round-trip-identical past the cap, which was
  already only order-of-magnitude meaningful.
- Not a directed bound: a follow-up could make the estimate directed
  (ceil for an upper interval end, floor for a lower) for sound printed
  intervals, and resolve the mantissa when `exp-log` is available. Filed
  as an enhancement; out of scope for the honesty fix.

### Inversion (failure paragraphs considered)

- *"Keep `inf`/`0`; it round-trips with parse."* Rejected by the chosen
  option: the round-trip was only order-of-magnitude meaningful past the
  cap (the value is unrepresentable in-budget either way), and printing a
  finite as `inf` is the actual defect, materially so for ball intervals.
- *"Render the exact leading digits."* Impossible at bounded cost under
  `big` alone (`exp10` gated out); the approximate mantissa is the
  accepted saturation.
- *"Make it a directed bound now."* Deferred: the honesty fix (finite,
  not `inf`) is separable from interval-print soundness, and the latter
  wants an `exp-log`-backed mantissa to be worth the sign/direction
  matrix.

## References

- pf-nyfz (P2), pf-f7vg (P3), epic pf-8iji; the user's display-format
  decision 2026-07-03.
- ADR-0051 (the format cap; its `inf`/`0` token superseded here),
  ADR-0117 (parse's mode-aware past-cap saturation).
- `src/fmt.rs` (`format_normal`, `format_shortest`,
  `approximate_magnitude_string`, `approximate_log10_floor`).
- `src/fmt.rs` tests `format_cap_renders_finite_huge_as_approximate_magnitude`,
  `format_cap_boundary_under_renders_over_saturates`.
