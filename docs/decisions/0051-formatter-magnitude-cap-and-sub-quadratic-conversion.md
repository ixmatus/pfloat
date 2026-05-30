# ADR-0051: Decimal formatter magnitude cap and sub-quadratic conversion

- **Status**: accepted
- **Date**: 2026-05-30

## Context

The 2026-05-29 correctness review found two reachable defects in the
decimal formatter (`src/fmt.rs`), both triggered by ordinary arithmetic
rather than by parse.

A `BigFloat` is `m · 2^E` with an `i64` exponent and no `emax`, so
repeated squaring drives `E` toward `i64::MAX` (for example `2^(2^40)`
after forty squarings has `E ≈ 1.1e12`). `format_normal` rendered the
full exact decimal by building a numerator and denominator scaled by
`5^p5 · 2^p2`. For a huge `E` the power-of-two width overflowed the `u32`
used for the shift (a debug panic at `fmt.rs:256`) and the `5^|shift|`
factor allocated without bound (a release out-of-memory). There is no
bounded exact rendering of such a value: the leading decimal digits
intrinsically need the full `5^|shift|` factor, whose size is linear in
`E` (`10^d = 2^d · 5^d`, and the five-part does not cancel). The bounded
scientific shortcut (`E · log10(2)` for the exponent, `exp10(frac)` for
the digits) is unavailable, because `log10` and `exp10` are behind the
`exp-log` feature while the formatter builds under `big` alone.

Separately, `int_to_decimal` extracted digits one at a time by repeated
division by ten, quadratic in the digit count. `Display` emits
`round_trip_digit_count(precision)` digits, which grows with the
caller-chosen precision, so `Display` at precision near `10^6` ran for
tens of seconds: a caller-reachable denial of service.

The threat model is that `Display` and `to_decimal_string` take an
arbitrary finite `BigFloat` and must never panic, exhaust memory, or
hang. A bounded wrong-or-error result past a documented cap is
permitted; an unbounded one is not.

ADR-0029 deferred the shortest-output (Dragon4) formatter as an output
aesthetic. The two defects here are robustness, not aesthetics, so the
review pulled them ahead of v1.0 as slice `pf-vbm2` while the
shortest-output half stays deferred under ADR-0029.

## Decision

1. Cap the formatter at `MAX_FORMAT_DECIMAL_EXPONENT = 1_000_000`,
   value-matched to parse's `MAX_DECIMAL_EXPONENT`. Within the cap the
   renderer is unchanged (exact, round-trip-safe). Past it a finite value
   saturates the way parse does: a magnitude above the cap renders as
   `inf` (or `-inf`), below it as `0` (or `-0`). The shared value makes
   the exactly-renderable range identical to parse's round-trippable
   range, so `parse(format(x)) == x` holds across the whole in-cap region
   with no new boundary, and a value past the cap reads back as the same
   saturated token parse would itself produce for a decimal that large or
   small.

2. Replace `int_to_decimal` with a divide-and-conquer base conversion
   (Brent & Zimmermann, *Modern Computer Arithmetic* §1.7). Over the
   sub-quadratic `divmod_limbs` (ADR-0052) it is `O(M(D)·log D)` in the
   digit count `D`. No digit or precision cap is needed: the cost now
   tracks the requested output size.

3. The cap is a fixed default, not a user-supplied knob. The renderable
   range coincides with parse's, so raising only the format cap would
   produce strings parse cannot round-trip; raising both is a separable,
   backward-compatible future addition tracked as `pf-hn9s`. The cap
   check, the rational scaling in `compute_scaled`, and the
   `int_to_decimal` primitive are factored so the deferred Dragon4 slice
   (ADR-0029) swaps only the digit-selection strategy.

The beyond-cap behaviour is deliberately lossy: a finite value renders as
`inf` or `0`. That is a resource bound, not a claim the value is infinite
or zero. It is documented on `to_decimal_string` and is the same place
parse puts any decimal past its own cap. A value past the cap has a
magnitude beyond every IEEE interchange format and beyond what parse can
represent, so no decimal string round-trips to it regardless.

## Consequences

- `Display` and `to_decimal_string` are now total over every finite
  `BigFloat`: bounded cost, no panic, no out-of-memory, no hang.
- Every value within the parse-round-trippable range still renders
  exactly and round-trips, unchanged from before.
- A finite value with `|decimal exponent| > 10^6` renders as a saturation
  token rather than its digits. This is lossy by design and documented.
  Callers needing the order of magnitude of such values, or a non-lossy
  "too large to print" signal, are served by the future override
  (`pf-hn9s`) or by the deferred shortest-output work, not by v1.0.
- The rational `log10(2)` magnitude estimate can disagree with the true
  decimal exponent by less than one only at the irrational product's
  integer boundaries, so a value at a true decimal exponent of exactly
  `±10^6` may saturate. The output stays within the threat model.
- The new conversion primitive, the cap guard, and `compute_scaled` are
  the seam for the deferred Dragon4 formatter; that slice does not
  relitigate the cap or the base conversion.

## Rejected alternatives

- **Order-of-magnitude scientific past the cap** (`1e1000001` with a
  placeholder mantissa). The decimal exponent from the rational estimate
  drifts past about `10^7`, and the leading digits are not obtainable at
  bounded cost without the `5^|shift|` factor, so the mantissa would be
  fabricated. Rejected: it invents partially-correct output.
- **An error or sentinel token.** A non-numeric `Display` surprises
  callers and breaks the `parse(format(x))` invariant the fmt fuzz target
  asserts. Rejected for no benefit over saturation.
- **A self-contained fixed-point `log10` / `exp10` inside `fmt`.** This
  would give correct leading digits at bounded cost for any magnitude,
  but it reimplements a transcendental kernel to dodge the `exp-log`
  feature gate and takes exactly the formatter-correctness risk ADR-0029
  defers, for a cosmetic gain. Rejected for v1.0.
- **A digit or precision cap for the quadratic conversion.** Rejected: it
  would refuse to render legitimately high-precision values, a
  completeness gap, where the sub-quadratic conversion removes the cost
  problem without a cap.

## Related

- Plan: `plans/buzzing-stargazing-elephant.md` (pf-vbm2).
- Commits: `423fced` (failing guard), `c028cba` (cap + saturation),
  `f2f79b0` (recursive division), `c98e77c` (divide-and-conquer
  conversion).
- Code: `src/fmt.rs` (`MAX_FORMAT_DECIMAL_EXPONENT`, `saturated_string`,
  `format_normal`, `compute_scaled`, `int_to_decimal`),
  `tests/regression_review_2026_05_29.rs`
  (`panic_fmt_large_exponent_does_not_overflow`).
- Issues: `pf-vbm2` (this slice), `pf-hn9s` (deferred cap override).
- Other ADRs: ADR-0029 (shortest-output deferral, the sibling), ADR-0031
  (parse's `MAX_DECIMAL_EXPONENT`), ADR-0052 (the recursive division this
  conversion rests on).
