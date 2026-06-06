# ADR-0029: Dragon4 / Steele-White shortest formatter deferred to 1.x

- **Status**: superseded by ADR-0071
- **Date**: 2026-05-18

## Context

`BigFloat`'s decimal formatter (`src/fmt.rs`) emits a fixed digit
count derived from `round_trip_digit_count(precision)`, which is
`ceil(p * log10(2)) + 1`. That output is correct and round-trip-safe:
`parse_str` at the same precision recovers the exact value. It is not
shortest: it emits enough digits to guarantee round-trip, not the
minimal number that round-trips. A shortest-output formatter
(Dragon4, or its Steele-White / Grisu / Ryu lineage) would emit fewer
digits while preserving round-trip.

The Phase 7 plan flagged slice 7e as deferral-gated: ship a shortest
formatter only if the v1.0 timeline has room, else document it
deferred to 1.x. This is a Parnell judgment call, taken at plan time.

The relevant scope: replacing the digit-count logic in
`extract_digits()` with a shortest algorithm is roughly 100 to 300
lines, on a formatter that is currently correct, and it introduces a
new round-trip-minimality property surface. The risk is a correctness
regression on a working component, taken on the v1.0 critical path,
for an output-aesthetics improvement that changes no numeric value.

## Decision

Defer the shortest-output formatter to 1.x.

The v1.0 formatter stays as it is: round-trip-correct, not
minimal-digit. No correctness risk is taken on a working component
for v1.0. The decision parallels the 7f.1 deferral (ADR-0028): a
real, scoped improvement is scheduled past the tag rather than
expanding the v1.0 surface and risking the verified behavior.

DESIGN.md previously asserted that the formatter "formats
shortest-round-trip-correct by default". That was inaccurate
(aspirational, not implemented). It is corrected to state the actual
behavior (round-trip-correct, not minimal-digit) with a pointer here,
under the honest-deviation posture: documentation describes what the
code does, not what it is intended to do. The `src/lib.rs` module
doc already states the accurate position ("the current output is
correct, just not always minimal"); it is left as is.

## Consequences

- v1.0 ships a correct, round-trip-safe formatter that emits a few
  more digits than strictly necessary for some values. No numeric
  output changes; only the digit count of the decimal rendering.
- 1.x slice 7e implements Dragon4 / Steele-White shortest output,
  preserving the `to_decimal_string(digits, mode)` and `Display`
  API and extending the `parse(format(x)) == x` property to assert
  minimality. Tracked in the issue graph (`discovered-from` the 7e
  decision) so it is not relitigated from zero.
- A documentation inaccuracy in DESIGN.md is removed rather than
  carried into v1.0; Phase 8 slice 8b's conformance evidence cites
  this ADR as a stated v1.0 deferral alongside ADR-0027 (7d) and
  ADR-0028 (7f).

## Related

- Plan: `~/.claude/plans/quizzical-prancing-lighthouse.md` (Phase 7).
- Code: `src/fmt.rs` (`round_trip_digit_count`, `extract_digits`),
  `src/lib.rs` module doc (already accurate), DESIGN.md "String I/O"
  (corrected with this slice).
- Other ADRs: ADR-0028 (7f; the parallel data-or-risk-backed v1.0
  deferral), ADR-0010 (the defer-invasive-work-past-v1.0 posture).
- References: Steele, G. L., and White, J. L. "How to Print
  Floating-Point Numbers Accurately." PLDI 1990. Burger, R. G., and
  Dybvig, R. K. "Printing Floating-Point Numbers Quickly and
  Accurately." PLDI 1996 (Dragon4).
