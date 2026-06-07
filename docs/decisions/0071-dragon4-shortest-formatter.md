# ADR-0071: Dragon4 shortest-output formatter

- **Status**: accepted
- **Date**: 2026-06-06

## Context

`BigFloat`'s decimal formatter emits a fixed digit count from
`round_trip_digit_count(precision)` (`ceil(p · log10 2) + 1`). That
output is correct and round-trip-safe but not minimal: it emits enough
digits to guarantee a round trip, not the fewest that round-trip.
ADR-0029 deferred the shortest-output formatter to 1.x to keep the v1.0
critical path free of correctness risk on a working component. Roadmap
Phase 3 (adoption polish) is that 1.x slot.

The robustness rewrite (ADR-0051) left the digit machinery factored
into reusable pieces (`compute_scaled`, `int_to_decimal`, `pow5`, the
magnitude cap), so the shortest algorithm is a digit-selection swap on
that seam rather than a rewrite.

## Decision

Add `BigFloat::to_shortest_decimal_string` (and the `FixedFloat<PREC>`
delegate) implementing the Steele-White / Burger-Dybvig free-format
algorithm: the shortest decimal that parses back to the value at its
own precision under round to nearest, ties to even. `Display` is left
unchanged (it keeps the round-trip-safe count), so this is an opt-in
method, not a behavior change; no working output is put at risk.

The algorithm runs on the existing big-integer limb routines
(`multiply_limbs`, `divmod_limbs`, `pow5` and shifts for powers of ten
and two, `cmp_limbs`). It sets up `R/S/M+/M-` from the mantissa integer
and exponent, with the unequal-gap adjustment when the value is a power
of two; estimates the decimal exponent from the existing rational
`log10` and fixes it up; then generates digits, stopping at the
round-to-nearest-even boundaries (inclusive when the mantissa is even,
the standard ties-to-even rule). The magnitude cap (ADR-0051) is
preserved: a finite value past it saturates to `inf` / `0`, since the
scale step would otherwise build an unbounded power of ten.

Tie-breaking on an exact dyadic halfway value (two equal-length outputs
equidistant from the value) rounds half to even on the decimal digit,
matching IEEE round-to-nearest-even. This is a measure-zero case where
other formatters (including Rust's own `f64` `Display`) may pick the
other neighbor; both are valid shortest round-trips.

## Consequences

Callers can now request the minimal-digit form while `Display` stays
stable. The implementation is verified in `tests/fmt_shortest.rs`: a
200000-value random `f64` sweep plus curated cases check, for each
value, that the output parses back to exactly that `f64` and has the
same significant-digit count as Rust's own shortest formatting (which
is minimal); arbitrary-precision cases (`1/3`, `1/7` at p up to 256)
check round-trip through pfloat's parser and that the digit count never
exceeds the round-trip-safe bound. Rust's shortest `f64` output is a
strong independent oracle for the precision-53 path.

f64 subnormals are out of scope for the Rust oracle: they carry fewer
than 53 bits, but pfloat has no subnormals, so a precision-53 value
always carries 53 bits and its shortest form is correct but longer than
f64's subnormal shortest. This is a precision-model difference, not a
formatter disagreement.

The cost is a new digit-selection path (~200 lines) on a previously
frozen-correct formatter; the opt-in method and the test coverage keep
that bounded. The shortest method is round to nearest only; directed-
mode shortest output is exotic and deferred.

## Related

- Plan: `~/.claude/plans/melodic-knitting-ripple.md` (Phase 3, slice C3d)
- Issue: `pf-bwm1`
- Other ADRs: supersedes ADR-0029 (the deferral this resolves); builds
  on ADR-0051 (the formatter seam) and ADR-0001 / ADR-0002 (mantissa
  representation)
- References: Steele & White, "How to Print Floating-Point Numbers
  Accurately" (PLDI 1990); Burger & Dybvig, "Printing Floating-Point
  Numbers Quickly and Accurately" (PLDI 1996)
