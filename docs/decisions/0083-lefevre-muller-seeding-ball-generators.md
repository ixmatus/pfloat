# ADR-0083: Lefèvre–Muller hardest-to-round seeding of the ball generators

- **Status**: accepted
- **Date**: 2026-06-07

## Context

The `pfloat-ball` verification lanes (the FTIA self-consistency lane and the
independent Arb containment lane) drew their inputs from `random_ball`: integer
midpoints scaled by a power of two. Those midpoints are arithmetically bland.
A hard-to-round scalar, where `f(mid)` sits within sub-ULP of a rounding
boundary, stresses the directed endpoints that define the ball radius far
harder than a random integer does, because the radius round-up has to be
correct exactly where the kernel's own residual is largest. ADR-0078 named the
Lefèvre–Muller seeding as a deferred enrichment of the generators.

pfloat-core already carries the corpus: `tests/differential/lefevre_muller_data.rs`
holds CORE-MATH hard-to-round binary64 inputs (MIT, transcribed under
attribution) with outputs recomputed independently by mpmath. The ball crate
did not yet draw on it.

## Decision

Seed the ball generators from the committed CORE-MATH corpus, reused verbatim,
and wire it into both lanes additively.

**Reuse, do not duplicate.** `pfloat-ball/tests/common/mod.rs` pulls the corpus
in through `include!` of pfloat-core's committed data file, so there is one
source of truth and the MIT attribution header travels with it. Only the input
field of each case is used; the expected-output field is ignored (the ball
lanes verify enclosure, not the scalar value).

**A bit-exact midpoint builder.** `bf_of_f64_bits` constructs the `BigFloat`
midpoint as the integer significand times an exact power-of-two scale, at
precision `p >= 53` where the 53-bit value stays exact. This is the bit-exact
route; a shortest-decimal round-trip would round at the finer `BigFloat`
precision and soften the hard-to-round property the seed exists to exercise. A
test round-trips a spread of seeds back to `f64` to pin the builder.

**Per-function corpus, additive wiring.** `lm_cases_for` maps a ball function
id to its own corpus, because an input hard to round for `f` makes `f(mid)`
boundary-close, which is what stresses `f`'s ball endpoints. `seeded_ball`
places a corpus midpoint at precision 53 or 113 with a small radius bound to
the magnitude. New tests in both lanes (`ftia_unary_self_consistency_hard_to_round`
and `arb_containment_unary_hard_to_round`) run the existing soundness checks on
seeded balls without changing any existing test, so the slice is purely
additive.

## Consequences

- Both ball lanes now exercise the hardest available inputs: the seeded Arb
  lane checked roughly 3600 hard-to-round witness brackets across 16
  corpus-backed unary functions with zero unsoundness.
- One corpus, one builder, no duplication: a future correction to the corpus or
  the conversion improves every lane at once, and the license attribution is
  not copied.
- The seeding is limited to `p` in `{53, 113}`, where the 53-bit binary64 value
  is exact; `cbrt` and `sqrt`, which have no corpus, and the inverse-trig edge
  functions keep random inputs. This is a deliberate scope, not a gap: the
  corpus is an elementary-function hard-to-round table.

## Related

- Plan: `plans/nested-prancing-lovelace.md` (S2).
- Commits: `9f202e6`.
- Beads: `pf-vcqh` (discovered from `pf-fe5f`).
- Other ADRs: feeds the lanes ADR-0082 strengthens; the corpus provenance
  discipline mirrors pfloat-core's `extract_lefevre_muller.py`.
