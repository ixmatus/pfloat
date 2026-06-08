# ADR-0082: `pfloat-ball` interval-input Arb bracket: range soundness and tightness

- **Status**: accepted
- **Date**: 2026-06-07

## Context

The independent Arb containment lane (ADR-0078, pf-fe5f) brackets each ball
operation at sampled witness POINTS inside the input ball. Point sampling is
structurally blind to two things ADR-0078 deferred to pf-fe5f.7. First, range
soundness for the non-monotonic functions: five witnesses cannot see a result
ball that fails to enclose an interior extremum of `sin`/`cos`/`tan` (no
witness reaches `sin`'s peak at `pi/2`, but the true image does). Second,
tightness: the witness-image span is rounding-dominated for near-exact inputs
and so measures nothing about how loose the enclosure is.

Both need Arb's rigorous image of `f` over the whole input INTERVAL rather than
at sampled points. The worker's `BRACKET` verb took a zero-radius point only.

## Decision

Extend the worker and the lane with an interval-input bracket, and assert
soundness only where the assertion is provably clean.

**The `BRACKETI` verb (`scripts/arb_oracle_worker.py`).** A new verb takes an
input radius per operand and brackets `f` over `[mid - rad, mid + rad]`. It
builds a rigorous Arb interval ball as the union of the two exact dyadic
endpoints (rigorous even if an endpoint is inexact at the oracle precision,
since the union of two balls contains both and the span between them), then
runs the existing `dispatch_elementary` and shares one reply encoder with
`BRACKET`. A zero radius collapses to the exact midpoint, so a degenerate
`BRACKETI` reduces to the point `BRACKET` bit-for-bit (pinned by a test). The
`cbrt` dispatch gains a straddle case: the real odd cube root is monotone
through `0`, so for a `0`-containing ball its image is
`[cbrt(lower), cbrt(upper)]` (Arb's principal `root(3)` is NaN on the negative
half).

**The interval predicate is SUPERSET, the opposite of the point lane.** For
input interval `B` the worker returns Arb's rigorous image enclosure `J` with
`J superset-or-equal f(B)`. Range soundness requires `result superset-or-equal
f(B)`. The point lane uses overlap, which is correct for a point oracle but
UNSOUND for interval input: a ball that misses an interior extremum still
overlaps `J`. The sound check is `result superset-or-equal J` (passing it
proves `result superset-or-equal J superset-or-equal f(B)`). The direction was
confirmed by verify-the-verdict, not assumed.

**The hard superset assertion is scoped to extremum straddles.** Arb's interval
image carries an outward overshoot, because the input ball's radius is an
inflated roughly 30-bit `mag` and a steep `f` propagates that overshoot. A
correct result ball can therefore be TIGHTER than `J` away from the extrema:
measured, at `p = 113` the great majority of general balls have `result` not a
superset of `J`, and `exp` fails it on every sample. So the lane asserts the
superset only at an extremum straddle, where `|f'| -> 0` makes `J` tight while
pfloat-ball's Lipschitz radius stays wide (measured: zero false-fails across
240 straddles at `p` in `{24, 53, 113}`). A negative control builds the naive
endpoint-only enclosure (treating `sin`/`cos` as monotonic) and confirms the
superset check rejects the missed extremum. This scope was set by reproducing
the result-versus-image relationship empirically before asserting, rather than
shipping a general superset lane that false-fails correct code.

**Tightness is MEASURED per bucket, not asserted.** The slack
`log2(width(ball) / width(J))` is recorded as an expected COUNT PER BUCKET
(`(function, precision, magnitude)`) in a checked-in baseline, compared cell by
cell so a regression in one bucket cannot hide behind another's improvement.
The baseline reads as a fingerprint of the three enclosure shapes: monotonic
functions cluster at slack `0` (their endpoint enclosure matches Arb's image),
the Lipschitz functions spread over the `1/|f'|` looseness of their radius, and
the composed `tan` sits between. The lane is `#[ignore]` because the image
width is the worker's and so the exact counts are Arb-version-sensitive; the
per-release sweep runs it with `--ignored` under the pinned venv, and a
regeneration mode rewrites the baseline for an intended enclosure change.

## Consequences

- The missed-interior-extremum blind spot is now directly checked for `sin` and
  `cos`, with a negative control proving the check has teeth. This is the first
  lane that is sound against a class of bug the exhaustive scalar sweep guards
  only at the scalar level.
- The honest scoping is load-bearing and recorded in the predicate's doc
  comment: the superset is asserted only at extrema, and the general width
  relationship is a measurement, precisely because a correct ball can be tighter
  than Arb's inflated interval image. A reader is not told the lane proves more
  than it does.
- The tightness baseline is a per-release regression artifact, not a per-push
  gate, and it is tied to the pinned Arb build. The cost is that an Arb upgrade
  requires a deliberate regeneration; the benefit is that an enclosure-tightness
  regression in any single bucket is caught.
- `pfloat-ball` stays `0.1.0`: the surface is test-only (the `differential-arb`
  feature), and pfloat-core is untouched.

## Related

- Plan: `plans/nested-prancing-lovelace.md` (S1, S3, S4).
- Commits: `5a1a063` (S1 verb and driver), `2cd7342` (S3 range soundness),
  `52804ed` (S4 tightness table).
- Beads: `pf-fe5f.7` (under epic `pf-fe5f`).
- Other ADRs: extends ADR-0078 (the verification design and the deferral);
  the radius soundness it leans on is ADR-0077; the seeding that feeds these
  lanes the hardest inputs is ADR-0083.
