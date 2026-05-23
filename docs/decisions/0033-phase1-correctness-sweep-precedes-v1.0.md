# ADR-0033: Phase 1 correctness sweep runs to completion before the v1.0 tag

- **Status**: accepted
- **Date**: 2026-05-22

## Context

The Phase 1 correctness-sweep plan, accepted on 2026-05-22 and
filed at `docs/decisions/plans/phase-1-correctness-sweep.md`, was
originally scheduled as v1.0 → v1.1 work. The plan's "Sequencing
relative to the v1.0 tag" section spelled out the reconciliation:

> v1.0 ships with the count-only `## Conformance evidence` section
> (the Phase 8b artifact). The disclosure text says explicitly:
> "exhaustive `f32` correctness audit in flight; the per function
> rounding status table is the next minor."

> Phase 1 sweep is v1.0 → v1.1 work, not pre v1.0. v1.0 stands
> alone as "implementation complete, conformance evidence
> describes the kinds of verification we run." v1.1 is
> "exhaustive `f32` audit complete, per function status table
> replaces the count-only section."

> `has-errors` findings in the sweep become v1.1 fixes, not v1.0
> blockers; v1.0 ships under the stake "we run these kinds of
> verification, the exhaustive audit is the next thing."

The motivating concern was timeline: the plan's own CPU-budget
estimate put a complete sweep at six to twelve months, with no
published 0.x consumers waiting through. The trade-off was
"published 1.0 sooner, sharper conformance evidence in 1.1" against
"ship under a less rigorous v1.0 claim and accept that the audit
might surface defects later."

Slice 8b (the conformance-evidence preparation) surfaced two
findings during execution that re-weight that trade-off.

**The first finding: pfloat's `exp` kernel mis-rounds inputs at
the binary64 underflow boundary.** The Lefèvre-Muller corpus
(now sourced from CORE-MATH under MIT per the slice-8b
provenance document) opens its `exp.wc` with an "exercise
underflow or overflow" block. On input `0xc0874385446d71c3`
(approximately −744.44), CORE-MATH and mpmath agree the
correctly-rounded `NearestEven` result is `0x0000000000000001`
(`f64::MIN_POSITIVE_SUBNORMAL`); pfloat returns a different
result. The corpus as shipped at slice 8b excludes that block
to keep the suite green; the underflow path is a documented
known limitation in `src/math/exp.rs` lines 23-30. This is a
real wrong-rounding defect on the public surface of a library
whose README headline claims correctly-rounded transcendentals.
The Phase 1 exhaustive `f32` sweep over `exp` would surface
this kind of defect on every input it touches in the underflow
region — it is precisely the class of finding the audit is
designed to catch.

**The second finding: the tier-2 specials' differential verification
is structurally weaker than the README posture implies.** The
slice 8b L-M corpus added 450 hard-to-round cases across nine
elementary functions (`exp`, `ln`, `sin`, `cos`, `tan`, `atan`,
`asin`, `acos`, `exp2`), all of which pfloat passes. CORE-MATH
covers 41 binary64 functions, including the tier-2 specials for
which MPFR has no primitive at all (modified Bessel `I`/`K`,
Airy, `Si`/`Ci`/`li`). Those functions' existing differential
lanes use an mpmath table at `p ≤ 256` plus identity cross-ties
(Wronskian, I/K reciprocity) plus dyadic self-consistency. That
posture catches gross errors but does not prove correct rounding
the way the bit-exact-vs-MPFR tier does for the standard
elementary surface. The slice 8b artifact made the gap legible;
extending the L-M corpus to the full pfloat-supported surface
would surface more.

The combination changes the cost-benefit. Shipping a v1.0 that
claims "correctly-rounded arbitrary-precision floats" while
carrying a known wrong-rounding case in `exp` and a structurally
weaker rigor posture on half the special surface than the README
implies is a credibility error. Once 1.0 is published on
crates.io it is immutable (yank-only, no replace), so the first
externally-visible version of the library is the version users
will quote, link to in their own changelogs, and judge the
permacomputing-horizon discipline against. The credibility cost
of "we shipped a wrong-rounding case in our headline kernel" is
permanent in a way that the timeline-cost of the audit is not.

## Decision

Phase 1 runs to completion before the v1.0 tag. The Phase 1
plan's "Sequencing relative to the v1.0 tag" section is updated
in this same commit set to record the new ordering. Slice 8c
(the v1.0 tag + `cargo publish` slice) is parked indefinitely;
its beads remain on file but their dependencies extend to the
Phase 1 exit criteria.

The v1.0 ship criterion changes:

- Every frozen unary function (the surface enumerated in the
  Phase 1 plan) has a status table row with a definitive
  `rounding_status`.
- No function in the v1.0 surface is `has-errors`. Each is
  either `correctly-rounded` (the headline claim, earned per
  function) or documented `faithful` with rationale (the
  honestly-downgraded claim, surfaced rather than buried).
- The L-M corpus is extended to cover every pfloat-supported
  binary64 function CORE-MATH has a `.wc` file for, not just
  the nine slice 8b sampled.
- The Arb backend (via `python-flint` as a long-lived
  subprocess) is integrated alongside the MPFR backend behind
  the `Enclosure` trait the plan already specifies, with
  per-function routing per the plan's table (MPFR-primary for
  functions MPFR has a primitive for; Arb-primary plus identity
  cross-check for the tier-2 specials).
- The per-function status table is published in README and docs;
  the count-only `## Conformance evidence` section ships only as
  long as Phase 1 is mid-execution and is replaced by the table
  at v1.0.

The Phase 1 exit criteria from the plan (lines 360-375 of
`docs/decisions/plans/phase-1-correctness-sweep.md`) become the
v1.0 ship criteria with the substitution above.

## Consequences

**Honest framing.** The published v1.0 is the version that
substantiates the README's correct-rounding claim with an
exhaustive `f32` audit and a per-function table, not the version
that promises the audit is forthcoming. A future stranger reads
the README and the published status table and gets a per-function
answer to "can I depend on this for `gamma`/`bessel_i`/`Airy`?"
without reading the source.

**Timeline cost, accepted.** The Phase 1 plan's CPU-budget
estimate (hours for cheap functions, days for trig and the
expensive specials, weeks total once Arb-subprocess overhead is
folded in) means v1.0 ships months later than the slice 8c plan
contemplated. pfloat has no published 0.x consumers; this cost
is borne entirely by the project, not by users routed through a
v1.x migration.

**Slice 8c parked.** The five unsigned commits and the YubiKey
touch sequence drafted in `~/.claude/plans/abundant-yawning-badger.md`
do not run until Phase 1 completes. The 8c beads (`pf-4fi`
through `pf-61n`) stay on file; the disclosure-correction diff
artifact at `docs/disclosure-correction-v1.0.diff` stays in
tree (the two factual corrections it makes are still required
at the eventual tag).

**The exp underflow defect is a v1.0 blocker.** Per the new
ship criterion `has-errors` rows are blockers. The Phase 1
sweep is the right surface to catch this defect, but it is
known today; the kernel-side fix is bounded (audit the
`src/math/exp.rs` underflow path against the binary64 subnormal
boundary semantics) and can land before the sweep does, with
the sweep's job becoming regression-verification.

**The CORE-MATH corpus expands.** The L-M differential tier
shipped at slice 8b covers nine functions. Coverage extends to
every pfloat-supported binary64 function CORE-MATH carries a
`.wc` file for: `exp10`, `expm1`, `log1p`, `log2`, `log10`,
`sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`, `erf`,
`erfc`, `tgamma`, `lgamma`. Pre-sweep work; bounded; ~5-10×
the slice 8b corpus.

**Arb backend integration.** The Phase 1 plan specifies the
`python-flint` subprocess worker with a streaming pipe protocol
and a binding caveat about exposing the ball's radius. That
work moves from "v1.0 → v1.1" to v1.0-blocking; it is the
load-bearing piece for upgrading the tier-2 specials' rigor
from "matches mpmath at p ≤ 256 plus identity cross-ties" to
"the certified ball bracket determined the rounded `f32` under
the rounding mode." The licensing posture (FLINT and Arb are
LGPL; oracle is a CI lane only, never linked into the shipped
crate) is unchanged from the plan.

**ROADMAP.md updates.** The "Currently in flight" section is
revised: Phase 1 is in flight, not queued; in-tree Phase 8
slice 8c is parked behind Phase 1, not ahead of it. The Phase 3
language about being sequenced "after the in repo Phase 8 v1.0
tag and after the Phase 1 status table" simplifies to "after
the v1.0 tag" — the status table now lands at v1.0.

## Alternatives considered

**Ship v1.0 with the count-only conformance section, audit
during 1.x.** The original plan. Rejected on the
credibility-cost argument above: a published 1.0 is immutable
on crates.io and is what users judge the project against
permanently.

**Fix only the exp underflow defect; ship slice 8c without
the broader sweep.** Closes the one known blocker without
absorbing the full audit cost, but does not address the
structurally weaker tier-2 specials verification. v1.0 would
ship with the same "matches mpmath at p ≤ 256" rigor for the
half of the special surface MPFR has no primitive for, while
the README claims correctly-rounded transcendentals across the
board. The credibility gap stays.

**Run the sweep partially.** Sweep the elementary functions
that have MPFR primitives; ship v1.0 with the tier-2 specials
on their existing posture. Same gap as the second alternative
in different words. Rejected.

## References

- `docs/decisions/plans/phase-1-correctness-sweep.md` — the work
  breakdown (whose "Sequencing relative to the v1.0 tag" section
  is updated in this same commit set).
- `docs/ROADMAP.md` — Track B sequence (also updated in this
  commit set to reflect the re-sequencing).
- `docs/lefevre-muller-corpus-provenance.md` — the CORE-MATH
  source identification and MIT-attribution posture that the
  expanded corpus inherits.
- `tests/differential_lefevre_muller.rs` — the slice 8b
  differential tier that surfaced the `exp` underflow defect
  and motivated the re-sequencing.
- `src/math/exp.rs` lines 23-30 — the durable note documenting
  the underflow defect, written at slice 8b.8 and now upgraded
  from "v1.x finding" to "Phase 1 blocker."
- ADR-0029 — Dragon4 shortest formatter deferred to v1.x; the
  deferral stands and is now post-v1.0-Phase-1 (Phase 3 of the
  roadmap).
- ADR-0032 — libm reciprocal and root kernels stay direct
  primary not aliased; the only other Phase 1/Phase 2 decision
  discharged early.

## Commits

- The ADR itself, the plan-document update, the roadmap update,
  and the memory pivot land as one signed merge of unsigned
  branch commits on `slice-phase1-presequence`. The implementation
  work (corpus expansion, Arb backend, exhaustive harness)
  follows in subsequent slices once the next session opens.