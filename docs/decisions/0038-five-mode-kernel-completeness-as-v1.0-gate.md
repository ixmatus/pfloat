# ADR-0038: Five-mode kernel completeness as the v1.0 strong-claim gate

- **Status**: accepted; scaffolding shipped at slice p1.23 (the
  `mpfr_oracle_for_mode` and `certified_round_bf_to_f32` helpers,
  the per-mode status TOML schema, the already-five-mode cohort's
  differential lanes widened to `BIT_EXACT_ROUNDING_MODES`); the
  per-family kernel migrations (p1.24 through p1.34) execute
  against the audit doc.
- **Date**: 2026-05-25

## Context

The pfloat marketing claim is uniformly strong. The Cargo manifest,
the README headline, the keyword list, the `src/lib.rs` crate doc,
and the v1.0 surface document all describe pfloat as "correctly
rounded arbitrary precision floats". The IEEE 754-2019 standard
defines five rounding modes; the unqualified "correctly rounded"
phrase claims agreement with each of them on every input the surface
covers.

The implementation, today, is two tier against that claim.

The arithmetic core (`add`, `sub`, `mul`, `div`, `sqrt`, `fma`,
parse) is correctly rounded under every IEEE 754-2019 mode by
construction. The Ziv driven cohort lifted across slices p1.2 and
p1.4 is also five mode correct: `pow`, `exp`, `ln`, `tanh`,
`lgamma`, `erf`, `bessel_j`. The log2 and log10 entries inherit
through composition with `ln_round`. Everything else on the v1.0
surface is correctly rounded only under NearestEven and faithful
within 1 ULP at tie cases under the four directed modes. The
DESIGN.md "Caveats and open questions" §1 paragraph names the gap
explicitly:

> Kernels still on the fixed 64-bit guard (`exp`, `ln`, `sin`, …)
> carry the original caveat until a later slice extends the driver.

The §1 inventory predates the slice p1.4 `erf` and `bessel_j`
migrations; the structural shape stands. Twenty plus kernels on the
v1.0 surface today carry the directed mode caveat and the uniform
marketing claim does not.

CLAUDE.md frames the resolution as a frugality property:

> Completeness against a spec is a frugality property: half
> implementations consume the canonical name and force every future
> user into workarounds, vendoring, or abandonment, costing more
> total effort than completing the work would have.

The crates.io publish is irreversible (yank only, no replace). v1.0
is the version users link to in their own changelogs, judge the
project's permacomputing horizon discipline against, and depend on
for two decades. Shipping v1.0 with a uniform marketing claim that
is true on the arithmetic core plus seven kernels and faithful at
ties on the remaining ~30 surface entries is a credibility error of
the kind ADR-0033 already paid the timeline cost to avoid for
`has-errors` rows. Shipping the same v1.0 with a narrowed claim
("correctly rounded under NearestEven; faithful at directed mode
ties on the special function cohort") is the cohort labelling
failure mode CLAUDE.md names: the qualifier consumes a section of
the README, every future doc reader has to learn it, and every
downstream user who needed five mode rigor on `gamma` or `bessel_y`
is forced to vendor or wait for v2. Shipping under a delayed v1.0
that closes the gap by extending the implementation matches the
claim without consuming the canonical name on a half story.

The Ziv interval test driver landed at slice p1.2 (ADR-0022) and is
the structural piece this ADR depends on. Its signature today
(`src/math/ziv.rs:96-130`) takes a single argument eval closure, the
target precision, and any of the five `RoundingMode` values; the
interval test returns a correctly rounded value when both ends of
the bounded uncertainty interval round to the same target precision
value under the caller's mode. The infrastructure to migrate every
remaining v1.0 kernel exists; what is missing is the audit work to
identify the per kernel cancellation regimes (the slice p1.4 tanh
short circuit precedent), the per family differential lane widening
to bit exact across all five modes, and the oracle harness work to
verify the result against MPFR, Arb, and mpmath at non NearestEven
modes (the existing bf to f32 Display+parse output bridge is
silently NE only at any verification precision above p=24 per
`feedback_bf_to_f32_directed_mode`).

ADR-0033 sequenced Phase 1 before the v1.0 tag on the same
credibility argument; the Phase 1 sweep closed `has-errors` rows.
This ADR extends the same sequencing posture to directed mode
rounding: Phase 1f closes the directed mode rigor gap before v1.0
in the same way.

## Decision

Phase 1f runs to completion before the v1.0 tag. The v1.0 surface
ships when every kernel on the frozen 47 entry list (per
`docs/v1.0-surface.md`) is correctly rounded under every IEEE
754-2019 rounding mode (NearestEven, NearestAway, TowardZero,
TowardPositive, TowardNegative), verified per family by the
existing oracle and differential harnesses extended to iterate
across all five modes.

The strategic commitment is non negotiable:

1. **No "faithful" `rounding_status` verdict on any row.** Every
   row in `tests/oracle/status/` exits Phase 1f reading
   `correctly-rounded` under every mode. The transitional
   `unswept` state introduced by the per mode TOML schema
   migration (slice p1.23) is internal to the phase and never the
   resting state; no row exits Phase 1f reading `unswept`.

2. **No qualifier added to public docs.** The README's uniform
   claim, the Cargo.toml description, and the lib.rs crate doc do
   not gain a "correctly rounded under NearestEven; faithful at
   directed mode ties on the special function cohort" qualifier at
   any point in the phase. Either the doc reads silent on rounding
   mode (kernel not yet migrated, the README's existing prose
   remains pre 1.0 marketing) or it states the five mode claim
   (kernel migrated). "Pending" is tracker state, not doc state.

3. **No cohort labels in public docs.** Phrases like "five mode
   cohort", "Ziv driven family", and "fixed guard kernels" are
   internal planning vocabulary; they appear in this ADR, in the
   audit document, in bead descriptions, and in per family ADRs,
   but not in README, lib.rs, Cargo.toml, or kernel doc comments.

4. **The intractability fork is drop or extend, not qualify.** If
   a kernel resists the Ziv migration after honest engineering
   work (a regime whose cancellation cannot be soundly bounded
   within `ZIV_ERROR_GUARD=24` slack; an asymptotic branch whose
   reformulation cost exceeds reasonable scope; an oracle gap
   wider than expected), the choices are: **(A)** drop the kernel
   from the v1.0 surface, updating `docs/v1.0-surface.md` to
   record the deferral and the reason; **(B)** extend the phase
   with a kernel specific reformulation slice. Adding a faithful
   qualifier to the public claim is off the table. This rule
   binds the phase's family slices uniformly; per family
   AskUserQuestion checkpoints surface intractability findings as
   forks for the project lead, not as defaults to take.

5. **Phase 1f closes before slice 8c opens.** The pf-g8h bead
   (v1.0 version bump) gains a `→` dependency on every Phase 1f
   slice's bead. ADR-0033's slice 8c parking persists through
   Phase 1f.

## Consequences

**Honest framing.** The v1.0 the project publishes carries the
strongest claim a pure Rust arbitrary precision library has shipped:
correctly rounded across every IEEE 754-2019 mode on the entire
unary v1.0 surface, exhaustively f32 verified against three
independent oracles per ADR-0035. The README's existing
"Verification posture" prose stays exactly true. The two tier reality
in DESIGN.md "Caveats and open questions" §1 deletes at the end of
the phase (slice p1.37) because §1's content becomes empty as
families migrate.

**Timeline cost, accepted.** The phase plan at
`~/.claude/plans/phase-1f-dynamic-fog.md` estimates 3 to 5 months
of active session time across 16 slices (p1.22 through p1.37).
Closer to 3 if the audit's drop in identifications hold; closer to
5 if multiple families need per kernel reformulation. The cost is
borne by the project; pfloat has no published 0.x consumers, so the
delay does not propagate to a downstream migration tax. ADR-0033's
"published 1.0 is immutable" argument carries through: months of
delay against an immutable strong claim beats prompt publication of
a half story.

**v1.0 ship criterion sharpens.** ADR-0033 fixed the criterion at
"every row reads `correctly-rounded` or `faithful`, no `has-errors`".
This ADR removes the `faithful` allowance: every row reads
`correctly-rounded` under every mode. The change is durable;
subsequent releases inherit the standard.

**Per mode status TOML schema migration.** The existing schema
(`rounding_modes = "RNE"`, `rounding_status = "correctly-rounded"`)
migrates at slice p1.23 to a per mode table. Every row in
`tests/oracle/status/` (44 today, more after the parametric Bessel
expansion at slice p1.35) gains a per mode verdict. The schema is
the durable artifact; the transient sweep output is not.

**Audit derived migration plan, not recall.** Slice p1.22 produces
the in tree audit document at
`docs/decisions/plans/phase-1f-five-mode-completeness.md`, which
records per kernel: the `eval(w)` shape, the cancellation regimes,
the per regime Ziv strategy (drop in wrap, cancellation short
circuit, reformulation, basis change), the cited spec section for
any reformulation, the oracle coverage, the estimated Ziv
iterations at cap precision, and a worked numeric example at a
known hard input. The audit is the load bearing artifact; the per
family slices execute against it. The pf-jn1y misdiagnosis
precedent (`feedback_kernel_vs_harness_diagnosis`) and the 8a
case-4 O(m) DoS precedent
(`feedback_derive_dont_recall_coefficients` instance 6) anchor the
discipline: every kernel's strategy is derived from the kernel's
source, not recalled from a similar kernel's migration pattern.

**Existing scaffolding beads supersede.** pf-lw3l (differential
lane widening on the already five mode cohort) and pf-fwtz
(`certified_round_bf_to_f32` helper) carry the load bearing
scaffolding scope and close at slice p1.22 with their work folded
into slice p1.23 (Phase 1f's scaffolding slice). pf-xyaq (CI matrix
widening, the surviving slice 8c blocker pre Phase 1f) stays
separate and lands before or during the phase.

**Permacomputing horizon unchanged.** The phase adds no runtime
dependency. The Ziv driver, the bf to f32 mode aware helper, and
the differential MPFR oracle synthesis run inside the existing test
infrastructure. Per `feedback_descope_saas_on_permacomputing_projects`
no new SaaS integration enters the verification stack. The Arb
worker stays a `python-flint` subprocess (LGPL out of link graph);
mpmath remains pure Python BSD; Maxima remains invoked through
`nix-shell`. pfloat retains its zero runtime dependency posture per
ADR-0037.

**Per push CI gate stays NE only and Python free.** Per ADR-0035
the per push gate runs MPFR only and avoids Python entirely; this
stays. The five mode sweep runs at per release cadence on the
existing differential MPFR Linux lane. Per slice gates exercise the
Arb worker (the existing pattern). No expansion of per push
compute cost during the phase.

**Disclosure block stays untouched.** Per
`feedback_disclosure_update_under_explicit_permission` the README
"How pfloat is developed" disclosure block carries protected
invariants that bind regardless of phase work. No Phase 1f slice
edits the block. The diff regeneration at slice p1.37 is context
line refresh only; the block content stays bit identical.

**Amendment at Phase 1g closure (2026-05-26, slice p1g.5).**
ADR-0038's correctness claim is tightened by ADR-0039 (Phase 1g
closure): the `ZIV_ERROR_GUARD` assumption is calibrated per
kernel (`pf-yupm`, central table in `src/math/ziv_calibration.rs`),
actively guarded by the per-release oracle-sweep cross-check
(`pf-tqzz`, `tests/oracle_cross_check_smoke.rs` + the full
release sweep), and the driver's interval-test soundness theorem
is formally stated and Kani-scaffolded over the canonical
operand-bounding pattern (`pf-hdh8`,
`src/verify/ziv_soundness.rs`). The `BoundedBigFloat<80>`
fixed-array encoding for the universal-quantification proof
discharge is recorded as post-v1.0 follow-up in ADR-0039. Phase 1f's
no-narrowing-of-the-claim posture stays exactly intact; Phase 1g
adds the verification-architecture closure layer beneath it.

**Family ordering accounts for inter family dependencies.** The
audit's recommended ordering puts forward trig (slice p1.26) before
the gamma family (slice p1.29) so the gamma reflection through
`sin` lands on a five mode correct `sin`. Zeta (slice p1.34) runs
last because its functional equation branch composes `sin`, `pow`,
and `gamma`. The audit at p1.22 may revise the ordering after eyes
on kernel source review; per family AskUserQuestion checkpoints
surface any ordering changes before code is written.

## Alternatives considered

**Add a "faithful at directed mode ties" qualifier to the
README.** Documents the actual implementation state honestly; ships
v1.0 promptly. Rejected on the CLAUDE.md frugality argument: the
qualifier consumes the canonical name on a half story. Downstream
users who needed five mode rigor on `gamma` or `bessel_y` would
have to vendor pfloat, route around it, or wait for v2. The
qualifier also normalises the cohort labelling pattern, against
which CLAUDE.md is explicit. The honest framing of "implementation
state is two tier today" is exactly the position v1.0 should not
ship under.

**Narrow the v1.0 surface preemptively, dropping the kernels that
would need migration.** Closes the gap by shrinking the
implementation surface to match the existing five mode cohort.
Rejected because the affected functions (trig, hyperbolic, gamma
family, Bessel Y/I/K, Airy, integrals, zeta) are exactly the
high consumption use cases pfloat is built for; dropping them
makes v1.0 a smaller and less useful library than the existing 0.1.
The pre v1.0 surface freeze at slice p1.3 captured the right shape;
the implementation catches up rather than the surface narrowing.

**Ship v1.0 with the existing two tier reality, document it
clearly, plan v2 as the five mode pass.** Rejected for the same
reason as the qualifier alternative: the published 1.0 is what
users judge the project against permanently. ADR-0033 already
absorbed timeline cost on the same argument; doing so again for
directed mode rigor is consistent.

**Run a partial migration: arithmetic core plus the elementary
transcendentals only, leave specials at NearestEven only.**
Reduces phase wall clock but ships v1.0 with the half story
qualifier the no qualifier rule rejects. The qualifier would have
to live somewhere, and the cohort labelling failure mode CLAUDE.md
names applies wherever it lands. Rejected.

## References

- `~/.claude/plans/phase-1f-dynamic-fog.md` — the phase plan
  approved 2026-05-25 that this ADR ratifies.
- `docs/decisions/plans/phase-1f-five-mode-completeness.md` — the
  in tree audit document produced at slice p1.22, the load bearing
  per kernel derivation work.
- ADR-0022 — the Ziv interval test driver, the structural
  prerequisite this phase consumes unchanged.
- ADR-0033 — Phase 1 sweep precedes v1.0; the prior credibility
  cost argument this ADR extends.
- ADR-0034 — Oracle layer; the per mode TOML schema migration
  refines ADR-0034's verdict surface.
- ADR-0035 — Oracle worker protocol and three way agreement; the
  protocol extends to per mode verification at per release cadence.
- ADR-0036 — `property_jn` dyadic constraint; the per family
  property tests inherit the dyadic discipline uniformly.
- ADR-0037 — `SmallVec` inline storage rejected; the zero runtime
  deps posture binds Phase 1f.
- `DESIGN.md` "Caveats and open questions" §1 — the paragraph this
  phase deletes by end of slice p1.34.
- `docs/v1.0-surface.md` — the frozen 47 entry surface this phase
  brings into uniform five mode correctness.
- `src/math/ziv.rs` — the driver this phase wires every remaining
  kernel through.
- `tests/differential/mod.rs` `BIT_EXACT_ROUNDING_MODES` /
  `NEAREST_EVEN_ROUNDING_MODES` — the rounding mode constants
  this phase migrates lanes between.
- `feedback_ziv_interval_test_and_mpfr_rnda` — cancellation to
  zero defeats the interval test; the audit identifies every
  kernel that needs a short circuit before the Ziv driver sees it.
- `feedback_bf_to_f32_directed_mode` — the Display+parse output
  bridge is silently NE only at p>24; slice p1.23's
  `certified_round_bf_to_f32` lifts the constraint.
- `feedback_kernel_vs_harness_diagnosis` — multi precision probe
  before kernel attribution on every sweep `has-errors` finding.
- `feedback_derive_dont_recall_coefficients` — every Ziv strategy,
  cancellation threshold, and complexity bound derived from kernel
  source and cited spec, not recalled.

## Commits

- This ADR, the audit document, the ROADMAP.md Phase 1.5 entry, and
  the bead surgery summary land as one signed merge of unsigned
  branch commits on `slice-p1-22-five-mode-completeness-audit`.
  The implementation work (scaffolding at slice p1.23, per family
  migrations across slices p1.24 through p1.34, parametric Bessel
  sweep at p1.35, multi arg confirmation at p1.36, prose alignment
  at p1.37) follows in subsequent slices as the audit's per family
  strategies execute.
