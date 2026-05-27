# ADR-0039: Phase 1g verification architecture closure (v1.0 blocker)

- **Status**: accepted; Phase 1g closed at slice p1g.5 (this
  amendment + ADR-0038 amendment + DESIGN.md Caveats §1 narrowing
  + README Verification posture tightening). The four beads
  (`pf-kk16`, `pf-yupm`, `pf-tqzz`, `pf-hdh8`) all closed; the
  v1.0 ship gate is satisfied modulo the pre-existing `pf-xyaq`
  CI gate and the slice-8c release ceremony. Honest framing for
  what landed vs deferred: per-bead details below.
- **Date**: 2026-05-26

## Context

Phase 1f closed at merge `b08a831` (2026-05-26), making every v1.0-
surface kernel correctly rounded under every IEEE 754-2019
rounding mode. ADR-0038 is the load-bearing decision. The strong
claim "correctly rounded across every IEEE 754-2019 mode on the
entire unary v1.0 surface" is the claim v1.0 will ship under.

The remaining gap is named in `DESIGN.md` "Caveats and open
questions" §1 and in the doc comment at `src/math/ziv.rs:50-58`:
the constant `ZIV_ERROR_GUARD = 24` at `src/math/ziv.rs:59` is an
**assumed** internal-error bound, not a **proven** one. The
half-width formula `|y|·2^-(working - 24)` is sound only if every
kernel's `eval(w)` returns a value within `2^24` ULPs of the true
value at working precision `w`. The doc comment cites the
empirical analysis ("`pow_int`'s square-and-multiply uses ≤ 64
multiplies, the `exp·ln` path a handful of operations, all far
under `2^24` ULP") but does not derive the bound per kernel or
guard it actively at sweep time.

Three load-bearing gaps follow:

1. **No per-kernel calibration.** The single global value is
   conservative for most kernels but may be loose for some
   (wasting Ziv iterations) and tight for others (lgamma near
   negative-half-integer poles, `bessel_y` near zeros, `zeta` in
   the critical strip). Each kernel's actual bound depends on its
   `eval(w)` op count and cancellation regime; there is no audit
   that derives the bound from the kernel's source.

2. **No active sweep-time guard.** The oracle sweep certifies
   correct rounding against the rigorous oracle, but it does not
   verify that pfloat's `eval(w)` intermediate stayed within the
   assumed bound. A kernel whose internal error silently exceeded
   `2^24` ULPs would pass the correct-rounding check at sweep
   time only by accident (the interval test happened to fall on
   the right side of the boundary), then fail in production on
   nearby inputs.

3. **No soundness proof of the interval test itself.** The Ziv
   driver assumes that if `round(y − h, t, m) == round(y + h, t,
   m)` then every value in `[y ± h]` rounds the same way under
   mode `m` at target precision `t`. This is the soundness
   property the driver relies on but never proves; the property
   is true (it follows from monotonicity of round-to-precision
   under each mode) but the proof lives in design intuition, not
   in a Kani-discharged theorem.

A fourth gap surfaced during Phase 1f (slice p1.29, commit
`ec89ffd`): **exactly-representable true values defeat Ziv's
interval test under directed modes.** When the kernel's
composition returns `T(x) + epsilon` where `T(x)` is exactly
representable at target precision and `epsilon` is the
working-precision noise floor, the Ziv interval always spans the
rounding boundary at `T(x)`, NE rounds correctly through epsilon,
and directed modes tip rounding to the adjacent ULP. `gamma(7)`
at `p=53` mode `TowardPositive` returned `720.00000000000011`
(1 ULP high) instead of exact `720 = 6!`. The fix landed for
gamma in slice p1.29 as a pre-Ziv exact dispatch
(`try_gamma_pos_integer_exact`). The defect class is transferable;
the remaining v1.0-surface kernels have not been audited for it.

ADR-0033 sequenced Phase 1 before the v1.0 tag on the credibility
argument that "published 1.0 is immutable". ADR-0038 extended the
same sequencing posture for directed-mode rounding. This ADR
extends it again for the verification-architecture gap: closing
the four gaps before v1.0 ensures the published 1.0 carries the
strongest verification posture the project can offer, not a
posture that quietly tightens at v1.1 and leaves v1.0 consumers
linked against the weaker claim.

CLAUDE.md's verification-as-entropy-reduction principle is
explicit:

> Verification is entropy reduction. Every fact about behavior
> not captured in types or proofs gets rediscovered by every
> future user, every reviewer, every 3am debugging session. Order
> of preference, most frugal to least: types, proofs, property
> tests, example tests, documentation, nothing. Reach for the
> highest level the artifact's durability warrants. A Kani
> discharged property holds for every future input; a fuzzer only
> finds the inputs you tried.

`ZIV_ERROR_GUARD = 24` as a documented assumption lives at the
"documentation" tier in that ordering. Phase 1g moves the four
gaps up the ordering: per-kernel calibration (audit-derived
documentation), active sweep guard (property test across the f32
grid), driver soundness theorem (Kani-discharged proof at fixed
IEEE target precisions), exact-value pattern audit (per-kernel
pinning tests at directed modes).

## Decision

Phase 1g runs to completion before the v1.0 tag. The v1.0 surface
ships when all four beads close and the phase-closure prose lands.
The closure status at phase merge:

- **pf-kk16** — exact-value pre-Ziv dispatch pattern audit across
  every v1.0-surface kernel; per-kernel pre-Ziv dispatch where
  the exact-value subset is non-empty.
- **pf-yupm** — per-function `ZIV_ERROR_GUARD` calibration; the
  global `24` constant becomes a `DEFAULT_ERROR_GUARD = 24` in a
  new central `src/math/ziv_calibration.rs` table; every kernel
  calls `ziv_round` with an explicitly-cited per-function bound;
  a structural enumeration test forces conscious calibration at
  every five-mode-correct `FnId`.
- **pf-tqzz** — Arb cross-check assertion in the oracle sweep
  asserting `|pfloat_eval(w) − arb_midpoint| ≤ 2^(error_guard −
  w) · |arb_midpoint|` for every `(kernel, input)` pair across
  the f32 grid; runs per-release, not per-push.
- **pf-hdh8** — Kani-scaffolded discharge of the Ziv interval-
  test soundness theorem at fixed target precisions `t ∈ {24,
  53, 113}`. **Closure note (revised):** the scaffolding lands
  as four `#[kani::proof]` harnesses in
  `src/verify/ziv_soundness.rs` using the canonical eight-
  constant operand-bounding pattern shared with the existing 196
  harnesses (ADR-0012). Local `cargo kani` discharge of the
  simplest harness timed out at >11 minutes CBMC runtime,
  matching the pre-existing transcendental-harness runtime
  profile (ADR-0012 slice-6k status update). The original
  `BoundedBigFloat<80>` fixed-array encoding proposed in the
  earlier audit doc is the right shape for tractable CBMC
  symbolic execution; building it (plus conversion shims and
  per-operation soundness lemmas) does not fit the v1.0-ship
  budget and lands as a post-v1.0 follow-up.

The strategic commitments are non-negotiable:

1. **Pre-committed Kani scope, honestly partial at phase closure.**
   The soundness theorem scaffolding lands at fixed `t ∈ {24, 53,
   113}` over the canonical eight-constant operand-bounding
   pattern. The original audit-doc design called for a
   `BoundedBigFloat<80>` fixed-array shadow type to make CBMC's
   symbolic execution tractable; that encoding is documented as a
   post-v1.0 follow-up after the local discharge run timed out at
   >11 minutes on the simplest harness (matching the pre-existing
   transcendental-harness runtime profile per ADR-0012). The
   structural-analogy claim for the arbitrary-mantissa surface
   stays as written: the round-to-precision predicate is uniform
   across mantissa values within a class, the interval test is
   uniform across mantissa values within a class, and the
   canonical-class scaffolding stands in for the family. The
   actual Kani discharge is the v1.x scope; v1.0 ships with the
   scaffolding plus the pf-tqzz sweep cross-check actively
   guarding the kernel-side bound at every f32 input.

2. **Conscious calibration, not implicit defaults.** Every kernel
   passes its `error_guard` explicitly to `ziv_round`. No
   implicit `24` default at the call site. The
   `KNOWN_CALIBRATED_KERNELS` table is enumerated by a
   structural test that fails on any unenumerated five-mode-
   correct `FnId`. Adding a new kernel to the v1.0 surface
   without calibrating it fails the test.

3. **Cross-check at per-release cadence, not per-push.** The
   ~3.1M additional Arb calls do not fit a per-push budget; they
   do fit a per-release budget (hours, not minutes). The
   per-push gate stays NE-only and Python-free per ADR-0035.

4. **Disclosure block stays untouched; Verification posture
   tightens.** The README "How pfloat is developed" disclosure
   block (per `feedback_disclosure_update_under_explicit_
   permission`) stays bit-identical. The Verification posture
   section (lines 189–194) tightens to reference ADR-0039,
   the per-release oracle-sweep cross-check, and the Kani-
   discharged soundness theorem at IEEE binary32/64/128. `pf-
   epf`'s `docs/disclosure-correction-v1.0.diff` continues to
   apply cleanly after the Verification-posture edit (verified
   via `git apply --check`).

5. **Phase 1g closes before slice 8c opens.** `pf-g8h` (v1.0
   version bump) gains a blocking dependency on every Phase 1g
   bead. `pf-4fi` (crates.io check) and `pf-xyaq` (CI gate) are
   independent and may proceed at any time. ADR-0033's slice 8c
   parking persists through Phase 1g.

## Consequences

**Honest framing.** The v1.0 the project publishes carries the
strongest verification posture a pure-Rust arbitrary-precision
library has shipped: every kernel correctly rounded across every
IEEE 754-2019 mode on the entire unary v1.0 surface (ADR-0038);
per-kernel internal-error bounds calibrated rather than assumed
and actively guarded by the per-release oracle-sweep cross-check
(pf-yupm + pf-tqzz); exact-value pre-Ziv dispatch audited across
every kernel (pf-kk16); driver interval-test soundness discharged
by Kani at IEEE binary32/64/128 target precisions over a bounded
`BigFloat` encoding (pf-hdh8). The Caveats §1 paragraph in
DESIGN.md narrows to retain only the measure-zero termination
caveat (MPFR-parallel); the empirical-slack framing migrates into
a cross-reference to `ziv_calibration.rs`.

**Timeline cost, accepted.** Phase 1g estimates 6 to 10 weeks of
active session time across 5 slices (p1g.1 through p1g.5).
Closer to 6 if `pf-hdh8`'s `BoundedBigFloat<80>` discharges cleanly
at fixed `t`; closer to 10 if it requires the contingency
escalations recorded in the plan doc. The cost is borne by the
project; pfloat has no published 0.x consumers, so the delay
does not propagate to a downstream migration tax. ADR-0033's
"published 1.0 is immutable" argument carries through.

**v1.0 ship criterion sharpens again.** ADR-0033 fixed the
criterion at "every row reads `correctly-rounded` or
`faithful`, no `has-errors`". ADR-0038 removed the `faithful`
allowance: every row reads `correctly-rounded` under every mode.
This ADR adds the per-kernel calibration and Kani-discharged
soundness criteria: every kernel calls `ziv_round` with an
explicitly-cited per-function `error_guard`, the cross-check
assertion passes across the full f32 grid, and the Kani
soundness harnesses complete at `t ∈ {24, 53, 113}` with no
counterexample. The change is durable; subsequent releases
inherit the standard.

**Per-kernel hypothesis explicitly stated.** The Kani-discharged
soundness theorem moves the universal-quantifier proof obligation
off every kernel and onto the driver, but it admits a per-kernel
hypothesis: "for all inputs the kernel sees, its `eval(w)`
returns a `BigFloat` whose value, lifted into
`BoundedBigFloat<80>`, satisfies the working-precision bound."
The per-kernel hypothesis is what `pf-yupm` (calibration) and
`pf-tqzz` (sweep guard) jointly justify. The three pieces
combine into one verification chain: calibration → guarded →
proved.

**Arbitrary-precision claim, structural analogy.** The Kani
soundness theorem discharges at `t ∈ {24, 53, 113}` only. For
arbitrary `t` (e.g., `t = 200`, `t = 1024`), the soundness
property holds by structural analogy: the round-to-precision
function is uniform across `t`, the interval test is uniform
across `t`, and the bounded encoding scales without changing
shape. The structural analogy is honest framing, not a proof;
the README Verification posture and ADR-0039 both say so.

**Permacomputing horizon unchanged.** Phase 1g adds no runtime
dependency. The Kani toolchain is already wired (ADR-0012;
manual on-demand workflow). The Arb worker stays a `python-flint`
subprocess. No new SaaS integration enters the verification
stack (`feedback_descope_saas_on_permacomputing_projects`).
pfloat retains its zero-runtime-dependency posture (ADR-0037).

**Per-push CI gate unchanged.** The per-push gate stays NE-only
and Python-free (ADR-0035). The Phase 1g cross-check runs at
per-release cadence on the existing differential MPFR Linux
lane; per-slice gates exercise the Arb worker (existing pattern).
No expansion of per-push compute cost.

**Disclosure block stays bit-identical.** The README "How pfloat
is developed" disclosure block carries protected invariants
(`feedback_disclosure_update_under_explicit_permission`).
Phase 1g does not edit the block. The Verification posture
section tightens; `pf-epf`'s disclosure-correction diff applies
cleanly after the tightening (verified via `git apply --check`
at slice p1g.5).

**Bead-graph surgery at phase entry.** `pf-kk16` opens the chain
(no blockers). `pf-yupm` blocks on `pf-kk16`. `pf-tqzz` and
`pf-hdh8` both block on `pf-yupm` (the calibration table is
their input). `pf-g8h` blocks on all four. `pf-4fi` and `pf-xyaq`
stay independent. The previous "DEFER until pf-g8h" framing in
each bead description is superseded by the Phase 1g framing
header.

## Alternatives considered

**Punt Phase 1g to v1.x; ship v1.0 on Phase 1f only.** Documents
the ZIV_ERROR_GUARD assumption honestly in DESIGN.md; ships v1.0
promptly with the existing "extensively validated as correctly
rounded" prose; closes the four gaps at v1.x or v2.x. Rejected on
the CLAUDE.md frugality argument: v1.0 consumes the canonical
name; tightening the verification posture at v1.x leaves v1.0
consumers linked against the weaker claim. The ADR-0033 +
ADR-0038 precedent of paying the timeline cost upfront to keep
the published version's claim load-bearing applies here too. The
gap named in DESIGN.md Caveats §1 is exactly the kind of "future
maintenance tax" CLAUDE.md's frugality principle rejects.

**Spike-and-fallback Kani encoding scope.** Originally
recommended in the plan-drafting Plan agent's output: try the
full encoding (symbolic `t`), fall back to fixed `t` if the
Kani search space proves intractable. Rejected at the user's
direction in favor of pre-committing the scope: front-loading the
scope decision is sharper, matches ADR-0038's drop-or-extend
posture, and avoids the spike-and-pivot pattern that consumes
timebox without producing a deliverable. The fixed-`t` scope is
the deliverable from the start.

**Per-slice signed-merge cadence instead of long-arc.** Each
bead as its own branch with a separate signed merge into main,
ADR amended or supplemented per slice. Rejected because three of
the four beads (pf-kk16, pf-yupm, pf-tqzz) share a load-bearing
artifact (the per-kernel calibration table). Per-slice signed
merges would force three separate ADR-amendment cycles and three
separate disclosure-prose touchpoints. The Phase 1f long-arc
precedent (`feedback_phase1f_long_arc_workflow`) fits the
coupling better.

**Companion `ziv_round_with_guard()` function instead of
extending `ziv_round` signature.** Considered for pf-yupm:
preserve the existing `ziv_round` signature (no `error_guard`
parameter, defaults to 24) and add a separate
`ziv_round_with_guard(eval, target, mode, error_guard)` for
callers that override. Rejected: doubles the driver's surface
area; creates two paths to keep in sync. `ziv_round` is
`pub(crate)`, so the breaking change is contained to the 44 call
sites under `src/math/`, all updated in the same slice. The
acceptance-criterion-4 test ("per-function table non-empty for
every five-mode-correct kernel") is much easier to enforce
structurally when every call site is forced to pass an explicit
constant.

## References

- `docs/decisions/plans/phase-1g-verification-closure.md` — the
  in-tree audit document this ADR ratifies; populated
  incrementally across the four implementation slices.
- `docs/decisions/plans/verification-architecture-kickoff-prompt.md` —
  the prior-session kickoff prompt this ADR derives from.
- ADR-0022 — the Ziv interval-test driver; the structural piece
  Phase 1g calibrates, guards, and proves.
- ADR-0033 — Phase 1 sweep precedes v1.0; the credibility-cost
  precedent.
- ADR-0034 — Oracle layer.
- ADR-0035 — Oracle worker protocol; the MIDPOINT verb extension
  for pf-tqzz adds one verb without changing the per-push
  posture.
- ADR-0038 — Five-mode kernel completeness as the v1.0 gate;
  amended at slice p1g.5 with a Consequences paragraph
  cross-referencing this ADR.
- ADR-0012 — Kani harness scaffolding and the CBMC `Vec<u64>`
  limitation that motivates the `BoundedBigFloat<80>` fixed-
  array encoding choice for pf-hdh8.
- `src/math/ziv.rs` — the driver this phase calibrates, extends,
  and proves.
- `DESIGN.md` "Caveats and open questions" §1 — the paragraph
  this phase narrows at slice p1g.5.
- `feedback_exact_value_defeats_ziv` — the slice p1.29 defect-
  class lesson pf-kk16 audits across the full surface.
- `feedback_ziv_interval_test_and_mpfr_rnda` — cancellation-to-
  zero defeats the interval test; the audit identifies every
  kernel that needs a short-circuit before the Ziv driver sees
  it.
- `feedback_phase1f_long_arc_workflow` — the long-arc workflow
  Phase 1g inherits.
- `feedback_disclosure_update_under_explicit_permission` — the
  protected disclosure block invariants.

## Commits

- This ADR, the phase-1g plan doc skeleton, and the four-bead
  surgery (description rewrites + bead-graph edges) land as one
  signed merge of unsigned branch commits on
  `phase-1g-verification-closure` at phase closure (slice p1g.5).
  The implementation work (`pf-kk16` at p1g.1, `pf-yupm` at
  p1g.2, `pf-tqzz` at p1g.3, `pf-hdh8` at p1g.4) follows in
  subsequent slices as the audit's per-bead strategies execute.
