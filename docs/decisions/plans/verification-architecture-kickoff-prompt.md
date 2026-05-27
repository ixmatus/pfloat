# Kickoff prompt — verification architecture (Option 3 + Kani)

Paste this into a fresh session to resume the Option 3 + Kani
verification-architecture work for pfloat. The four beads
referenced below are the deliverable.

---

I'm continuing the verification-architecture work we discussed
after Phase 1f closed (merge `b08a831` on `main`, 2026-05-26).
Phase 1f made every v1.0-surface kernel correctly rounded under
all five IEEE 754-2019 rounding modes; the remaining gap is the
honest one named in `DESIGN.md` § Caveats: the runtime
`ZIV_ERROR_GUARD = 24` in `src/math/ziv.rs` is an assumed
internal-error bound, not a proven one.

We talked through four paths to close that gap. The plan is
**Option 3 + Kani**, scoped across four beads:

1. `pf-yupm` — per-function `ZIV_ERROR_GUARD` calibration (Option 3,
   part 1): replace the global 24 with a per-kernel tunable bound,
   pinned by error analysis or widened sweep at each kernel's
   worst-case input.
2. `pf-tqzz` — Arb cross-check assertion in oracle sweep (Option 3,
   part 2): plumb the working-precision intermediate out of
   `ziv_round` and assert it stays within the per-function bound on
   every f32 sweep input. Converts the existing oracle sweep into
   an active guard on the assumption.
3. `pf-hdh8` — Kani-discharge of Ziv interval-test soundness
   (Option 4, the formal-methods piece): state and prove the
   theorem `∀ y, half_width, mode: interval_test(y, half_width, t,
   mode) ⇒ ∀ y' ∈ [y ± half_width]: round(y', t, mode) =
   round(y, t, mode)`. Discharge once for the driver; per-kernel
   side becomes a stated hypothesis justified by `pf-yupm` +
   `pf-tqzz`.
4. `pf-kk16` — exact-value pre-Ziv dispatch pattern audit: enumerate
   the inputs across every v1.0-surface kernel where the true value
   is exactly representable in target precision (the gamma(n)
   defect class surfaced in slice p1.29). Each non-empty subset
   gets a pre-Ziv exact dispatch.

**Suggested ordering** (none of this is strictly forced; you may
re-sequence based on what surfaces during the work):

- `pf-kk16` first as a mechanical audit pass. It bounds the
  "exactly-representable-true-value" surface so subsequent error
  analysis doesn't have to special-case those.
- `pf-yupm` and `pf-tqzz` together as the Option 3 pair. `pf-yupm`
  is the calibration; `pf-tqzz` is the runtime guard that catches
  the assumption being violated on swept inputs. They feed each
  other: calibration sets the per-function bound, the cross-check
  guards it.
- `pf-hdh8` as a separate longer-arc workstream. Kani works on
  bounded models; the proof requires a bounded-`BigFloat`
  encoding (limbs, precision, exponent ranges) covering the
  operational envelope. Independent of the other three at the
  proof level but consumes their work as the kernel-side
  hypothesis.

**Bead descriptions** carry "DEFER until pf-g8h (v1.0 ship)"
framing from when I filed them. That framing reflects my
assumption about scheduling, not your direction — ignore the
defer language and treat the scheduling as open.

**Durable in-tree references** (no need for me to re-read anything
huge; these are the anchors):

- `src/math/ziv.rs` — the driver. `ZIV_ERROR_GUARD = 24`
  constant; `ziv_round(eval, target, mode)` signature; interval
  test via half-width `|y|·2^-(working-24)`; iteration cap 5.
- `tests/oracle/pfloat_kernels.rs` — per-FnId verification
  precision dispatch.
- `tests/oracle/arb.rs` — Arb worker subprocess wiring for
  `pf-tqzz`.
- `docs/decisions/0038-five-mode-kernel-completeness-as-v1.0-gate.md`
  — ADR-0038, Phase 1f's load-bearing decision.
- `docs/decisions/plans/phase-1f-five-mode-completeness.md` — the
  audit doc with per-kernel derivations.
- `MEMORY.md` entries `feedback_exact_value_defeats_ziv`,
  `feedback_ziv_interval_test_and_mpfr_rnda`,
  `feedback_irrational_constant_special_case_mode_aware`.

**Where to start a fresh session.** Run `bd ready` to confirm the
four bead IDs and check no one else has claimed them. Then `bd
update <bead-id> --claim` on the first one and start.

The same long-arc workflow that worked for Phase 1f applies here
if any of these slices want it (single phase branch, unsigned
commits, pause-to-debrief at slice boundaries, one signed merge at
phase closure). Otherwise the standard per-slice signed-merge
cadence is fine; each bead's scope is small enough to fit one
slice.
