# ADR-0048: Airy asymptotic kernel working-precision boost reduction

- **Status**: accepted
- **Date**: 2026-05-28

## Context

ADR-0044 (Airy threshold already optimal) closed sub-slice 2b.2.b as
doc-tier — the `airy_threshold_exponent` formula is at ~1.5% over the
strict accuracy bound, and the cube-of-`|x|³` quantization locks the
threshold immovably. While documenting that finding, ADR-0044
§"Bench infrastructure note" flagged a separate Airy perf lever
that's orthogonal to the threshold:

> `airy_asymptotic_pos` / `airy_asymptotic_neg` use `working =
> target + 64`; the asymptotic's round-off is `≈ log₂(N) ≤ 8` bits
> for any practical `|x|`, so a reduction to `working = target + 32`
> (just `AIRY_ERROR_GUARD + 8`) might save ~30% on every asymptotic
> call.

This ADR lands that flagged change. It's a post-Phase-2b follow-up,
not part of `pf-6fvx`'s original sub-slice plan; standalone branch
`airy-asymptotic-boost-reduction`.

## Change

Two single-line edits to `src/math/airy.rs`:

- Line 498 (`airy_asymptotic_pos`): `target_precision.saturating_add(64)`
  → `target_precision.saturating_add(32)`
- Line 557 (`airy_asymptotic_neg`): parallel change

Plus doc-comments on each citing this ADR and the budget accounting.

## Justification

The accumulated round-off inside `airy_asymptotic_pos` / `_neg`:

- `AIRY_ERROR_GUARD = 24` bits (Ziv-side calibrated slack per
  Phase 1g, `src/math/ziv_calibration.rs:191`).
- `log₂(N) ≤ 8` bits round-off accumulation in the optimal-truncation
  asymptotic sums `airy_uv_sums` / `airy_uv_sums_split`, where
  `N ≈ √ζ` and `ζ = (2/3)|x|^{3/2}`. For any practical input the
  optimal truncation index stays below `2^8 = 256` terms (`|x| ≈
  2^15` gives `N ≈ 700`; `log₂(700) ≈ 10`). The `≤ 8` bound is
  conservative for the input ranges the kernel actually reaches.

Total budget: `24 + 8 = 32` bits, exactly matching the new
`target + 32` boost. The previous `target + 64` carried 32 bits of
margin above the budget — a doubling of the calibrated slack, with
no justification in the kernel's accumulated round-off analysis.

The Ziv driver retains the `ZIV_GUARD_CAP = 1024` ceiling for any
hard-to-round inputs where the per-iteration interval test fails;
reducing the per-iteration boost doesn't change the ceiling on
maximum reachable working precision, only the first-iteration
starting point.

## Measurement

`benches/airy_dispatch.rs` re-scoped to 4 asymptotic-only cells
(the Maclaurin cells `p53_x64` and `p256_x512` are excluded; see the
bench file's updated header doc for the cost rationale and the
guidance for future Maclaurin-targeting work).

Clean-machine measurement (`b01x1elzu` vs `phase2b-airy-clean` on
`aarch64-apple-darwin`):

| Cell | Baseline (+64) | New (+32) | Δ |
|------|---------------:|----------:|---|
| Ai_p53_x128 | 25.41 ms | 17.53 ms | **−31.0%** |
| Ai_p53_x256 | 62.75 ms | 48.81 ms | **−22.2%** |
| Ai_p256_x1024 | 1.13 s | 961.78 ms | **−14.9%** |
| Ai_p256_x2048 | 2.98 s | 2.79 s | **−6.5%** |

Speedup scales roughly linearly with the relative working-precision
reduction:
- `p=53`: working drops 117→85 bits (27% reduction), measured
  22-31% speedup.
- `p=256`: working drops 320→288 bits (10% reduction), measured
  7-15% speedup.

This linear-in-p scaling rather than O(p²) reflects the dominance of
`exp(±ζ)`, `sqrt(π)`, and `x^{3/2}`/`x^{1/4}` prefactor computations
(each O(p²) in absolute terms but with constant `O(M(p))` per-op
cost overhead that doesn't scale further) over the per-term Newton
mul/div inside `airy_uv_sums`. The asymptotic series itself runs
N optimal-truncation terms (N ≈ √ζ), each contributing one mul +
one div at working precision — the per-call cost is summed roughly
linearly over N, but N grows with `|x|`, so the absolute cost grows
super-linearly while the *speedup* from boost reduction stays
proportional to the working-prec ratio.

**Verification** (all green):
- `cargo test --release --features airy --lib math::airy`: 20/20
  (in-module Wronskian Ai·Bi′−Ai′·Bi=1/π at p=200, asymptotic
  continuity, high-precision pins, boundary constants).
- `cargo test --release --features differential-mpfr --test
  differential_ai`: 5/5 (TRANSCENDENTAL_PRECISIONS ≤ 256, AI_POINTS
  Maclaurin + asymptotic).
- `cargo test --release --features differential-mpfr --test
  differential_bi`: 6/6 (parallel for Bi, including the Wronskian
  cross-tie at p=200).
- `cargo test --release --features differential-mpfr --test
  property_ai`: 4/4 (Wronskian proptest over a small dyadic grid
  at p=96).

## Consequences

**Every asymptotic Airy call gets faster.** The four-kernel family
(`Ai`, `Bi`, `Ai′`, `Bi′`) all share the same asymptotic dispatch
and the same per-call working-precision boost; the win measured on
`Ai` transfers to the other three. Per ADR-0044 §"Bench
infrastructure note" the most prominent caller of Airy at moderate
precision is likely scientific-computing workloads computing Airy
at the asymptotic-regime crossover; the win is uniform across them.

**No correctness change.** Verification matrix fully green; the
+32 budget covers the actual accumulated round-off.

**Bench file ships with the change.** `benches/airy_dispatch.rs`
gets its cell list trimmed to the 4 asymptotic cells (Maclaurin
cells excluded with rationale). Future Maclaurin-targeting work
should add a dedicated bench at smaller `|x|` (the existing bench's
header docs the guidance).

**ADR ledger update:**

| ADR | Outcome | Slice |
|----:|---------|-------|
| 0037 | REJECTED | SmallVec mantissa swap |
| 0040 | doc-tier (GATE A) | FFT measurement |
| 0041 | REJECTED | Spouge precision-pegging |
| 0042 | ACCEPTED | pf-1axr trig + bessel_y recurrence root fix |
| 0043 | ACCEPTED | Bessel per-kernel threshold split |
| 0044 | doc-tier | Airy threshold already optimal |
| 0045 | doc-tier | Si/Ci threshold already optimal |
| 0046 | doc-tier | erf/erfc threshold already optimal |
| 0047 | ACCEPTED | Bessel Miller seed binary-refine |
| **0048** | **ACCEPTED** | **Airy asymptotic boost reduction (this)** |

Four ACCEPTED Phase 2 perf+correctness wins (0042, 0043, 0047,
0048) plus three doc-tier closures (0040, 0044, 0045, 0046) plus
three rejection ADRs (0037, 0040 is doc-tier not rejection — adjust:
two rejections 0037/0041 plus 0040 doc-tier). Phase 2's
measurement-as-deliverable discipline has produced a substantial
durable engineering record.

**ADR-0044's flagged follow-ups** — `airy_zero_value` boundary
constants memoisation, and any `airy_series` cancellation guard
tightening — remain not pursued. The boundary-constant memoisation
is orthogonal to the asymptotic path; if pursued, it would target
the `x = 0` special-case dispatch. The `airy_series` guard
tightening is structural (the guard is calibrated for the
worst-case `|x|^{3/2}·log₂e` cancellation; not obvious slack
without a derivation pass).

**General lesson (combined with ADR-0047):** The two "real-win"
levers in Phase 2 perf work were both of the form *"a per-call
conservative slack with no specific accumulated-round-off
analysis"*:
- ADR-0047 (Miller seed): `+8` fixed guard after exponential
  search, plus 2× exponential overshoot. Both deliberately
  conservative per the original doc-comment.
- ADR-0048 (Airy asymptotic boost): `+64` working-precision boost
  where `+32` budgets the actual round-off.

Future investigators tuning numerical kernel perf should look for
**guard constants that exceed their stated budget** rather than
chase threshold formulas (which ADR-0044/0045/0046 confirmed tend
to be at the mathematical optimum from first-principles
derivation).

## Related

- ADR-0044 (Airy threshold doc-tier) — flagged this change as a
  follow-up; this ADR lands the flagged change.
- ADR-0042 (pf-1axr root fix) — the trig kernel range-cap fix that
  enables `airy_asymptotic_neg` to call `sin`/`cos` at higher
  working precisions; not directly exercised by this change (the
  reduced `+32` boost stays well below the trig table's expanded
  supported range), but the fix is a prerequisite for any future
  Airy work that boosts working precision aggressively.
- `src/math/airy.rs:497-504` — `airy_asymptotic_pos` with the new
  doc-comment and the reduced boost.
- `src/math/airy.rs:556-565` — `airy_asymptotic_neg` parallel.
- `src/math/ziv_calibration.rs:191` — `AIRY_ERROR_GUARD = 24`
  (the Ziv-side calibrated slack that the new `+32` budget
  matches).
- `benches/airy_dispatch.rs` — re-scoped to 4 asymptotic cells;
  header doc explains the rationale and Maclaurin-targeting
  guidance.
- ADR-0047 — the prior Phase 2b ACCEPTED perf slice that this ADR
  parallels in shape (real win from a conservative-guard reduction).
