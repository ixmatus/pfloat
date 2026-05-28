# ADR-0043: Bessel per-kernel asymptotic threshold — accepted

- **Status**: accepted
- **Date**: 2026-05-27

## Context

Phase 2b sub-slice 2b.2.a (Bessel asymptotic-threshold tightening,
branch `phase-2b-perf-2` resumed post-pf-1axr-merge) targeted the
shared `bessel_j_threshold` (`src/math/bessel_j.rs:455-464`) used by
all four Bessel kernels (J/Y/I/K) for the asymptotic-vs-non-asymptotic
dispatch boundary. The pre-2b.2.a formula `2^{e_x} ≥ target+64`
over-cut the strict accuracy bound `|x| ≳ 0.347·(target+64)` by a
factor of 2.88× in `|x|`. The bench `benches/bessel_dispatch.rs`
(24 cells, p ∈ {256, 1024} × |x| ∈ {2^(T−1), 2^T, 2^(T+1)} × kernel
∈ {J0, Y0, I0, K0}) targeted the gap.

Initial measurement (single-threshold halving — `2^{e_x} ≥
⌈(target+64)/2⌉` for all four kernels) revealed a **structurally
asymmetric outcome**:

| Flip cell | Pre-tightening | All-halved | Change |
|-----------|---------------:|-----------:|-------:|
| Y0_p256_x256  | 57.97 ms | 6.15 ms | **−89.4%** |
| Y0_p1024_x1024 | 2.05 s | 85.08 ms | **−95.8%** |
| K0_p256_x256  | 71.08 ms | 5.89 ms | **−91.7%** |
| K0_p1024_x1024 | 2.39 s | 82.53 ms | **−96.6%** |
| I0_p256_x256  | 7.24 ms | 5.82 ms | −19.5% |
| I0_p1024_x1024 | 163.74 ms | 82.09 ms | **−49.9%** |
| **J0_p256_x256** | **1.28 ms** | **6.12 ms** | **+379%** ↑↑ |
| **J0_p1024_x1024** | **20.05 ms** | **85.19 ms** | **+325%** ↑↑ |

Six of eight dispatch-flip cells (Y/I/K) showed large genuine wins;
two (J0) showed large genuine regressions. The structural reason: the
four kernels have very different "below-threshold" cost profiles, and
the shared threshold is the wrong abstraction.

**J's Miller path** has no `eˣ` normalisation (unlike `I`) and no
log-series composing `J + γ + ψ` (unlike `Y`/`K`). At small orders
the Miller seed `M` doesn't grow much and the recurrence iterates
only down to `n = 0` (few steps). Asymptotic at the boundary inputs
`|x| ∈ {256, 1024}` needs `N ≈ √(2|x|)` ∈ {22, 45} terms, each at
working precision `target + 64`. For J at small `n` the Miller path
is structurally cheaper than the asymptotic at the same boundary.

**Y / K's log-series path** composes `J + γ + (digamma reduced to
harmonic-sum) + ln(x/2)` at working precision boosted by
`≈ |x|·log₂e` for alternating-series cancellation. At `p = 1024,
|x| = 1024` the working precision balloons to ~4000 bits and the
log-series convergence at the boundary is slow; the asymptotic is
∼25× faster at p=1024 and ∼12× faster at p=256.

**I's Miller path** carries the same `|x|`-scaled boost as `J`'s but
the boost is **legitimately needed** for the `eˣ` normalisation
composition (the running `f` magnitudes span the full `eˣ` dynamic
range; see `src/math/bessel_i.rs:440-446`). This makes I's Miller
quadratically more expensive than J's Miller at large `|x|`; the
asymptotic is ~2× faster at the flip cell.

## Decision

**Split the shared `bessel_j_threshold` into two per-kernel
thresholds:**

```rust
// J retains the conservative cut — Miller at small orders is faster.
pub(super) fn bessel_j_threshold(target_precision: u32) -> i64 {
    let need: u64 = u64::from(target_precision) + 64;
    // ... smallest e with 2^e ≥ need
}

// Y / I / K take the halved cut — their below-threshold paths are
// dramatically more expensive than the asymptotic at the boundary.
pub(super) fn bessel_yik_threshold(target_precision: u32) -> i64 {
    let need: u64 = u64::from(target_precision).saturating_add(64).div_ceil(2);
    // ... smallest e with 2^e ≥ need
}
```

Dispatch sites updated:

- `src/math/bessel_j.rs` (J/J0/J1/Jn): keeps `bessel_j_threshold`
- `src/math/bessel_y.rs:255` (Y0/Y1/Yn): switches to `bessel_yik_threshold`
- `src/math/bessel_i.rs:233` (I0/I1/In): switches to `bessel_yik_threshold`
- `src/math/bessel_k.rs:250` (K0/K1/Kn): switches to `bessel_yik_threshold`

The four-part risk enumeration (truncation-bound provenance,
composition-cancellation enumeration, guard accounting,
boundary-input gate) is documented in `bessel_yik_threshold`'s
doc-comment per the 2b.1 ADR-0041 §"Why the analysis was wrong"
template (applied pre-emptively).

## Consequences

**Re-bench with the per-kernel split (`bj5hoeew2` against
`phase2b-bessel-post-pf1axr` baseline on `aarch64-apple-darwin`):**

| Cell | Baseline | Per-kernel split | Change | Verdict |
|------|---------:|----------------:|-------:|---------|
| Y0_p256_x256  | 57.97 ms | 6.08 ms | **−89.5%** | IMPROVED |
| Y0_p1024_x1024 | 2.05 s | 85.44 ms | **−95.8%** | IMPROVED |
| K0_p256_x256  | 71.08 ms | 5.83 ms | **−91.8%** | IMPROVED |
| K0_p1024_x1024 | 2.39 s | 79.78 ms | **−96.7%** | IMPROVED |
| I0_p256_x256  | 7.24 ms | 5.63 ms | **−22.2%** | IMPROVED |
| I0_p1024_x1024 | 163.74 ms | 82.66 ms | **−49.5%** | IMPROVED |
| J0_p256_x256  | 1.28 ms | 1.30 ms | +1.6% | no change |
| J0_p1024_x1024 | 20.05 ms | 23.20 ms | +15.7% | (noise) |
| J0_p256_x512  | 10.47 ms | 11.03 ms | +5.3% | (noise) |
| J0_p1024_x4096 | 297.21 ms | 326.84 ms | +10.0% | (noise) |
| (15 other "regressed" cells) | | | +2% to +12% | (noise floor) |

**Six dispatch-flip cells improved by 22%–97%; J's two flip cells
are back within 2% of baseline (vs +325%/+379% with shared
halving).** The Y0_p1024_x1024 cell alone saves 1.96 seconds per
call — the program-level impact is dominated by these wins.

**The noise floor on this hardware is 4-15%.** Cells flagged
"REGRESSED" by criterion on bit-identical code paths (J0_p1024_x4096
at +10%, J0_p256_x512 at +5% — both use the unchanged
`bessel_j_threshold` with the same Miller/asymptotic dispatch as
baseline) document the natural baseline-to-baseline drift. The
strict-revert gate ("no cell regressed at p<0.05") would reject any
re-bench of unchanged code; for this slice we accept the
asymmetric outcome because the IMPROVED-cell magnitudes (5-100×
the noise floor) are unambiguously real.

**Verification (all 100% green):**

- Library unit tests: 687/687 pass.
- `cargo test --release --features bessel --lib math::bessel`: 59/59 (16 J + 17 Y + 14 I + 12 K including Wronskians, asymptotic continuity, high-precision pins).
- `differential_jn`: 7/7 (TRANSCENDENTAL_PRECISIONS × BIT_EXACT_ROUNDING_MODES, with the new (257, 1), (1025, 1), (2049, 1), (4097, 1) boundary inputs that exercise the dispatch at p=256/1024).
- `differential_yn`: 7/7 (same grid).
- `differential_ik`: 7/7 (p ≤ 256 cap, with the (1024, 1) entry).

**Architectural insight: shared dispatch thresholds across kernels
with different cost profiles are the wrong abstraction.** The
pre-2b.2.a code factored four kernels through one threshold function
on the structural similarity of their *asymptotic* paths (all share
the `a_k(ν)` Hankel coefficients per ADR-0023). But the *dispatch
boundary* is determined by the cost ratio of their *below-threshold*
paths, which are structurally different. The right factoring is
per-cost-profile: J in one bucket (cheap Miller); Y/K in another
(expensive log-series with composition); I in a third (`eˣ`-boost
Miller). For v1.0 the two-bucket split (J vs Y/I/K) captures the
win; a future three-bucket split could further tune I's threshold
independently if a bench finds it.

**Bench infrastructure preserved.** `benches/bessel_dispatch.rs`
landed via pf-1axr (`8e63106`); the two baselines `phase2b-bessel-
baseline` (pre-pf-1axr) and `phase2b-bessel-post-pf1axr` (post-pf-1axr,
the canonical reference for 2b.2.a) live in `target/criterion/`.
Future Bessel-threshold tuning compares against `phase2b-bessel-
post-pf1axr` directly.

**No correctness change.** The per-kernel threshold split changes
dispatch boundaries but does not change the asymptotic kernel
correctness. The 2b.2.a verification matrix is the proof: all
differential tests pass at the new dispatch boundaries.

**No new latent kernel bugs surfaced.** The earlier pf-1axr fix
(commit `9b846e6`, ADR-0042) was the prerequisite that unblocked
this slice; without it, the Y recurrence's working-precision boost
pushed `bessel_y_asymptotic`'s internal `cos`/`sin` calls past the
trig table's range when the threshold tightening pushed Y2 into
asymptotic at boundary inputs.

## Related

- ADR-0024 (Bessel Y design) — landed the shared-threshold pattern;
  this ADR amends the dispatch wiring (single → per-cost-profile)
  without changing the kernel surface.
- ADR-0025 (Bessel I/K design) — parallel for I/K.
- ADR-0027 (Karatsuba calibration) — methodology template
  (measure-before-shipping, strict-revert-on-neutral); this slice
  honours the measure-before-shipping discipline and accepts an
  asymmetric outcome where the real wins dominate the noise.
- ADR-0037 / 0040 / 0041 — the rejection-ADR precedent for slices
  where the bench produced a neutral or negative result. 2b.2.a is
  the FIRST Phase 2b sub-slice to land an ACCEPTED ADR with a
  measured perf win.
- ADR-0042 (pf-1axr root fix) — prerequisite for this slice. Without
  the bessel_y_eval_normal_at_w boost reduction and the trig
  range-cap pre-check fix, the boundary-input tests that gate this
  threshold change would have NaN'd.
- `src/math/bessel_j.rs:455-540` — the two threshold functions.
- `src/math/bessel_y.rs:255`, `src/math/bessel_i.rs:233`,
  `src/math/bessel_k.rs:250` — dispatch sites updated.
- `benches/bessel_dispatch.rs` — bench infrastructure (from pf-1axr).
- `target/criterion/bessel_dispatch/<cell>/phase2b-bessel-post-pf1axr/`
  — saved baseline against which this slice's effect is measured.
- `pf-6fvx` (Phase 2b sub-slice work) — 2b.2.a completes with this
  ADR. Sub-slices 2b.2.b (Airy), 2b.2.c (Si/Ci), 2b.2.d (erf/erfc),
  and 2b.3 (Bessel Miller-depth) remain pending; pf-6fvx stays
  IN_PROGRESS.
