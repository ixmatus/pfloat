# ADR-0126: the fixed-cap series boost gives way to cancellation_boosted (the pf-6naq crate sweep)

- **Status**: accepted
- **Date**: 2026-07-05

## Context

The R4.12 completeness probe filed pf-6naq alongside pf-1vzg (ADR-0125).
Where pf-1vzg concerned the DIVERGENT asymptotic paths, pf-6naq concerns
their CONVERGENT series siblings. Across the special-function surface a
single idiom sized the working-precision boost for a cancelling series:

```rust
let extra = if e_x <= 0 { 64 } else {
    let shift = (e_x + 1).min(20) as u32;
    let mag: u64 = 1u64 << shift;
    (mag.saturating_mul(23) / 16).min(4096) as u32   // <- FIXED cap
};
let working = target.saturating_add(64).saturating_add(extra)
    .min(target.saturating_add(4096));
```

The estimate `mag·23/16 ≈ |x|·log₂ e` is the alternating cancellation an
argument of magnitude `|x|` produces (the peak partial term is `≈ 2^{|x|·
log₂ e}` and cancels back to the `O(1)` result). Two things break it:

1. **The `.min(4096)` cap.** When the realised cancellation `C` exceeds
   `≈ 4096 + guard`, the boost undershoots. The series is then evaluated
   with fewer than `target` accurate bits, and — the decisive point — the
   Ziv half-width model `|y|·2^-(working − guard)` becomes UNSOUND: the
   true error exceeds the claimed interval, so `ziv_round` certifies a
   WRONG value. This is the pf-1vzg class one tier shallower: the fix is
   to charge the REALISED cancellation, not a fixed constant. Reachable
   through the ordinary public API at `target > ~4096` with `|x|` near the
   Maclaurin/asymptotic regime boundary — no near-zero input required.

2. **The `shift.min(20)` and the estimate itself** are moot once the boost
   is realised rather than estimated; `cancellation_boosted` (ADR-0110)
   already measures the peak partial term and grows working by exactly the
   cancellation depth. R5.1 fixed the first three sites (`ci_series` and
   the Ci/J/Y convergent fallbacks) as a side effect of ADR-0125; pf-6naq
   is the crate-wide sweep of the remaining sites.

Reproducers (`tests/regression_review_2026_07_05_r52.rs`, mpmath 1.4.1 on
bit-identical dyadics, release-gated deep rows):

| kernel | deep input | realised C | pre-fix |
|---|---|---|---|
| `erf` (Maclaurin) | `erf(64)`, target 6000 | `≈ 5909` | wrong |
| `erfc` (`1 − erf`) | `erfc(64)`, target 6000 | `≈ 2·5909` | wrong |
| `Ei` (`x < 0` series) | `Ei(−1900)`, target 2900 | `≈ 5482` | wrong |
| `Si` | `Si(4000)`, target 5000 | `≈ 5771` | wrong |
| `K₀` (DLMF 10.31.1 log series) | `K₀(4000)`, target 6000 | `≈ 11542` | wrong |
| `Y₀` (DLMF 10.8.1 log series) | `Y₀(4000)`, target 6000 | `≈ 5771` | wrong |
| `J₀` (Miller, middle regime) | near `j_{0,40}`, `D ≈ 4998` | `≈ D` | wrong |

Two findings sharpened the picture:

- **erfc double-cancels.** `erfc(x) = 1 − erf(x)` stacks the erf Maclaurin
  peak (`≈ 2^C`) AND the `1 − erf` subtraction (`erfc ≈ 2^{−C}` at the
  boundary), so `op_scale − result_exp ≈ 2C`. Its old inner
  `.min(w + 512)` cap bit at HALF erf's argument: `erfc(40)` was already
  wrong pre-fix where `erf(40)` was correct. The shallow control moved to
  `erfc(20)`.

- **`bessel_j_miller` is pf-1vzg-class, not pf-6naq-class.** A GENERIC
  middle-regime `J` barely cancels: the huge sum-rule normalisation
  `c = 1/J_M` divides out in `J_m = f_m / S`, leaving only `≈ ½·log₂|x|`
  bits. The cap bites ONLY for a deep near-zero `J` (`f_m → 0`), where the
  cancellation is the proximity depth `D`. Because the middle-regime
  estimate `≈ |x|·log₂ e` is small there (`≈ 184` at `|x| ≈ 125`), the
  effective cap is low and `D` past `≈ 1300` already certifies wrong —
  moderate depth, the pf-1vzg near-zero family one regime in from the
  asymptotic R5.1 handled. Unlike the asymptotic, Miller has NO truncation
  floor, so `cancellation_boosted` always resolves it (no reliability
  test needed).

## Decision

Replace the fixed-cap boost with the realised-cancellation boost at every
genuinely-cancelling site, uniform with R5.1:

1. **Each series returns `(value, op_scale)`,** where `op_scale` is the
   largest partial term's exponent lifted to the result scale, and its
   caller wraps the evaluation in `super::ziv::cancellation_boosted`. The
   internal working precision drops to a small fixed `target + 64` guard;
   the boost is driven by the returned `op_scale`, so no fixed cap sits
   between the series and its correct rounding.
   - `erf_maclaurin` (shared by `erf` and `erfc`; `erfc`'s op_scale is the
     erf peak term, its own `1 − erf` subtraction charged by the tiny
     result exponent).
   - `ei_series` — now tracks the peak term (not just the leading
     `γ + ln|x|`), and the `in_zero_window` special branch is REMOVED: the
     general path's `cancellation_boosted` subsumes both Ei's near-zero
     (ADR-0110) and the large-`|x|` alternating cancellation.
   - `si_series`, `bessel_y_series` (its plain call in `bessel_y01`'s
     moderate-`x` branch now boosted, matching the R5.1 fallback path).
   - `bessel_j_miller` returns `op_scale = f_max_exp − S_exp` (the peak
     recurrence value on the `J_m = f_m / S` scale).

2. **`bessel_k_series` self-contains the boost.** It has multiple callers
   (`bessel_k01` directly, `bessel_k_climb` for its seeds), so wrapping at
   each site would duplicate; instead `bessel_k_series` wraps its own body
   (`bessel_k_series_at` returning `(value, op_scale)`) so every caller
   gets a resolved value.

3. **The dead `ci.rs` cap is removed** (superseded by R5.1's op_scale
   return).

4. **Audited and deliberately unchanged:**
   - `bessel_i_miller` and `bessel_k_climb` — all-positive recurrences
     (`(2k/x)·f_k ± f_{k∓1}`) with NO subtractive cancellation; their
     `|x|·log₂ e` boost is conservative dynamic-range provisioning, and the
     true need is `≈ log₂ M`, far below the cap. `I₀(3000)` at target 6000
     is bit-exact pre-fix. The `bessel_k_climb` seeds now come from the
     fixed self-boosting `bessel_k_series`, so `K_{n≥2}` is correct
     transitively.
   - `beta.rs` `exact_sum` — the `+4096` is DoS-budget slack in an exact-
     addition span check, not a cancellation cap.
   - `zeta.rs` `zeta_fe` — the boost is the input-proportional
     `pole_proximity_depth`; the `+4096` sits only in a non-binding
     backstop cap that scales with `x.precision()`.

## Consequences

- The fixed-cap-undershoot family is closed for `erf`, `erfc`, `Ei`, `li`
  (via `Ei`), `Si`, `K₀`/`K₁`/`K_n`, `Y₀`/`Y₁`, and `J`'s middle regime.
  Each is bit-exact against mpmath at a deep target where the pre-fix
  kernel was wrong, with a shallow control that was already correct and
  stays so.
- **Frugality, not a tax.** Removing the fixed `+extra` (up to `+4096` on
  EVERY call) makes the common small-argument path CHEAPER — the boost now
  provisions only what the argument needs, two small evaluations rather
  than one over-wide one. The deep near-boundary rows cost more (input-
  proportional, DoS-budget posture, release-gated), which is the price of
  correctness there.
- **Inversion: the reproducer must exceed the CURRENT cap, and the cap is
  argument-dependent.** For `erf`/`Si`/`K`/`Y` that is `target > ~4096`
  with a large regime-boundary `|x|`. For `bessel_j_miller` the cap is the
  SMALL middle-regime estimate (`≈ 184`), so a moderate `D ≈ 1300`
  near-zero already exceeds it — labelling it a `> 4096` edge would have
  mis-scoped it as rare. Grep the whole crate for the shape (it was in
  eight kernels); do not trust the `4096` literal as the effective bound.
- **Inversion: not every `min(4096)` is a cancellation cap.** The bead
  named `bessel_i` and `bessel_k_climb`; both are all-positive recurrences
  where the boost is dynamic range, not cancellation, and
  `cancellation_boosted` (which measures `op_scale − result_exp ≈ 0`
  there) would add nothing. The verdict was verified in both directions:
  `I₀` is bit-exact pre-fix, so it was left alone. Pattern-matching the
  idiom without the cancellation analysis would have churned correct code.
- **Inversion: `bessel_j_miller`'s generic case must NOT be over-boosted.**
  The op_scale `f_max_exp − S_exp` tracks `result_exp` for a generic `J`
  (the normalisation divides out), so `cancellation_boosted` charges
  `≈ 0` and the boost fires only near a zero. Using `f_max_exp` alone
  (without the `− S_exp` normalisation) would have inflated every
  middle-regime `J` by the seed magnitude.
- The differential lanes (`differential_erf`, `_ei`, `_si`, `_kn`, `_yn`,
  `_jn`) and the property suites guard the common path; the R5.1
  regression rows (Ci/Airy/J₀/Y₀ near-zero) still pass, confirming the
  `ci_series`/`bessel_y_series` cap removals did not perturb them.

## References

- pf-6naq (the bug), pf-1vzg / ADR-0125 (the asymptotic sibling), pf-tzqz
  (the R4.12 probe that filed it), ADR-0110 (`cancellation_boosted`),
  ADR-0097 (the realised-cancellation / operand-scale posture).
- DLMF §6.6 (erf/Ei/Si series), §10.31.1 (K log series), §10.8.1 (Y log
  series), §10.6.1 / §10.12.4 (J Miller recurrence + sum rule).
- `tests/regression_review_2026_07_05_r52.rs`; oracle generators
  `scratchpad/gen_r52_oracles.py`, `scratchpad/gen_r52_besselj.py`.
