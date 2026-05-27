# ADR-0042: pf-1axr trig range-cap pre-check and bessel_y recurrence boost — root fix

- **Status**: accepted
- **Date**: 2026-05-27

## Context

Sub-slice 2b.2.a (Bessel asymptotic-threshold tightening, branch
`phase-2b-perf-2`) added boundary-input cases `(1025, 1)` and
`(2049, 1)` to `tests/differential_yn.rs`'s `DYADIC` table.
`yn_negative_order_and_dyadic_matches_mpfr` then failed:

```
Y2(1025/1) at p=53
  left:  NaN
  right: -1.4253945762591734e-3 (MPFR)
```

The failure was traced (`tests/bessel_y::tests::probe_pf1axr_*`,
`#[ignore]`'d in `src/math/bessel_y.rs`) through three layers:

1. `bessel_y_kernel(n=2)` → `ziv_round` → `bessel_y_eval_normal_at_w`.
2. `bessel_y_eval_normal_at_w` boosted working precision to 3061 via
   `extra = mag·23/16` with `mag = 2^(e_x+1) = 2048` (`e_x = 10` for
   `|x| = 1025`), yielding `extra = 2944`.
3. The recurrence's `bessel_y01(0, x, 3061, true)` →
   `bessel_y_asymptotic(0, x, 3061)` internally bumped working to
   3125 and called `omega.cos(NearestEven)` / `omega.sin(NearestEven)`.
4. `sin_kernel` / `cos_kernel`'s range-cap pre-check at
   `ziv_max_working = target_precision + 1024 = 4085` invoked
   `reduce(omega, 4085)`. With `e_omega = 10`, the reduce condition
   `e_x + working + 64 < 4096` evaluates to `10 + 4085 + 64 = 4159 ≥
   4096`, so reduce returned `None` and the kernel returned
   `NaN + INVALID`.

The bug was independent of the threshold tightening (confirmed by
stashing the change and re-running). It had been latent since the
slice that landed `bessel_y` (slice 6p, ADR-0024): the
`next_i64_in(state, 1, 40)` random sweep at `p = 53` in the prior
test grid never reached `|x| ≥ 128` (the asymptotic threshold at
`p = 53`), so `Y_n` asymptotic at small target precision was
untested for any `|x|`. The new boundary inputs were the first to
hit the regime.

A direct trig sweep confirms the structural defect:

| target_p | 2800 | 2900 | 2950 | **3000** | 3050 | 3072 | 3100 | 3125 |
|----------|------|------|------|----------|------|------|------|------|
| `sin(1025)` / `cos(1025)` | ✓ | ✓ | ✓ | NaN | NaN | NaN | NaN | NaN |

The cliff at `target_precision ≥ 3008 − e_x` (10 in this case)
arises because the 4096-bit `2/π` reduction table satisfies
`reduce` iff `e_x + working_prec + 64 < 4096`; the pre-check at
`target + 1024` (the Ziv worst-case ceiling) fails when
`target + 1024 + e_x + 64 ≥ 4096`. Any caller asking for
`sin`/`cos` at target precision ≥ 3000 with `|x| > 1` hit it.

The blast radius across the six asymptotic kernels using
`sin`/`cos` is narrower than the cliff suggests: `bessel_j_asymptotic`
(`+64`/`+512` boost), `si_asymptotic` (`+64`/`+512`),
`ci_asymptotic` (`+64`/`+512`), and `airy_asymptotic_neg` (`+64`)
all use modest boosts that stay below the cliff at any
input. Only `bessel_y_eval_normal_at_w`'s `|x|`-scaled
recurrence boost (`extra ≈ 2.875·|x|`, capped at 4096 bits) reached
the cliff, and only via the recurrence path's internal call to
`bessel_y_asymptotic` at the boosted precision.

## Decision

Two fixes, applied together:

**Fix 1 — `sin.rs` / `cos.rs` / `tan.rs` range-cap pre-check.**
Change the pre-check from `ziv_max_working = target + 1024` (the Ziv
worst-case ceiling) to `ziv_first_working = target +
ZIV_BASE_GUARD` (the Ziv first-iteration working precision, where
`ZIV_BASE_GUARD = 64` is now `pub(super)` in `ziv.rs`). If reduce
fails at the first iteration the input is fundamentally out of
range and no Ziv iteration could recover; if it succeeds the first
iteration runs but higher Ziv iterations may fail. To handle the
latter without panicking, change each closure's
`reduce(x, w).expect("range-cap pre-checked")` to graceful match:
return `NaN` at the working precision when reduce returns `None`.
The Ziv driver propagates the `NaN`; a post-Ziv check raises
`INVALID` if the final result is `NaN` and the status didn't
already capture it.

The pre-check still catches genuinely-out-of-range inputs (those
where even the first iteration can't reduce). It no longer fires
spuriously for inputs in the range `target + 64 ≤ working_max`
where Ziv would only lift to higher working if the interval test
fails. The cliff for genuine failure moves from
`target ≥ 3008 − e_x` to `target ≥ 4032 − 64 − e_x = 3968 − e_x`,
gaining ~960 bits of supported precision.

**Fix 2 — `bessel_y_eval_normal_at_w` recurrence boost.**
Change `extra` from `(mag·23/16).min(4096)` (the alternating-
series-cancellation budget, matching `bessel_y_series` correctly)
to `32 + 4·m` (matching the recurrence's actual per-step error
amplification bound). The recurrence
`Y_{k+1} = (2k/x)·Y_k − Y_{k−1}` cancels only when `(2k/x)·Y_k`
and `Y_{k−1}` have the same sign and similar magnitude; for
`x ≫ k` the `(2k/x)` factor is small and `|Y_{k+1}| ≈ |Y_{k−1}|`
with no cancellation. For `x ≲ k` the amplification per step is
bounded by `1 + x/(2k)` (≤ 1 bit per step). For typical orders
`m ≤ 20`, `32 + 4·m ≤ 112` bits is comfortably above the worst
case; the `|x|`-scaled budget over-provisioned by factors of 30+
in the asymptotic regime.

## Consequences

**`Y_n` at large `|x|` and small `p` is now correct.** Verified:
- `Y2(1025)` at `p=53`: returns `-1.4253945762591734e-3` (was `NaN`)
- `Y2(2049)` at `p=53`: returns `1.8104618091279389e-3` (was `NaN`)
- `Y2(4097)` at `p=53`: returns correct value (was `NaN`)
- All differential tests pass under both fixes:

| Test | Before | After |
|------|--------|-------|
| `differential_jn` | 7/7 | 7/7 |
| `differential_yn` | 6/7 (Y2 NaN) | 7/7 |
| `differential_ik` | 7/7 | 7/7 |
| `differential_sin`/`cos`/`tan` | 4/4 each | 4/4 each |
| `differential_si`/`ci` | 5/5 each | 5/5 each |
| `differential_ai` / `bi` | 5/5 + 6/6 | 5/5 + 6/6 |
| `differential_erf`/`erfc` | 4/4 each | 4/4 each |
| Library unit tests | 687/687 | 687/687 |

**The trig kernel's supported precision range expands.** Callers can
now compute `sin`/`cos`/`tan` at target precisions up to
`3968 − e_x` (was `3008 − e_x`); for `|x| ≤ 1` the supported range
extends to `~3968 bits`. Beyond this, the kernels still return
`NaN + INVALID` cleanly. The 4096-bit `2/π` reduction table remains
the fundamental limit; expanding it (to e.g. 8192 bits) is a
separate v1.x perf/audience-expansion item, not blocking for v1.0.

**`Y_n` recurrence is faster.** For large `|x|`, the working
precision drops from `target + 64 + min(2.875·|x|, 4096)` to
`target + 64 + 32 + 4·m ≤ target + 176`. At `|x| = 2049, m = 3,
target = 1024`: working drops from 4149 → 1200, a ~3.5× reduction
in working precision for the Y0/Y1 base computation. Bench
infrastructure (`benches/bessel_dispatch.rs`) is in tree from
sub-slice 2b.2.a; the recurrence speedup is a side benefit not
measured here (the 2b.2.a measurement is the threshold-tightening
bench, separate from this correctness fix).

**Boundary-input audit confirmed narrow blast radius.** Of the six
asymptotic kernels using `sin`/`cos`/`tan`, only `bessel_y` had
the over-boost pattern. `bessel_j_asymptotic`, `si_asymptotic`,
`ci_asymptotic`, and `airy_asymptotic_neg` use modest `+64`
(capped `+512`) boosts that stay well below the table's supported
range at any input. `bessel_i_miller` and `bessel_k_climb` use the
same `|x|`-scaled boost but legitimately need it (the `eˣ`
normalization composition in I; could be tightened in K's case
since `K` uses `exp`/`cosh` not `sin`/`cos`, but doesn't cause
correctness bugs and is out of pf-1axr scope).

**Latent-bug risk for future high-precision callers.** A
caller asking for `sin`/`cos`/`tan` at `target ≥ 3008 − e_x` would
now succeed where before it would have spuriously failed. This
includes any caller that boosts working precision for cancellation
compensation. For v1.0 the existing test surface confirms this is
safe; the post-v1.0 work item is to expand the reduction table
(deferred to v1.x).

**Test coverage gap closed.** The latent bug had been undetected
since slice 6p (Y kernel landing) because no test exercised
`Y_n` asymptotic at small `p` and large `|x|`. The new boundary
inputs `(257, 1)`, `(1025, 1)`, `(2049, 1)`, `(4097, 1)` in
`differential_jn.rs` and `differential_yn.rs`, plus `(1024, 1)`
in `differential_ik.rs`, close this gap. The audit-equivalent for
`Si`/`Ci` shows their tables already cover `x ≤ 5000` at
`TRANSCENDENTAL_PRECISIONS`; for Airy the existing
`(±150, ±180)` at p=53 covers asymptotic and the Wronskian at
p=200 binds Ai/Bi/Ai′/Bi′.

**Diagnostic test artifacts.** `src/math/bessel_y.rs::tests`
carries three `#[ignore]`'d probe tests (`probe_pf1axr_y2_1025_p53`,
`probe_pf1axr_asymptotic_internals`,
`probe_pf1axr_trig_range_cap_sweep`) that document the diagnosis
trace. They don't run by default but are immediately invokable for
future regression investigation.

**Bench infrastructure carried.** `benches/bessel_dispatch.rs`
(24 cells, 2 precisions × 3 `|x|` straddles × 4 kernels) and its
`[[bench]]` entry in `Cargo.toml` land here as part of the sub-slice
2b.2.a setup that is now unblocked. The `phase2b-bessel-baseline`
criterion baseline on disk is preserved; the resumed 2b.2.a perf
work compares against it.

## Related

- `pf-1axr` (the bead this slice closes), discovered-from `pf-6fvx`.
- `pf-6fvx` (Phase 2b sub-slice work) — was blocked by pf-1axr;
  unblocks after this lands. Sub-slice 2b.2.a resumes from the
  threshold-tightening step with the bench infrastructure already
  in tree.
- ADR-0024 (Bessel Y design) — the slice that landed the over-
  boost pattern in `bessel_y_eval_normal_at_w`. ADR-0042 amends
  the recurrence boost calibration without changing the surface.
- ADR-0022 (Ziv envelope, slice p1.4) and ADR-0038 (Phase 1f
  Ziv-driver expansion to directed modes) — the architectural
  context for the trig kernel's pre-check.
- ADR-0039 (Phase 1g per-kernel error_guard calibration) — the
  precedent for `SIN_ERROR_GUARD`/`COS_ERROR_GUARD`/`TAN_ERROR_GUARD`
  named constants the kernels supply to `ziv_round`.
- `src/math/sin.rs:117-145` — sin_kernel range-cap pre-check (fixed)
- `src/math/cos.rs:112-142` — cos_kernel parallel (fixed)
- `src/math/tan.rs:111-148` — tan_kernel parallel (fixed)
- `src/math/bessel_y.rs:306-321` — `bessel_y_eval_normal_at_w` boost
  (fixed)
- `src/math/ziv.rs:35-42` — `ZIV_BASE_GUARD` now `pub(super)`
- `tests/differential_jn.rs:36-49`, `tests/differential_yn.rs:44-55`
  — boundary inputs
- `tests/differential_ik.rs:78-90` — `IK_TABLE` (1024, 1) entry
- `feedback_precision_gated_verification_surface` (2b.1 lesson)
  — argument-gated dispatch analog applies: a boundary-input
  addition surfaced a latent bug independent of the threshold
  change. The corpus-extension defect-rollback pattern was
  considered but rejected in favor of root-fix-and-audit per user
  direction.

## Bench infrastructure note

`benches/spouge_lgamma.rs` (sub-slice 2b.1, ADR-0041) is the
template; `benches/bessel_dispatch.rs` mirrors its criterion-group
shape (`measurement_time = 20s`, `sample_size = 20`,
`warm_up_time = 2s`). The 24-cell baseline at
`phase2b-bessel-baseline` is on disk; the future 2b.2.a perf
slice diffs the proposed-tightening run against it via
`--baseline phase2b-bessel-baseline`. The bench file isn't
strictly part of the pf-1axr correctness fix but lands together
because the sub-slice 2b.2.a setup that produced it is paused on
this commit and resumes after merge.
