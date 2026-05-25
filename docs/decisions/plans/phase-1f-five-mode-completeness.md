# Phase 1f audit: five-mode kernel completeness

This is the in-tree audit document ratifying ADR-0038's strategic
commitment. It walks every kernel on the frozen v1.0 unary surface
plus the multi-arg surface, derives the `eval(w)` shape from the
kernel source (not recalled), identifies cancellation regimes that
would defeat the Ziv interval test, and records the per-regime
strategy each per-family slice will execute against.

The audit is the load-bearing artifact for Phase 1f. Per
ADR-0038's no-narrowing principle and CLAUDE.md's derive-don't-recall
discipline, each kernel entry derives its `eval(w)` shape from the
kernel source citation, its cancellation analysis from the eval
shape, and its Ziv strategy from the cancellation regime. The
8a case-4 O(m) DoS precedent and the tanh subnormal cancellation
precedent stay live: recalled cancellation soundness ships as
latent has-errors.

This document is built incrementally across slice p1.22's audit
phase. The structural skeleton, the per-kernel template, and the
already-Ziv-wrapped cohort (sanity entries) plus the elementary
completions family land in the first p1.22 session. Subsequent
sessions populate one family at a time per the task list.

## Per-kernel entry template

Each kernel block uses the following structure:

```
### <kernel_name>

- **Source**: `src/math/<file>.rs:<lines>`
- **Status today**: <ziv-wrapped | fixed-64-bit-guard | composition>
- **eval(w) shape**: <derivation from source; identify the
  intermediate computations that run at working precision w>
- **Cancellation regimes**: <each input range that produces
  pathological cancellation; for each, whether eval(w) can return
  exact zero at some w while the true value is non-zero>
- **Per-regime Ziv strategy**:
  - <regime A>: <drop-in wrap | short-circuit BEFORE Ziv |
    reformulation | basis change>
  - <regime B>: ...
- **Cited spec**: <DLMF section, IEEE 754-2019 §, paper reference>
- **Oracle coverage**: <MPFR-primary | Arb-primary | mpmath cross-check>
- **Estimated Ziv iterations at cap**: <number; based on the
  representative input set>
- **Worked example**: <known hard input + the kernel's behavior at
  it, derived from the source not the runtime>
- **Migration commit shape**: <number of commits the family slice
  ships for this kernel, per the one-concern-per-commit rule>
```

Entries marked **AUDIT TBD** are populated in subsequent slice
p1.22 sessions; the structural skeleton stands so the per-family
slices can claim entries as they open.

## Already-Ziv-wrapped cohort (sanity entries)

These kernels are correctly rounded under every IEEE 754-2019 mode
today via the existing `ziv_round` wrap (slices p1.2 and p1.4).
Phase 1f's work on this cohort is differential lane widening
(slice p1.23 scaffolding folds pf-lw3l) plus per-mode status TOML
sweep coverage; no kernel-side changes.

### pow

- **Source**: `src/math/pow.rs` (ADR-0022, slice 7c)
- **Status today**: Ziv-wrapped. Five-mode correct.
- **eval(w) shape**: integer exponent fast path via
  square-and-multiply (≤ 64 multiplies); general case via
  `exp(y · ln(x))` at working precision w. Each multiplication
  accumulates ≤ 2⁶ ULP of NE rounding error, comfortably below
  the `ZIV_ERROR_GUARD = 24` slack.
- **Cancellation regimes**: y · ln(x) near zero when x near 1 (the
  log1p analogue); the existing pow path boosts working precision
  by the lost-bit count before invoking the exp·ln composition. No
  collapse-to-exact-zero regime on the documented domain.
- **Per-regime Ziv strategy**: drop-in wrap (already shipped).
- **Cited spec**: IEEE 754-2019 §9.2 `pown`/`pow`; ADR-0022 for the
  Ziv interval test.
- **Oracle coverage**: MPFR-primary (with NearestAway synthesis at
  p+128 per `differential_pow.rs:45-58`).
- **Estimated Ziv iterations at cap**: 1-2 on most inputs; up to 5
  on measure-zero exact-tie inputs (the documented cap caveat).
- **Worked example**: `pow(63, -3)` at p=53 under NearestAway —
  `differential_pow.rs` exercises this; the integer fast path
  evaluates `63^(-3) = 1/250047` at working precision and the
  interval test certifies.
- **Migration commit shape**: no new commits (already shipped).

### exp

- **Source**: `src/math/exp.rs:1-66` (slice p1.2, ADR-0022).
- **Status today**: Ziv-wrapped. Five-mode correct.
- **eval(w) shape**: range reduce `k = round(x / ln(2))`,
  `r = x − k · ln(2)` with `|r| ≤ ln(2)/2`; Taylor series
  `exp(r) = 1 + r + r²/2! + …` summed at working precision w with
  ~4w iterations to reach below `2^-w` per-term; compose
  `exp(x) = exp(r) · 2^k` via free exponent shift.
- **Cancellation regimes**: none on the documented domain. Slice
  p1.2 closed the slice-8b underflow defect.
- **Per-regime Ziv strategy**: drop-in wrap (already shipped).
- **Cited spec**: standard exp algorithm; LN2_LIMBS_1024 (slice 7b2
  audit-corrected, project_ln2_constant_defect).
- **Oracle coverage**: MPFR-primary.
- **Estimated Ziv iterations at cap**: 1-2.
- **Worked example**: `exp(0xc0874385446d71c3)` ≈ `exp(-744.44)` →
  `f64::MIN_POSITIVE_SUBNORMAL` (the slice-8b corpus's underflow
  block; closed in p1.2 by the Ziv envelope).
- **Migration commit shape**: no new commits.

### ln

- **Source**: `src/math/ln.rs` (slice p1.2, ADR-0022).
- **Status today**: Ziv-wrapped. Five-mode correct.
- **eval(w) shape**: range reduce by extracting the binary exponent
  e: `x = m · 2^e` with `m ∈ [1, 2)`; `ln(x) = ln(m) + e · ln(2)`;
  `ln(m)` via `atanh` series on `(m−1)/(m+1)` ∈ `(-1/3, 1/3]`.
- **Cancellation regimes**: x near 1, where ln(x) ≈ x − 1; the
  `atanh` reformulation absorbs this. No collapse-to-exact-zero.
- **Per-regime Ziv strategy**: drop-in wrap (already shipped).
- **Cited spec**: standard ln algorithm via atanh series.
- **Oracle coverage**: MPFR-primary.
- **Estimated Ziv iterations at cap**: 1-2.
- **Worked example**: the log2 / log10 inheritance path — pfloat's
  `log2(x) = ln(x)/ln(2)` and `log10(x) = ln(x)/ln(10)`
  compositions inherit five-mode correctness because the divisor is
  exact and `ln(x)` is Ziv-driven.
- **Migration commit shape**: no new commits.

### tanh

- **Source**: `src/math/tanh.rs:1-67` (slice p1.2 Ziv envelope;
  slice p1.4 tiny-x short-circuit).
- **Status today**: Ziv-wrapped with tiny-x short-circuit. Five-mode
  correct.
- **eval(w) shape**: `tanh(|x|) = (1 − e^{−2|x|}) / (1 + e^{−2|x|})`
  on `|x|` to avoid the `−∞/+∞` indeterminate at `x = −∞`. For `|x|`
  small enough that the cubic Taylor correction falls below the Ziv
  driver's error guard, the kernel short-circuits to `|x|` with the
  input's sign.
- **Cancellation regimes**: tiny |x|, where `e^{−2|x|}` rounds to 1
  exactly at every working precision and the numerator collapses to
  0. This is the `feedback_ziv_interval_test_and_mpfr_rnda` §3
  archetype: `half_width(0) = 0` certifies 0 as the wrong answer.
  The short-circuit BEFORE the Ziv driver sees the cancellation
  path is the fix.
- **Per-regime Ziv strategy**: short-circuit at small |x|, drop-in
  wrap on the composition path (already shipped).
- **Cited spec**: standard tanh algorithm; slice p1.4 closes pf-7d7.
- **Oracle coverage**: MPFR-primary.
- **Estimated Ziv iterations at cap**: 1-2 on the composition path.
- **Worked example**: `tanh(2^-149)` (smallest positive f32
  subnormal) — short-circuit returns the input.
- **Migration commit shape**: no new commits.

### lgamma

- **Source**: `src/math/lgamma.rs` (slice p1.2, ADR-0022).
- **Status today**: Ziv-wrapped. Five-mode correct.
- **eval(w) shape**: Stirling asymptotic series + reflection
  `lgamma(x) = ln|π / (sin(π·x) · Γ(1−x))|` for x < 0.5. Working
  precision w drives both branches; the Stirling truncation is
  picked so the residual is below `2^-w`.
- **Cancellation regimes**: x near positive integer roots of
  `sin(π·x)` (the reflection branch) — handled by the existing pole
  detection in the kernel. No collapse-to-exact-zero in the
  composition.
- **Per-regime Ziv strategy**: drop-in wrap (already shipped). The
  reflection branch composes through `sin`; when `sin` migrates at
  p1.26 the composition inherits five-mode correctness already,
  because `sin` is called at a high working precision under
  NearestEven and the outer `ziv_round` envelope handles the final
  rounding mode.
- **Cited spec**: DLMF 5.11 (Stirling), 5.5.3 (reflection).
- **Oracle coverage**: MPFR-primary.
- **Estimated Ziv iterations at cap**: 1-3 depending on regime.
- **Worked example**: `lgamma(0.5) = ln(√π) ≈ 0.5723649...`.
- **Migration commit shape**: no new commits.

### erf

- **Source**: `src/math/erf.rs` (slice p1.4, ADR-0022 Ziv envelope).
- **Status today**: Ziv-wrapped. Five-mode correct.
- **eval(w) shape**: Maclaurin series `erf(x) = (2/√π) Σ
  (-1)^k x^(2k+1) / (k! (2k+1))` for small |x|; continued-fraction
  or composition with `erfc` for larger |x|. The Maclaurin path
  runs at `verification_precision = 53` (per slice p1.4) to handle
  the f32-subnormal-grid midpoint trap (pf-z0f). The 2/√π constant
  is from `TWO_OVER_SQRT_PI_LIMBS_1024` (slice 6m-audit corrected
  for the 1-ULP truncation defect).
- **Cancellation regimes**: large |x| where `erf(x) → ±1` and the
  remaining bits are in the deep tail; the composition through
  `erfc` is the cancellation-resistant path. No collapse-to-exact-zero.
- **Per-regime Ziv strategy**: drop-in wrap (already shipped).
- **Cited spec**: DLMF 7.6.1 (Maclaurin); IEEE 754-2019 §9.4 erf.
- **Oracle coverage**: MPFR-primary.
- **Estimated Ziv iterations at cap**: 1-3.
- **Worked example**: `erf(0x33800000)` (input 2^-24) — slice p1.4
  pf-z0f f32-subnormal-grid midpoint case.
- **Migration commit shape**: no new commits.

### bessel_j (J0, J1, Jn)

- **Source**: `src/math/bessel_j.rs` (slice p1.4 Ziv envelope).
- **Status today**: Ziv-wrapped at `verification_precision = 320`
  per slice p1.4 (BesselJ family's sub-midpoint cubic Maclaurin
  correction lives below p=24 ULP, the f32-grid harness needs the
  higher precision to retain the correction; per
  `feedback_bf_to_f32_directed_mode` the bumped path is NE-only in
  the f32 sweep). Five-mode correct on the kernel side.
- **eval(w) shape**: Maclaurin series for small |x|; Miller's
  backward recurrence for moderate orders; Hankel asymptotic for
  large |x|. Each branch runs at working precision w; the regime
  dispatch is in the kernel.
- **Cancellation regimes**: x near zeros of J_ν (the Miller regime);
  the existing kernel handles by retaining ample working precision.
  The slice-p1.4 sub-midpoint Maclaurin correction is the documented
  case.
- **Per-regime Ziv strategy**: drop-in wrap (already shipped).
  Slice p1.23 will need to confirm that the f32 sweep at
  bumped-precision can route through `certified_round_bf_to_f32`
  under directed modes (the helper's deliverable lifts the
  NE-only constraint on bumped precisions, per the pf-fwtz bead).
- **Cited spec**: DLMF 10.2 (Maclaurin), 10.17 (Hankel
  asymptotic); ADR-0036 for the property_jn dyadic constraint.
- **Oracle coverage**: MPFR-primary for J0, J1, Jn; differential +
  property + Wronskian cross-tie.
- **Estimated Ziv iterations at cap**: 1-3 on most inputs; the
  near-zero amplification regime may push to 4-5.
- **Worked example**: `J0(0)` = 1 (exact); `J1(0)` = 0 (exact);
  the Wronskian `J_n · Y_{n-1} − J_{n-1} · Y_n = 2/(π·x)` cross-ties.
- **Migration commit shape**: no new commits on the kernel; slice
  p1.23 ships the `certified_round_bf_to_f32` helper that lifts the
  NE-only constraint on the bumped precision path; the f32 sweep
  under directed modes then runs at p=320 without losing the
  directed-rounding signal in the bridge.

## Family p1.24 (1f.2): Elementary completions

### expm1

- **Source**: `src/math/expm1.rs:1-147`.
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard
  composition with `exp`.
- **eval(w) shape**: Short-circuit at `x.exponent ≤ -target - 8`
  (returns x rounded under `mode`). General path: `working_prec =
  target + 64 + cancellation`, where `cancellation = -exponent(x)`
  if `x.exponent < 0` (capped at `target + 1024`); round x to
  working_prec under NE; compute `e_x = exp(x_w)` under NE
  (internally Ziv-driven so e_x is correctly rounded at
  working_prec under NE); compute `diff = e_x - 1` under NE at
  working_prec; round diff to target under `mode`.
- **Cancellation regimes**: x near 0 — the `e_x - 1` cancellation
  loses `|exponent(x)|` bits, the existing cancellation boost
  exactly recovers them. The short-circuit handles the tiniest x
  where `x²/2` falls below half a ULP of x. No collapse-to-exact-zero.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. The
  existing cancellation boost moves INSIDE the eval(w) closure
  (eval(w) = `e^x_w - 1` at working_prec = `w + cancellation`).
  The current short-circuit stays before Ziv (per the slice p1.4
  precedent: short-circuits go BEFORE the driver, not inside).
- **Cited spec**: IEEE 754-2019 §9.2 expm1; standard cancellation
  boost technique (Higham, Accuracy and Stability of Numerical
  Algorithms §1.14).
- **Oracle coverage**: MPFR-primary (mpfr_expm1).
- **Estimated Ziv iterations at cap**: 1-2. The cancellation boost
  already compensates for the lost bits; Ziv adds the
  rounding-mode certificate.
- **Worked example**: `expm1(2^-50)` at p=53 — current path
  computes at working_prec = 53 + 64 + 50 = 167, returns
  `2^-50 + 2^-101 + …` rounded to NE. Under TowardPositive the
  current kernel rounds the same 167-bit working result to p=53
  under TP at the final step; this is correct to within 1 ULP at
  ties (the §1 caveat). Ziv wrap grows the working precision when
  the TP boundary lies in the uncertainty interval, certifying the
  result.
- **Migration commit shape**: 1 commit (kernel wrap + the
  short-circuit's documentation; unit tests pin against
  mpmath-derived values at p=53 under all five modes).

### log1p

- **Source**: `src/math/log1p.rs:1-167`.
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard
  composition with `ln`.
- **eval(w) shape**: Domain checks (`x = -1` → -∞ + DIV_BY_ZERO;
  `x < -1` → NaN + INVALID; `x = -∞` → NaN). Short-circuit at
  `x.exponent ≤ -target - 8`. General path: `working_prec = target
  + 64 + cancellation`, where `cancellation = -exponent(x)` if
  `x.exponent < 0`; round x to working_prec under NE; compute
  `one_plus_x = 1 + x_w` under NE; compute `ln_val = ln(one_plus_x)`
  under NE (Ziv-driven so ln is correctly rounded at working_prec);
  round ln_val to target under `mode`.
- **Cancellation regimes**: x near 0 — the `1 + x` cancellation
  loses `|exponent(x)|` bits, the existing cancellation boost
  exactly recovers them. Short-circuit handles the tiniest x. No
  collapse-to-exact-zero.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. eval(w) =
  `ln(1 + x_w)` at working_prec = `w + cancellation`. Short-circuit
  stays before Ziv.
- **Cited spec**: IEEE 754-2019 §9.2 log1p; standard cancellation
  boost.
- **Oracle coverage**: MPFR-primary (mpfr_log1p).
- **Estimated Ziv iterations at cap**: 1-2.
- **Worked example**: `log1p(2^-50)` at p=53 — current path
  computes at working_prec ≈ 167, returns ln(1 + 2^-50) ≈ 2^-50 -
  2^-101 + … under NE. Directed-mode pinning the same as expm1.
- **Migration commit shape**: 1 commit (kernel wrap + unit tests).

### exp2

- **Source**: `src/math/exp2.rs:1-121`.
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard composition
  with `exp` and `ln(2)`.
- **eval(w) shape**: `working_prec = target + 64`; round x to
  working_prec under NE; `ln_2 = ln_2_at(working_prec)` (the
  on-demand 1024-bit-capped constant from `agm_constants`); compute
  `product = x_w · ln_2` under NE; compute `result = exp(product)`
  under NE (Ziv-driven so correctly rounded at working_prec); round
  result to target under `mode`.
- **Cancellation regimes**: None on the documented domain — the
  composition `exp(x · ln(2))` has no cancellation. The exp path
  has no collapse-to-exact-zero on the documented domain (slice p1.2
  closed the underflow case).
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. eval(w) =
  `exp(x_w · ln_2_at(w))` at working precision w. No short-circuit
  needed.
- **Cited spec**: IEEE 754-2019 §9.2 exp2; standard composition.
- **Oracle coverage**: MPFR-primary (mpfr_exp2).
- **Estimated Ziv iterations at cap**: 1-2.
- **Worked example**: `exp2(0.5) = √2 ≈ 1.41421356...`; `exp2(10) =
  1024` (exact).
- **Migration commit shape**: 1 commit.

### exp10

- **Source**: `src/math/exp10.rs` (not yet read in this session;
  same shape as exp2 with `ln_10_at` instead of `ln_2_at`).
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard composition
  with `exp` and `ln(10)`.
- **eval(w) shape**: working_prec = target + 64; round x; ln_10 =
  ln_10_at(working_prec); product = x_w · ln_10 under NE; result =
  exp(product); round to target under mode. **AUDIT TBD on the
  exact source citation at slice p1.22 second session** — confirm
  the kernel shape matches the exp2 pattern.
- **Cancellation regimes**: None.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap.
- **Cited spec**: IEEE 754-2019 §9.2 exp10; standard composition.
- **Oracle coverage**: MPFR-primary (mpfr_exp10).
- **Estimated Ziv iterations at cap**: 1-2.
- **Worked example**: `exp10(3) = 1000` (exact).
- **Migration commit shape**: 1 commit.

### Family p1.24 ADR posture

The four migrations are all drop-in `ziv_round` wraps with the
existing short-circuits and cancellation boosts preserved inside
the eval(w) closure. No per-family ADR needed; the slice's kernel
doc comments record the change and the differential lanes confirm
the result. Slice p1.24's commit shape: 4 kernel-migration commits
+ 1 differential-lane-widening commit + 1 status-TOML-row update
commit + 1 caveats-§1-narrowing commit + 1 doc-comment-qualifier
commit = 8 commits total. Wall-clock estimate: 2-3 days.

## Family p1.25 (1f.3): Inverse trig

### asin

- **Source**: `src/math/asin.rs:70-151`.
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard.
- **eval(w) shape**: special-case dispatch (NaN, ±0, ±∞, |x|>1, x=±1
  exact-value short-circuits to ±π/2_at(target)); general path at
  working_prec = target + 64 under NE: round `|x|` to working_prec;
  compute `x_sq = |x|²`, `one_minus_sq = 1 − x_sq`, `s = sqrt(one_minus_sq)`,
  `denom = 1 + s`, `y = |x| / denom`; call `atan_finite_unsigned(y,
  working_prec)` (returns atan(|x|/denom) at working_prec under NE);
  compute `twice = 2 · atan_y`; sign-flip if x < 0; round to target
  under `mode`. The identity `asin(x) = 2·atan(x/(1+sqrt(1−x²)))`
  keeps `denom` bounded below by 1 for all `|x| ≤ 1`.
- **Cancellation regimes**: |x| near 1. `1 − x²` loses up to
  ~`|x|.precision` bits as |x| → 1, but the identity's algebraic
  structure absorbs the worst case: the divisor 1 + sqrt(1 − x²) is
  bounded between 1 (at x = ±1) and 2 (at x = 0), so y = x/denom
  stays bounded. The current fixed-64-bit guard handles this without
  collapse-to-exact-zero. The x = ±1 exact-value case short-circuits
  to ±π/2 before any cancellation occurs.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. eval(w) =
  the existing composition at working precision w. The x = ±1
  exact-value short-circuit stays before Ziv (returns the
  precomputed ±π/2 at target precision under mode). The special
  cases (NaN, ±0, ±∞, |x|>1) stay before Ziv too.
- **Cited spec**: IEEE 754-2019 §9.2 asin; standard half-angle
  identity (Abramowitz & Stegun 4.4.45).
- **Oracle coverage**: MPFR-primary (mpfr_asin).
- **Estimated Ziv iterations at cap**: 1-2. The cancellation in
  `1 − x²` is absorbed by the identity; the Taylor series in
  atan_finite_unsigned converges fast.
- **Worked example**: `asin(0.5) = π/6 ≈ 0.5235987...` —
  asin_kernel computes 2·atan(0.5/(1+sqrt(0.75))) = 2·atan(0.5/(1+
  0.8660...)) = 2·atan(0.2679...) = π/6 at working_prec = target +
  64. Under TowardPositive at p=24 the result rounds toward the
  upper f32 neighbor; the current kernel applies TP at the final
  step. Ziv certifies whether the working-precision result's
  uncertainty interval lies entirely on one side of the f32 tie.
- **Migration commit shape**: 1 commit (kernel wrap; unit tests pin
  against mpmath at p=53 under all five modes; the
  asin_sin_round_trip test stays).

### acos

- **Source**: `src/math/acos.rs:74-110` (head section read; full
  body identical pattern).
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard with
  two-branch composition.
- **eval(w) shape**: special-case dispatch (NaN, 0 → π/2, ±∞
  invalid). For x ≥ 0: `acos(x) = 2·atan(sqrt((1−x)/(1+x)))` at
  working_prec = target + 64. For x < 0: `acos(x) = π − 2·atan(sqrt(
  (1+x)/(1−x)))`. The two-branch form is the cancellation-resistant
  one: `π/2 − asin(x)` near x = 1 would lose ~`target` bits to
  cancellation; this form keeps the atan argument bounded.
- **Cancellation regimes**: x near +1 (first branch's `1 − x`
  loses bits but the algebraic structure recovers); x near −1
  (second branch's `1 + x` similar). The `acos(+1) = +0` and
  `acos(−1) = π` exact-value short-circuits handle the tightest
  cases. No collapse-to-exact-zero — the branch dispatch keeps
  each composition stable.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap, with the
  branch dispatch staying inside eval(w) (the branch choice is on
  the sign of x, which is a Class-level discriminant that doesn't
  depend on working precision). Special cases stay before Ziv.
- **Cited spec**: IEEE 754-2019 §9.2 acos; A&S 4.4.46 (the
  two-branch identity).
- **Oracle coverage**: MPFR-primary (mpfr_acos).
- **Estimated Ziv iterations at cap**: 1-2.
- **Worked example**: `acos(0.5) = π/3 ≈ 1.04719...` — first
  branch (x ≥ 0): 2·atan(sqrt(0.5/1.5)) = 2·atan(sqrt(1/3)) =
  2·atan(0.5773...) = 2·(π/6) = π/3.
- **Migration commit shape**: 1 commit.

### atan

- **Source**: `src/math/atan.rs:67-200` (header + atan_finite_unsigned).
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard with
  half-angle reduction + Taylor.
- **eval(w) shape**: special cases (NaN, ±0, ±∞ → ±π/2). General
  path at working_prec = target + 64 under NE: call
  `atan_finite_unsigned(|x|, working_prec)` which (i) inverts via
  `atan(|x|) = π/2 − atan(1/|x|)` if |x| > 1, (ii) applies half-
  angle `y ← y/(1 + sqrt(1+y²))` up to 64 times (cap line 147)
  until y.exponent < −4, (iii) sums the Taylor series `y − y³/3 +
  y⁵/5 − …` with `max_iter = 2·working_prec` terms (line 191), (iv)
  multiplies the sum by 2^k to reverse the half-angles. Sign-flip
  for negative x; round to target under `mode`.
- **Cancellation regimes**: |x| very small — atan(x) ≈ x with cubic
  correction `−x³/3`. For |x| < 2^(−working_prec/2) the cubic
  correction falls below the working-precision ULP; the Taylor
  series's first term is the answer at working precision and the
  remaining terms sum to ≤ |x|³/3 ≤ working_prec_ULP. No
  collapse-to-exact-zero — `atan_taylor`'s first term is `y`
  itself, which equals the input at zeroth half-angle iterations.
  The `should_halve(y)` predicate at line 174 returns false for
  exponent < −4, so tiny-x inputs skip the half-angle reduction
  entirely and the Taylor sum starts with `y` exactly.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. eval(w) =
  `atan_finite_unsigned(|x|, w)` with the sign applied outside.
  The Taylor max_iter cap of `2·working_prec` is the existing
  guarantee that the series converges below the working-precision
  ULP; under Ziv the working precision grows so the cap grows with
  it. The half-angle 64-iteration cap is fine for any reasonable
  target precision (each step shrinks |y| by ~2, so 64 iterations
  shrink by 2^64 from any starting |y| ≤ 1, far past any realistic
  precision target).
- **Cited spec**: IEEE 754-2019 §9.2 atan; standard half-angle
  reduction (atan(y) = 2·atan(y/(1+sqrt(1+y²)))).
- **Oracle coverage**: MPFR-primary (mpfr_atan).
- **Estimated Ziv iterations at cap**: 1-2.
- **Worked example**: `atan(1) = π/4 ≈ 0.7853981...` —
  atan_finite_unsigned(1, w): |1| = 1 (no inversion); half-angle
  shrinks until |y| < 1/16; Taylor sums; doubles back through 2^k.
- **Migration commit shape**: 1 commit.

### atan2

- **Source**: `src/math/atan2.rs:1-80` (header + dispatch table).
- **Status today**: NOT Ziv-wrapped. Quadrant dispatch + atan(y/x)
  composition at fixed-64-bit guard.
- **eval(w) shape**: Quadrant dispatch per the IEEE §9.2.1 table
  (lines 1-13 of the source). Special cases (NaN, signed zeros,
  ±∞ combinations) resolve to exact values at target precision
  (0, ±π/2, ±π, ±π/4, ±3π/4) using `pi_at` / `pi_over_2_at` at
  target. General finite case: working_prec = max(y.prec, x.prec)
  + 64; ratio = y/x at working_prec under NE; alpha =
  atan(ratio_w) at working_prec under NE (internally not Ziv-driven
  today; that's the dependency on p1.25 also migrating atan); for
  x < 0 add π or subtract π depending on sign of y. Round to target
  under `mode`.
- **Cancellation regimes**: x near 0 with y nonzero — the quadrant
  dispatch handles x = ±0 exact via the table; for x near 0 the
  division y/x has high magnitude and atan(y/x) → ±π/2, which the
  atan kernel handles smoothly. The `x < 0` branch's `π − atan(y/|x|)`
  composition can lose bits when atan(y/|x|) is near π, but only at
  large y/|x| where atan saturates to ±π/2; the subtraction π −
  (π/2) = π/2 is well-conditioned. No collapse-to-exact-zero on the
  finite-x finite-y path.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. eval(w) =
  the existing finite-case composition at working precision w
  (depends on atan being Ziv-driven at w under NE, which p1.25
  family migration achieves in the same slice). Exact-value
  short-circuits for special inputs stay before Ziv. Multi-arg eval:
  the closure captures both `y` and `x` as immutable references; no
  extension to the Ziv driver itself.
- **Cited spec**: IEEE 754-2019 §9.2.1 atan2 (the dispatch table).
- **Oracle coverage**: MPFR-primary for differential lane;
  property + worst-case-vector for the multi-arg surface (no f32
  sweep per docs/v1.0-surface.md "multi-argument defers").
- **Estimated Ziv iterations at cap**: 1-2.
- **Worked example**: `atan2(1, 1) = π/4`, `atan2(1, -1) = 3π/4`,
  `atan2(-1, -1) = -3π/4`. The quadrant dispatch + composed atan
  gives the correctly-signed-and-quadranted angle.
- **Migration commit shape**: 1 commit (multi-arg eval closure; the
  differential lane folds into the multi-arg surface confirmation
  at slice p1.36 — atan2's lane lives there per the multi-arg
  defers in the surface document).

### Family p1.25 ADR posture

All four kernels are drop-in `ziv_round` wraps with their existing
identities, branch dispatches, and special-case short-circuits
preserved inside eval(w) or before-Ziv. atan2 is multi-arg via
closure capture, no driver extension. No per-family ADR needed;
kernel doc comments record the change. Slice p1.25 commit shape:
4 kernel-migration commits + 1 differential-lane-widening commit
(asin/acos/atan single-arg lanes) + 1 status-TOML-row update
commit + 1 caveats-§1-narrowing commit + 1 doc-comment-qualifier
commit = 8 commits total. atan2's lane confirmation moves to
slice p1.36. Wall-clock estimate: 3-4 days (each kernel needs the
unit-test directed-mode pins against mpmath; the Taylor + half-
angle structure of atan is the trickiest convergence audit).

## Family p1.26 (1f.4): Forward trig

### sin

- **Source**: `src/math/sin.rs:76-160` (+ `sin_taylor` body).
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard with
  Payne-Hanek-style argument reduction.
- **eval(w) shape**: special cases (NaN, ±0, ±∞ → NaN+INVALID,
  range cap → NaN+INVALID). General path at working_prec = target +
  64: `Reduction { quadrant, r } = reduce(x, working_prec)` via the
  hardcoded 4096-bit `TWO_OVER_PI_LIMBS_4096` table (`trig_reduce.rs:
  11-29`) — multiplies x by 2/π, rounds to nearest integer to get
  quadrant mod 4, scales the fractional remainder by π/2 to get
  r ∈ [−π/4, π/4]. Quadrant dispatch: q=0 → sin_taylor(r), q=1 →
  cos_taylor(r), q=2 → −sin_taylor(r), q=3 → −cos_taylor(r).
  sin_taylor sums `r − r³/3! + r⁵/5! − …` with first term = r,
  iterating `term_{n+1} = −term_n · r²/((2n)(2n+1))` until the
  term's exponent < `−working_prec − 4` (max_iter cap of
  `2·working_prec`). Round to target under `mode`.
- **Cancellation regimes**: x near integer multiples of π. The
  reduction's r is non-zero because π is irrational and the 4096-bit
  table supplies enough bits; r ≈ x − k·π is well-defined at
  working_prec. sin_taylor's first term is `r` itself, which equals
  the input at zeroth iterations — does NOT collapse to exact zero
  unless r is itself exact zero (which only happens at x = ±0,
  handled by the Zero special case). At very large |x| past the
  4096-bit table budget the reduction returns None and the kernel
  short-circuits to NaN+INVALID before any Taylor evaluation;
  Ziv-driven working precision growth is bounded by the 1024-bit
  ZIV_GUARD_CAP, so max working = target + 1024 << 4096 and the
  cap is non-binding for any realistic v1.0 target precision.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. eval(w) =
  reduce(x, w) followed by the quadrant-dispatched Taylor at w.
  The range cap stays before Ziv (returns NaN+INVALID at target
  precision under any mode; no Ziv interval test possible on NaN).
  The special cases stay before Ziv.
- **Cited spec**: IEEE 754-2019 §9.2 sin; Payne–Hanek reduction
  (ACM TOMS 1983). The 4096-bit 2/π table is the existing
  TWO_OVER_PI_LIMBS_4096 constant.
- **Oracle coverage**: MPFR-primary (mpfr_sin).
- **Estimated Ziv iterations at cap**: 1-3. Near zeros of sin
  (r → 0 mod π) the amplification |f'/f| = |cot(r)| → ∞ which
  could push to 3-4 iterations; the 1024-bit cap accommodates this.
- **Worked example**: `sin(π/2) = 1` — reduce(π/2, w) returns
  quadrant=1, r=0; cos_taylor(0, w) = 1. `sin(2π)` — quadrant=0,
  r = 2π − 2·(π/2 rounded to nearest integer) ≈ 0 (small, non-zero
  at working precision); sin_taylor(r) = r at first order, true
  value 0; the Ziv interval test certifies whether the working-
  precision residual is below the rounding tie.
- **Migration commit shape**: 1 commit (kernel wrap; the range-cap
  early-return stays; unit tests pin against mpmath at p=53 under
  all five modes, including a near-2π input to exercise the
  reduction's residual r path).

### cos

- **Source**: `src/math/cos.rs` (same shape as sin via quadrant +1
  offset).
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard.
- **eval(w) shape**: identical to sin but quadrant dispatch is
  shifted by +1: q=0 → cos_taylor(r), q=1 → −sin_taylor(r), q=2 →
  −cos_taylor(r), q=3 → sin_taylor(r). cos_taylor sums `1 − r²/2! +
  r⁴/4! − …` starting at term_0 = 1. cos(±0) = 1; cos(±∞) =
  NaN+INVALID.
- **Cancellation regimes**: x near (2k+1)π/2 — cos = 0
  mathematically; reduced r near ±π/4 in some quadrant, and
  cos_taylor returns a non-zero value at working_prec near the true
  zero. The amplification |f'/f| = |tan(r)| → ∞ near r = ±π/2, but
  the reduction puts r in [−π/4, π/4] so cos_taylor's argument is
  bounded and the series converges. No collapse-to-exact-zero on
  the cos_taylor path (cos_taylor's first term is 1, not the
  potentially-tiny r²/2!).
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap.
- **Cited spec**: IEEE 754-2019 §9.2 cos.
- **Oracle coverage**: MPFR-primary (mpfr_cos).
- **Estimated Ziv iterations at cap**: 1-3.
- **Worked example**: `cos(0) = 1` (exact short-circuit via Zero);
  `cos(π/2) = 0` mathematically — reduce gives quadrant=1, r=0;
  −sin_taylor(0) = 0; the kernel returns 0 at target precision.
  At a non-zero-but-near-(2k+1)π/2 input the kernel returns a tiny
  non-zero value; Ziv certifies the rounding direction.
- **Migration commit shape**: 1 commit.

### tan

- **Source**: `src/math/tan.rs:73-120`.
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard ratio of
  sin and cos Taylor sums.
- **eval(w) shape**: special cases identical to sin/cos. General:
  reduce(x, working_prec); compute s = sin_taylor(r, working_prec),
  c = cos_taylor(r, working_prec); for q ∈ {0, 2} return s/c, for
  q ∈ {1, 3} return −c/s. Round to target under mode.
- **Cancellation regimes**: x near odd multiples of π/2 — tan
  diverges; the kernel returns a very large finite value (matches
  MPFR). The reduction places r near ±π/4 in quadrants 1 and 3
  where the denominator s = sin_taylor(r) is bounded away from
  zero (|sin(±π/4)| = √2/2). x near integer multiples of π —
  tan = 0 mathematically; the reduction gives r near 0 in quadrants
  0 and 2, and s = sin_taylor(r) returns r (small, non-zero at
  working precision); c = cos_taylor(r) returns ≈ 1; tan = s/c ≈
  r, well-conditioned. No collapse-to-exact-zero (the small-r case
  has s = r ≠ 0 at working precision, so s/c is the small but
  non-zero true value).
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. The ratio
  s/c runs at working precision under NE inside eval(w); the outer
  Ziv envelope handles the rounding mode. The near-pole regime (q
  ∈ {1, 3} with r near ±π/4 ⇒ sin near √2/2, cos near √2/2, ratio
  near ±1) is well-conditioned.
- **Cited spec**: IEEE 754-2019 §9.2 tan; standard ratio
  identity.
- **Oracle coverage**: MPFR-primary (mpfr_tan).
- **Estimated Ziv iterations at cap**: 1-3. Near asymptotes the
  amplification grows; the cap may bind on measure-zero exact-tie
  inputs (the documented Ziv cap caveat from ADR-0022).
- **Worked example**: `tan(π/4) = 1` mathematically — reduce gives
  quadrant=0, r=π/4 (at working precision); sin_taylor(π/4) =
  √2/2, cos_taylor(π/4) = √2/2; ratio = 1. `tan(π/3) ≈ 1.732…` —
  quadrant=0, r=π/3 (slightly more than π/4 mathematically but the
  reduction wraps it into [-π/4, π/4] via quadrant=1 offset, with
  r=-π/6 modulo the dispatch); actually wait, π/3 ≈ 1.047 which is
  > π/4 ≈ 0.785 — reduce would put it in quadrant=1 with r = π/3
  − π/2 = −π/6 ≈ −0.524. q=1 → tan = −c/s = −cos(−π/6)/sin(−π/6)
  = −(√3/2)/(−1/2) = √3 ≈ 1.732. ✓
- **Migration commit shape**: 1 commit.

### Family p1.26 ADR posture

All three kernels are drop-in `ziv_round` wraps. The
Payne–Hanek-style 4096-bit reduction table stays unchanged; the
1024-bit ZIV_GUARD_CAP is non-binding against the 4096-bit budget
for any realistic target precision. The range-cap special case
(input too large for the table) stays before Ziv as a NaN+INVALID
short-circuit. No per-family ADR needed. INTER-FAMILY DEPENDENCY:
slice p1.26 must land before slice p1.29 (gamma family) because
gamma's reflection branch composes through sin. Slice p1.26
commit shape: 3 kernel-migration commits + 1 differential-lane-
widening commit (sin/cos/tan) + 1 status-TOML-row update + 1
caveats-§1-narrowing + 1 doc-comment-qualifier = 7 commits.
Wall-clock estimate: 3-4 days (the near-zero amplification regime
needs careful unit-test pins; the 4096-bit reduction's edge
behavior at extreme |x| deserves a directed-mode property test).

## Family p1.27 (1f.5): Hyperbolic + inverse

### sinh

- **Source**: `src/math/sinh.rs:64-110`.
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard with
  `(expm1(x) − expm1(−x)) / 2` composition.
- **eval(w) shape**: special cases (NaN, ±0 → ±0, ±∞ → ±∞).
  General path at working_prec = target + 64: round x to working_prec;
  compute em1_pos = expm1(x_w), em1_neg = expm1(−x_w) under NE
  (each internally precision-boosted by the expm1 cancellation rule);
  diff = em1_pos − em1_neg; result = diff / 2; round to target under
  mode. The expm1-based form avoids the catastrophic cancellation
  the naive (exp(x) − exp(−x))/2 would suffer at |x| < 1, because
  expm1's internal cancellation boost recovers the leading bits in
  each individual call.
- **Cancellation regimes**: x near 0 — expm1(x) ≈ x and expm1(−x)
  ≈ −x, so diff ≈ 2x and result ≈ x. Each expm1 call handles its
  own cancellation regime; the outer subtraction is em1_pos −
  em1_neg = (x − x²/2 + …) − (−x − x²/2 + …) = 2x + … which has no
  leading-bit cancellation. The Zero special case handles x = ±0.
  No collapse-to-exact-zero: at any non-zero x, each expm1 returns
  a value with magnitude ≥ |x|, and the subtraction preserves
  ≥ |x| magnitude in diff.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. eval(w) =
  the existing composition at working precision w; the inner expm1
  calls run under NE at w (they need ≤ 2^24-ULP accuracy at w for
  the outer Ziv interval test, which the slice-3a fixed-guard NE
  expm1 already meets; p1.24's Ziv migration of expm1 makes the
  inner calls strictly NE-correctly-rounded but is not a strict
  prerequisite). No short-circuit before Ziv needed; the Zero case
  is handled by Class dispatch before eval(w) ever runs.
- **Cited spec**: IEEE 754-2019 §9.2 sinh; standard expm1-based
  identity (Higham §1.14).
- **Oracle coverage**: MPFR-primary (mpfr_sinh).
- **Estimated Ziv iterations at cap**: 1-2.
- **Worked example**: `sinh(1) ≈ 1.1752011936438014…` — em1_pos =
  e − 1 ≈ 1.7182…, em1_neg = e^(−1) − 1 ≈ −0.6321…, diff ≈ 2.3504…,
  result ≈ 1.1752…
- **Migration commit shape**: 1 commit.

### cosh

- **Source**: `src/math/cosh.rs:63-112`.
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard with
  `(exp(x) + exp(−x)) / 2` direct evaluation.
- **eval(w) shape**: special cases (NaN, ±0 → 1, ±∞ → +∞). General
  path at working_prec = target + 64: round x; compute e_pos =
  exp(x_w), e_neg = exp(−x_w) under NE; sum = e_pos + e_neg; result
  = sum / 2; round to target under mode. Both summands are
  positive, so the addition has no cancellation regardless of x's
  sign or magnitude.
- **Cancellation regimes**: NONE. cosh's identity is the
  cancellation-free composition (additivity of positive terms).
  cosh(±0) = 1 handled by Zero special case before eval(w) runs.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. eval(w) =
  the existing composition at working precision w; the inner exp
  calls inherit the Ziv-driven five-mode correctness from slice p1.2.
- **Cited spec**: IEEE 754-2019 §9.2 cosh.
- **Oracle coverage**: MPFR-primary (mpfr_cosh).
- **Estimated Ziv iterations at cap**: 1-2.
- **Worked example**: `cosh(0) = 1` (exact); `cosh(1) = (e + 1/e)/2
  ≈ 1.5430806348…`.
- **Migration commit shape**: 1 commit.

### asinh

- **Source**: `src/math/asinh.rs:66-121`.
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard with
  `log1p(|x| + |x|²/(sqrt(|x|² + 1) + 1)) · sign(x)`.
- **eval(w) shape**: special cases (NaN, ±0, ±∞). General path at
  working_prec = target + 64 on |x|: compute x_sq = |x|², x_sq_plus_one
  = |x|² + 1, s = sqrt(x_sq_plus_one), s_plus_one = s + 1,
  correction = |x|² / s_plus_one, arg = |x| + correction, lp =
  log1p(arg); sign-flip if x < 0; round to target under mode. The
  composition's algebraic structure ensures every term in the
  log1p argument is non-negative and s_plus_one ≥ 2 > 0; no
  cancellation.
- **Cancellation regimes**: NONE on the documented path. The
  identity was chosen to be cancellation-free for all |x|. The
  naive `ln(x + sqrt(x² + 1))` form for large negative x would
  cancel because x + sqrt(x² + 1) → 0⁺; computing on |x| and
  applying sign avoids it. At x = 0 the Zero special case handles
  the exact-value short-circuit.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. eval(w) =
  the existing composition at working precision w; the inner log1p
  call provides ≤ 2^24-ULP accuracy at w (sufficient pre- and
  post-p1.24).
- **Cited spec**: IEEE 754-2019 §9.2 asinh; standard log1p-based
  asinh identity.
- **Oracle coverage**: MPFR-primary (mpfr_asinh).
- **Estimated Ziv iterations at cap**: 1-2.
- **Worked example**: `asinh(1) = ln(1 + √2) ≈ 0.8813735…`.
- **Migration commit shape**: 1 commit.

### acosh

- **Source**: `src/math/acosh.rs:68-130` (header + body up to
  working_prec setup).
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard with
  `log1p((x − 1) + sqrt((x − 1)(x + 1)))`.
- **eval(w) shape**: special cases (NaN, ±∞ behaviors, x = 0
  invalid because 0 < 1, x = 1 → +0, x < 1 → NaN+INVALID). General
  path at working_prec = target + 64: compute xm1 = x − 1, xp1 = x
  + 1, prod = xm1 · xp1 (which is x² − 1), s = sqrt(prod), arg =
  xm1 + s, lp = log1p(arg); round to target under mode. The
  log1p-of-(xm1+s) form keeps the argument bounded above zero for
  x > 1 and going to 0 as x → 1⁺, where log1p(0) = 0 matches
  acosh(1) = 0.
- **Cancellation regimes**: x near 1 — the naive ln(x + sqrt(x² −
  1)) form collapses (because x + sqrt(x² − 1) → 1⁺ and ln(1) =
  0). The log1p((x − 1) + sqrt((x − 1)(x + 1))) form's argument
  goes to 0 smoothly, and log1p(0) = 0 is the limit. No
  collapse-to-exact-zero in the inner steps: xm1 is non-zero
  whenever x > 1 (x = 1 is the special-case short-circuit), and
  sqrt(prod) ≥ 0. arg is the sum of two non-negative quantities,
  one of which is positive.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. eval(w) =
  the existing composition at working precision w. The x = 1
  short-circuit stays before Ziv.
- **Cited spec**: IEEE 754-2019 §9.2 acosh; standard log1p-based
  acosh identity (avoids the near-1 cancellation of the naive
  form).
- **Oracle coverage**: MPFR-primary (mpfr_acosh).
- **Estimated Ziv iterations at cap**: 1-2.
- **Worked example**: `acosh(1) = 0` (short-circuit); `acosh(2) =
  ln(2 + √3) ≈ 1.3169579…`.
- **Migration commit shape**: 1 commit.

### atanh

- **Source**: `src/math/atanh.rs:68-130` (header + body up to
  working_prec setup).
- **Status today**: NOT Ziv-wrapped. Fixed-64-bit-guard with
  `(log1p(x) − log1p(−x)) / 2`.
- **eval(w) shape**: special cases (NaN, ±0, ±1 → ±∞+DIV_BY_ZERO,
  |x|>1 → NaN+INVALID, ±∞ invalid). General path at working_prec =
  target + 64: round x; lp_pos = log1p(x_w), lp_neg = log1p(−x_w);
  diff = lp_pos − lp_neg; result = diff / 2; round to target under
  mode. The log1p-based form handles small-x cancellation
  internally in each log1p call.
- **Cancellation regimes**: x near 0 — log1p(x) ≈ x − x²/2 + …,
  log1p(−x) ≈ −x − x²/2 + …, diff ≈ 2x + …, result ≈ x. No
  leading-bit cancellation in the subtraction (the constant and
  quadratic terms add rather than subtract). x near ±1 — diverges
  to ±∞; the special-case short-circuit at x = ±1 handles the
  exact case. For x near but not at ±1, log1p(−1 + ε) is large and
  negative (→ −∞), so diff is large; no cancellation. No
  collapse-to-exact-zero on the documented finite-x path.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap.
- **Cited spec**: IEEE 754-2019 §9.2 atanh; standard log1p-based
  identity.
- **Oracle coverage**: MPFR-primary (mpfr_atanh).
- **Estimated Ziv iterations at cap**: 1-2.
- **Worked example**: `atanh(0.5) = (1/2) · ln(3) ≈ 0.5493061…`.
- **Migration commit shape**: 1 commit.

### Family p1.27 ADR posture

All five kernels are drop-in `ziv_round` wraps. The existing
identities were chosen for cancellation resistance; no
short-circuits before Ziv are needed (the Zero/±1/±∞ special cases
stay before Ziv as exact-value dispatches). The inner expm1, exp,
log1p calls provide ≤ 2^24-ULP accuracy at working precision, well
under the Ziv error guard; the family is independent of p1.24's
elementary-completions migration (neither blocks the other). No
per-family ADR needed; kernel doc comments record the change.
Slice p1.27 commit shape: 5 kernel-migration commits + 1
differential-lane-widening commit + 1 status-TOML-row update +
1 caveats-§1-narrowing + 1 doc-comment-qualifier = 9 commits.
Wall-clock estimate: 4-5 days (five kernels, each with directed-
mode unit-test pins against mpmath at near-boundary inputs).

## Family p1.28 (1f.6): erf family (erfc migration)

`erf` is already Ziv-wrapped (covered in the already-Ziv cohort).
`erfc` is the remaining migration target.

### erfc

- **Source**: `src/math/erfc.rs:77-166` plus `erfc_asymptotic`
  body at `:178-` and the shared `erf_maclaurin` /
  `asymptotic_threshold_exponent` helpers from `erf.rs`.
- **Status today**: NOT Ziv-wrapped. Fixed-guard with two-regime
  dispatch and a regime-aware working-precision boost.
- **eval(w) shape**: special cases (NaN, ±0 → 1, +∞ → +0, −∞ → 2).
  For x < 0: reflection `erfc(−x) = 2 − erfc(|x|)` at working_prec
  = target + 8 (the subtraction loses at most one bit because
  erfc(|x|) ≤ 2). For x > 0: regime dispatch on x's binary exponent
  against `asymptotic_threshold_exponent(target)`. **Asymptotic
  regime** (x large enough that the divergent expansion's smallest
  term is below target ULP): `erfc(x) = (e^(−x²)/(x√π)) · Σ (−1)^k
  (2k−1)!!/(2x²)^k` truncated at the smallest term; working_prec =
  min(target + 64, target + 512). **Maclaurin regime** (small/moderate
  x): `1 − erf_maclaurin(x)` with working_prec = min(target + 192,
  target + 512); the +128 boost absorbs the cancellation of
  1 − erf(x) when erf(x) is close to 1 (e.g. at x=4, erf(x) ≈
  0.99999998 and 1 − erf(x) ≈ 1.5e−8 loses ~25 bits).
- **Cancellation regimes**: moderate x where erf(x) is close to 1
  (the Maclaurin regime's cancellation). The existing +128-bit
  boost handles erf(x) within 2^−128 of 1, which corresponds to
  x ≈ √(64·ln(2)) ≈ 6.6 for binary precisions up to ~128. Beyond
  that the asymptotic regime takes over (its truncation error is
  bounded by the smallest series term, independent of cancellation).
  Negative x: reflection at +8 bits absorbs the 2 − erfc(|x|)
  subtraction since erfc(|x|) ∈ [0, 2]. **No collapse-to-exact-zero**
  on the documented domain: the asymptotic envelope e^(−x²) is
  exponentially small but non-zero; the Maclaurin 1 − erf_maclaurin
  diff is the actual small value the kernel must compute, and
  erf_maclaurin is bounded above by 1 strictly for finite x.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. eval(w)
  carries the existing regime dispatch inside, with the +128-bit
  Maclaurin sub-boost preserved (the boost is INSIDE eval(w), not
  outside it; this is the slice p1.4 erf precedent). The
  asymptotic regime's truncation bound depends on x's magnitude;
  at very large x the smallest term sits far below working_prec
  ULP and the truncation is automatic. Negative-x reflection
  stays as the same composition (now nested through a Ziv-driven
  erfc(|x|) call).
- **Cited spec**: DLMF 7.6 (Maclaurin), 7.12 (asymptotic
  expansion); IEEE 754-2019 §9.4.3 erfc; the
  `asymptotic_threshold_exponent` heuristic per the existing
  `erf.rs` documentation.
- **Oracle coverage**: MPFR-primary (mpfr_erfc).
- **Estimated Ziv iterations at cap**: 1-3. The asymptotic regime
  near the threshold and the Maclaurin regime near the deepest
  cancellation are the higher-iteration cases.
- **Worked example**: `erfc(0) = 1` (exact short-circuit);
  `erfc(1) ≈ 0.1572992…` (small-x Maclaurin); `erfc(5) ≈
  1.5374598e−12` (asymptotic).
- **Migration commit shape**: 1 commit (kernel wrap; the
  reflection branch and regime dispatch preserved inside eval(w);
  unit tests pin against mpmath at strategic x values: x=0 exact,
  x=1 Maclaurin mid-range, x=4 deep-cancellation Maclaurin, x=5
  asymptotic threshold, x=10 deep-asymptotic).

### Family p1.28 ADR posture

Single-kernel slice (erf already five-mode). Drop-in Ziv wrap with
the existing regime dispatch and cancellation boost preserved
inside eval(w). No per-family ADR needed. Slice p1.28 commit
shape: 1 kernel-migration + 1 differential-lane-widening + 1
status-TOML-row update + 1 caveats-§1-narrowing (drops erfc from
the list) + 1 doc-comment-qualifier = 5 commits. Wall-clock
estimate: 1-2 days.

## Family p1.29 (1f.7): Gamma family

### gamma

- **Source**: `src/math/gamma.rs:78-149` (+ `gamma_sign_of` at
  `:156-`).
- **Status today**: NOT Ziv-wrapped. Fixed-guard composition
  `sign · exp(lgamma(x))`.
- **eval(w) shape**: special cases (NaN, ±0 → ±∞+DIV_BY_ZERO,
  negative integer pole → NaN+INVALID, +∞ → +∞, −∞ → NaN+INVALID).
  Positive-integer exact fast path: gamma(n) = (n−1)! when it fits
  in target precision. General path at working_prec = min(target +
  64, target + 512): ln_abs_gamma = lgamma_round(x, working_prec, NE)
  (Ziv-driven internally per the already-Ziv cohort); abs_gamma =
  exp(ln_abs_gamma) (Ziv-driven); result_sign = gamma_sign_of(x,
  working_prec) (composes sin(πx) for x < 0 non-integer); apply
  sign; round to target under mode. The lgamma+exp composition's
  errors are bounded by ≤ 2·2^−working_prec absolute on the
  positive-x path, well under ZIV_ERROR_GUARD.
- **Cancellation regimes**: x near non-positive integers — handled
  by special-case dispatch before eval(w). x where sin(πx) is near
  zero (negative non-integer near a pole) — gamma_sign_of's sign
  determination depends on `sin(πx)`'s sign, NOT its magnitude;
  the binary sign is robust away from the pole (where the
  negative-integer special case takes over). No
  collapse-to-exact-zero: positive-x branch has lgamma finite,
  exp finite-and-positive; negative-non-integer branch composes the
  same with a known-correct sign. Positive-integer exact fast path
  returns (n−1)! when it fits, which is exact and non-zero.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. eval(w) =
  the existing composition at working precision w. The positive-
  integer exact fast path stays before Ziv (returns the factorial
  at target precision under mode, exact for n−1)! ≤ target ULP). The
  gamma_sign_of call stays inside eval(w) at working precision w,
  but since it returns only a Sign enum (not a BigFloat), it
  doesn't participate in the Ziv interval test directly; the sign
  is applied AFTER the magnitude is rounded.
- **Cited spec**: DLMF 5.2 (poles and residues), 5.5.3 (reflection
  for sign), 5.7-5.11 (Stirling); IEEE 754-2019 §9.4 gamma.
- **Oracle coverage**: MPFR-primary (mpfr_gamma).
- **Estimated Ziv iterations at cap**: 1-3. Near gamma's zeros
  (which on the positive reals don't exist — gamma > 0 everywhere
  positive; on the negative reals there are no zeros either, just
  the alternating-sign minima between poles) the magnitude doesn't
  amplify the rounding mode boundary the way trigonometric kernels
  do. Most inputs converge in 1-2 iterations.
- **Worked example**: `gamma(5) = 24` (exact via the positive-integer
  fast path); `gamma(0.5) = √π ≈ 1.7724539…`; `gamma(−0.5) = −2√π
  ≈ −3.5449077…` (sign from gamma_sign_of: sin(−π/2) < 0 negates
  the positive |gamma(−0.5)|).
- **Migration commit shape**: 1 commit (kernel wrap; the positive-
  integer fast path and gamma_sign_of helper stay; unit tests pin
  against mpmath at positive integers, half-integers, and a near-
  pole negative-non-integer to exercise the sign path).

### digamma

- **Source**: `src/math/digamma.rs:70-195` (+ `z_min_for_target` at
  `:201-`).
- **Status today**: NOT Ziv-wrapped. Multi-branch fixed-guard.
- **eval(w) shape**: special cases (NaN, ±0 → −∞+DIV_BY_ZERO,
  negative integer pole → −∞+DIV_BY_ZERO, +∞ → +∞, −∞ → NaN+INVALID).
  General path at working_prec = min(target + 64, target + 512).
  **Negative non-integer branch**: reflection ψ(x) = ψ(1−x) − π·cot(πx);
  computes π, x, sin(πx), cos(πx), cot = cos/sin, π·cot, and the
  recursive digamma_round(1−x, working_prec, NE); subtracts. **Positive
  branch**: regime dispatch on x's exponent against z_min_for_target
  (target=53 yields z_min=64 roughly). For x ≥ z_min: direct
  stirling_digamma(x, working_prec). For 0 < x < z_min: shift by n
  = z_min − approx_x via the recurrence ψ(x) = ψ(x+n) − Σ_{k=0}^{n−1}
  1/(x+k); stirling_digamma on the shifted argument.
- **Cancellation regimes**: x near positive integer where ψ(x) = 0
  + small (e.g., ψ(1) = −γ ≈ −0.5772, ψ(2) = 1 − γ ≈ 0.4228) — no
  cancellation; the recurrence shift gives a well-conditioned
  result. Negative non-integer near a pole — the π·cot(πx) term has
  sin(πx) → 0, so cot blows up; the negative-integer special case
  handles the exact pole, but x slightly off-integer gives a large
  finite cot. The subtraction ψ(1−x) − π·cot(πx) doesn't cancel
  because the two terms have different magnitudes near the pole
  (ψ(1−x) ≈ ψ(large) ≈ ln(large), while π·cot is large; the
  divergence comes from π·cot, not from cancellation). No
  collapse-to-exact-zero on the documented domain.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap with the
  full branch dispatch inside eval(w). The recursive digamma_round
  call in the reflection branch needs careful handling: it routes
  through the positive branch (since 1−x > 1 when x < 0), so it
  doesn't recurse indefinitely. The recursive call's mode is NE at
  working_prec; the outer Ziv envelope drives the final mode.
- **Cited spec**: DLMF 5.5 (digamma); 5.11 (asymptotic series); the
  shift-recurrence form (Higham §16.3) and the reflection
  ψ(1−x) − ψ(x) = π·cot(πx) (DLMF 5.5.4).
- **Oracle coverage**: MPFR-primary (mpfr_digamma).
- **Estimated Ziv iterations at cap**: 1-3.
- **Worked example**: `digamma(1) = −γ ≈ −0.5772156649…`;
  `digamma(2) = 1 − γ ≈ 0.4227843350…`; the recurrence
  `ψ(x+1) = ψ(x) + 1/x` checks against `ψ(2) − ψ(1) = 1` exactly.
- **Migration commit shape**: 1 commit. Unit tests pin against
  mpmath at integer points (1, 2, 10), half-integers, and a
  near-pole negative-non-integer (e.g. −0.5 + 2^−40) to exercise
  the reflection path's large-cot magnitude regime.

### beta

- **Source**: `src/math/beta.rs:1-60` (header + impl signature). The
  case 4 O(1) factorial form per ADR-0030 / slice 8a.4c stays.
- **Status today**: NOT Ziv-wrapped. Multi-arg via composition
  `sign · exp(lgamma(a) + lgamma(b) − lgamma(a+b))` with full
  negative-domain case dispatch per ADR-0030.
- **eval(w) shape**: extensive special-case dispatch (NaN, infinite
  operands, both operands in Zle = {0, −1, −2, …}, a+b ∈ Zle with
  one operand negative integer and the other positive integer →
  the case 4 (−1)^m (m−1)!(n−m)!/n! evaluated through lgamma of
  three positive-integer factorials, case 4 O(1) closed form per
  ADR-0030 case-4 algorithm-cost lesson). General finite path at
  working_prec = max(a.prec, b.prec) + 64: compute ln_a =
  lgamma(a), ln_b = lgamma(b), ln_apb = lgamma(a+b) under NE
  (each Ziv-driven internally); ln_beta = ln_a + ln_b − ln_apb;
  abs_beta = exp(ln_beta); sign = product of three gamma_sign_of
  calls; apply sign; round to target under mode.
- **Cancellation regimes**: a+b near non-positive integer where
  lgamma diverges; handled by the a+b ∈ Zle special case. a+b such
  that lgamma(a) + lgamma(b) ≈ lgamma(a+b) (which happens generally
  — the gamma identity (a+b−1)!/((a−1)!(b−1)!) = (a+b−1)·…·a·1/b!
  approximately, with logarithmic terms canceling at scale) —
  there's cancellation in the ln_beta computation, but the magnitude
  of ln_beta itself is small (order ~ln(target)) so working_prec +
  64 is sufficient. The case-4 O(1) factorial form per ADR-0030
  closed-form lessons stays unchanged.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap with a
  multi-arg eval closure that captures both a and b. The Ziv
  driver's single-u32 eval signature is unchanged; the closure
  captures both BigFloats. The full ADR-0030 case dispatch stays
  before Ziv (special cases return exact values at target precision
  under mode; only the case-1 finite general path enters eval(w)).
- **Cited spec**: ADR-0030 (the negative-domain extension);
  DLMF 5.12 (beta); 5.5.3 (reflection for sign); 5.2 (poles and
  residues).
- **Oracle coverage**: differential against MPFR (no f32 sweep;
  multi-arg defers per docs/v1.0-surface.md); property + worst-case
  vectors at p1.36 multi-arg confirmation.
- **Estimated Ziv iterations at cap**: 1-3.
- **Worked example**: `beta(2, 3) = 1/12` (exact via positive-integer
  composition); `beta(0.5, 0.5) = π` (mathematically); `beta(−10, 4)
  = 1/840` (case-4 closed form per ADR-0030 hand-derivation pin).
- **Migration commit shape**: 1 commit. The ADR-0030 case dispatch
  is preserved unchanged; unit tests stay green (the
  beta_case4_factorial_exact_rational test and the
  beta_case4_large_m_terminates DoS-prevention test from slice
  8a.4c).

### Family p1.29 ADR posture

Three kernels, all drop-in `ziv_round` wraps. INTER-FAMILY
DEPENDENCY: gamma_sign_of for negative non-integer x calls sin(πx)
for sign determination; the SIGN is binary-stable away from poles,
so this is a soft dependency on p1.26 (sin), not a strict
prerequisite. The recursive digamma_round call in the reflection
branch routes to the positive branch (no infinite recursion). beta
preserves the ADR-0030 case dispatch and the case-4 O(1) closed
form unchanged. No per-family ADR needed; the three kernel doc
comments record the five-mode claim. Slice p1.29 commit shape:
3 kernel-migration commits + 1 differential-lane-widening (gamma,
digamma single-arg) + 1 status-TOML-row update + 1
caveats-§1-narrowing + 1 doc-comment-qualifier = 7 commits. beta's
multi-arg lane confirmation moves to slice p1.36. Wall-clock
estimate: 4-5 days (digamma's multi-branch dispatch needs the most
unit-test pinning).

## Family p1.30 (1f.8): Integrals

### Ei (exponential integral)

- **Source**: `src/math/ei.rs:1-80` (header + impl signature).
- **Status today**: NOT Ziv-wrapped. Two-regime fixed-guard.
- **eval(w) shape**: special cases (NaN, ±0 → −∞+DIV_BY_ZERO,
  +∞ → +∞, −∞ → −0). Regime dispatch on |x|'s binary exponent
  against an erf-style asymptotic threshold. **Small |x|**:
  convergent series `Ei(x) = γ + ln|x| + Σ_{k≥1} x^k/(k·k!)` with
  working precision boosted by ≈ |x|·log₂ e to absorb the
  alternating cancellation that dominates for x < 0 (the terms
  alternate while Ei(x) → 0⁻). **Large |x|**: divergent asymptotic
  `Ei(x) ∼ (e^x/x) · Σ_{k≥0} k!/x^k` summed to smallest term
  (optimal truncation near k ≈ |x|). Both regimes use γ (the
  Euler-Mascheroni constant from `euler_gamma_at`).
- **Cancellation regimes**: small x < 0 — alternating series whose
  partial sums lose bits as Ei(x) crosses zero between consecutive
  partial sums. The existing |x|·log₂ e working-precision boost
  recovers the lost bits (for x = −1, ~1.5 extra bits; for x =
  −10, ~14 extra bits). For x near 0, Ei → −∞; the +0 special
  case handles the exact pole. **No collapse-to-exact-zero on the
  finite-x path**: the partial-sum cancellation is bounded by the
  boost; the result is a well-defined finite value at working
  precision.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap with regime
  dispatch and cancellation boost preserved inside eval(w) per the
  erfc precedent. The inner ln|x|, e^x calls are Ziv-driven from
  earlier cohort migrations (correctly rounded under NE at working
  precision); the outer Ziv envelope drives the rounding mode.
- **Cited spec**: DLMF 6.2 (definition), 6.6.2 (convergent
  series), 6.12.2 (asymptotic).
- **Oracle coverage**: Arb-primary (no MPFR primitive for Ei in
  the public f32 sweep; the existing `differential_ei` runs
  against MPFR's mpfr_eint via the bindings).
- **Estimated Ziv iterations at cap**: 1-3. Near Ei's positive
  zero at x ≈ 0.3725 (`Ei(x) = 0` for x ≈ 0.37250741…) the
  amplification |f'/f| could push to 3-4 iterations on
  pathological inputs; the cap accommodates this.
- **Worked example**: `Ei(1) ≈ 1.8951178…`; `Ei(−1) ≈ −0.2193839…`
  (small-x series with cancellation boost); `Ei(10) ≈ 2492.2289…`
  (asymptotic).
- **Migration commit shape**: 1 commit.

### Si (sine integral)

- **Source**: `src/math/si.rs:1-60` (header + impl signature). The
  auxiliary functions `si_ci_f` / `si_ci_g` and the shared
  `asymptotic_threshold_exponent` live in this file too.
- **Status today**: NOT Ziv-wrapped. Two-regime fixed-guard, odd
  function computed on |x| with sign reapplied.
- **eval(w) shape**: special cases (NaN, ±0 → ±0, ±∞ → ±π/2).
  Regime dispatch on |x|. **Small |x|**: convergent alternating
  series `Si(x) = Σ_{k≥0} (−1)^k x^(2k+1)/((2k+1)·(2k+1)!)` with
  working precision boosted by ≈ |x|·log₂ e for the alternating
  cancellation. **Large |x|**: auxiliary-function form `Si(x) =
  π/2 − f(x)·cos(x) − g(x)·sin(x)` with f, g the shared asymptotic
  auxiliaries `si_ci_f`/`si_ci_g`, summed to their smallest term.
  Sign reapplied (Si is odd).
- **Cancellation regimes**: large |x| where Si(x) oscillates near
  ±π/2 — the asymptotic composition π/2 − f·cos − g·sin can
  cancel near the oscillation peaks (where Si(x) ≈ π/2 + small
  oscillation), but the magnitude of (f·cos + g·sin) is bounded
  by O(1/x) so the residual `π/2 − (π/2 ± O(1/x))` retains target
  precision. Small |x| alternating cancellation absorbed by the
  boost. No collapse-to-exact-zero: Si is entire and its zeros
  are isolated points (the first positive zero is at x = π);
  Si(π) ≈ 1.8519, not zero; Si has no real zeros.
  Wait — Si(x) ≠ 0 for any real x ≠ 0, because Si is monotonic
  increasing... actually let me re-check: Si(x) approaches π/2
  oscillating; it crosses through π/2 at certain points but isn't
  zero. Si(0) = 0 is the only zero, handled by the Zero special
  case.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap with regime
  dispatch and cancellation boost preserved inside eval(w). Inner
  sin, cos calls are Ziv-driven after p1.26 lands (or
  faithful-within-1-ULP under NE before, still under the Ziv error
  guard). SOFT inter-family dependency on p1.26 for the asymptotic
  regime's sin/cos composition.
- **Cited spec**: DLMF 6.2 (definition), 6.6.5 (convergent series),
  6.12.3 (asymptotic).
- **Oracle coverage**: Arb-primary (no MPFR primitive for Si).
- **Estimated Ziv iterations at cap**: 1-3.
- **Worked example**: `Si(1) ≈ 0.9460830…`; `Si(π) ≈ 1.8519370…`
  (close to the global max ≈ 1.8519+); `Si(10) ≈ 1.6583470…`
  (asymptotic regime).
- **Migration commit shape**: 1 commit.

### Ci (cosine integral)

- **Source**: `src/math/ci.rs:1-60` (header + impl signature).
- **Status today**: NOT Ziv-wrapped. Two-regime fixed-guard.
- **eval(w) shape**: special cases (NaN, +0 → −∞+DIV_BY_ZERO,
  +∞ → +0, x ≤ 0 → NaN+INVALID since Ci(-x) = Ci(x) - iπ is
  complex). Regime dispatch on x. **Small x**: convergent
  alternating series `Ci(x) = γ + ln(x) + Σ_{k≥1} (−1)^k
  x^(2k)/((2k)·(2k)!)` with working precision boost. **Large x**:
  auxiliary-function form `Ci(x) = f(x)·sin(x) − g(x)·cos(x)`
  using the shared si_ci_f/si_ci_g (no π/2 baseline like Si).
- **Cancellation regimes**: x near zeros of Ci (the first positive
  zero is x ≈ 0.6165) — the small-x branch's series and the
  large-x branch's f·sin − g·cos can both have cancellation near
  zero-crossings. The +0 special case handles the divergent end.
  No collapse-to-exact-zero on the small-x branch (the γ + ln(x)
  baseline is non-zero for any positive x ≠ 1; Ci(1) = γ + 0 + Σ…
  = γ - 1/2·1/2 + ... is non-zero anyway).
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap with regime
  dispatch and cancellation boost preserved inside eval(w). SOFT
  inter-family dependency on p1.26 for the asymptotic regime's
  sin/cos.
- **Cited spec**: DLMF 6.2 (definition), 6.6.6 (convergent series),
  6.12.4 (asymptotic).
- **Oracle coverage**: Arb-primary (no MPFR primitive for Ci).
- **Estimated Ziv iterations at cap**: 1-3. Near zeros of Ci (e.g.
  x ≈ 0.6165, 3.3842, 6.4272, …) the amplification could push to
  3-4 iterations.
- **Worked example**: `Ci(1) ≈ 0.3374039…`; `Ci(π) ≈ 0.0736680…`;
  `Ci(10) ≈ −0.0454564…`.
- **Migration commit shape**: 1 commit.

### li (logarithmic integral)

- **Source**: `src/math/li.rs:1-60` (header + impl signature).
- **Status today**: NOT Ziv-wrapped. Composition `li(x) = Ei(ln(x))`
  at boosted working precision.
- **eval(w) shape**: special cases (NaN, +0 → +0 with no flag (the
  defining integral over the empty interval; slice p1.6 fix), x = 1
  → −∞+DIV_BY_ZERO, +∞ → +∞, x ≤ 0 → NaN+INVALID). Slice p1.6's
  li-at-zero exact-bracket handling stays: the verifier's pinned-
  corpus entry for li(+0) = 0 certifies via the worker skipping the
  ±1 mantissa-unit padding when rad == 0. General path: t = ln(x)
  at working_prec; result = Ei(t) at working_prec; round to target
  under mode.
- **Cancellation regimes**: x near 1 — t = ln(1) = 0, Ei(0) = -∞;
  the x=1 special case handles the pole. x near li's positive
  zero at x ≈ 1.4514 (where li(x) = 0) — the ln(x) → 0.3725
  approximately, and Ei(0.3725) ≈ 0; the existing kernel boosts
  the inner Ei call's working precision sufficiently. The
  cancellation in Ei's small-x branch on the inner result handles
  this. Cross-cancellation between ln and Ei is bounded by the
  composition's working_prec boost.
- **Per-regime Ziv strategy**: drop-in `ziv_round` wrap. eval(w) =
  Ei(ln(x), w) at working precision w. Inner ln is Ziv-driven from
  the cohort; inner Ei is Ziv-driven after this family's migration
  (within slice p1.30 itself; ordering inside the slice: Ei first,
  then li). The slice p1.6 li-at-zero exact-bracket handling
  stays unchanged at the harness level (not in the kernel).
- **Cited spec**: DLMF 6.2.8 (li = Ei∘ln definition).
- **Oracle coverage**: Arb-primary (no MPFR primitive for li).
- **Estimated Ziv iterations at cap**: 1-3.
- **Worked example**: `li(2) ≈ 1.0451637…`; `li(10) ≈ 6.1655599…`;
  `li(1.4514) ≈ 0` (the Ramanujan-Soldner constant, the unique
  positive zero of li).
- **Migration commit shape**: 1 commit.

### Family p1.30 ADR posture

Four kernels; Ei, Si, Ci are independent two-regime drop-ins; li
composes through Ei within the same slice. SOFT inter-family
dependency: Si and Ci's asymptotic regime composes through sin/cos
which become five-mode-correct after p1.26 lands, but the
fixed-guard NE sin/cos already meets ZIV_ERROR_GUARD so the family
is not strictly blocked. The slice p1.6 li-at-+0 exact-bracket
handling stays at the harness level. No per-family ADR needed.
Slice p1.30 commit shape: 4 kernel-migration commits + 1
differential-lane-widening (Ei, Si, Ci, li) + 1 status-TOML-row
update + 1 caveats-§1-narrowing + 1 doc-comment-qualifier = 8
commits. Ordering inside the slice: Ei migrates before li (li
depends on Ei). Wall-clock estimate: 4-5 days (the asymptotic
auxiliary functions' truncation analysis needs care; the
cross-oracle agreement between Arb and mpmath on the alternating-
series cancellation regime deserves a directed-mode property
test).

## Family p1.31 (1f.9): Airy

**AUDIT TBD** — populated in slice p1.22 next session. Kernels:
`Ai`, `Bi`, `Ai_prime`, `Bi_prime`. Expected shape: Maclaurin for
small |x|; DLMF 9.7 asymptotic for large |x| (the recurrence sign
correction per slice 6n / ADR-0021 stands; per
`feedback_derive_dont_recall_coefficients` instance 3, the
`u_k` recurrence with the `(2k-1)·216·k` divisor is the
non-recalled form, with `u_1 = 5/72`, `u_2 = 3465/93312`). Wronskian
`Ai · Bi' − Ai' · Bi = 1/π` cross-tie reused.

## Family p1.32 (1f.10): Bessel J/Y

`bessel_j` family is already Ziv-wrapped (covered in the
already-Ziv cohort). `bessel_y` family is the remaining migration
target.

### bessel_y (Y0, Y1, Yn)

**AUDIT TBD** — populated in slice p1.22 next session. Expected
shape: forward recurrence for J then `Y_n = (J_n · cos(nπ) -
J_{-n}) / sin(nπ)` for non-integer order (n/a for the v1.0 surface
which restricts to integer n); for integer order the Hankel
asymptotic and the limit form at integer ν. J/Y Wronskian
`J_n · Y_{n-1} − J_{n-1} · Y_n = 2/(π·x)` cross-tie reused.

## Family p1.33 (1f.11): Bessel I/K

**AUDIT TBD** — populated in slice p1.22 next session. Kernels:
`I0`, `I1`, `In`, `K0`, `K1`, `Kn`. Expected shape: Miller's
backward recurrence for I (reuses J's `a_k(n)` coefficients);
asymptotic for K. **K-recurrence sign convention** per
`feedback_derive_dont_recall_coefficients` instance 4 / ADR-0025:
`K_{k+1} = (2k/x) K_k + K_{k-1}` (PLUS, opposite I, opposite the
naive read of the unified DLMF 10.29.1 cylinder relation). The
`recurrence_spot_check` test stays as the durable artifact.
**I-family** sub-midpoint cubic-Maclaurin precision bump
(BESSEL_TINY_VERIFICATION_PRECISION = 320 per slice p1.8) stands
and pre-dates Phase 1f. I/K cross-tie `I_ν · K'_ν − I'_ν · K_ν =
1/x` cross-tie reused.

## Family p1.34 (1f.12): Zeta (LAST)

**AUDIT TBD** — populated in slice p1.22 next session. Kernel:
`zeta`. **Inter-family dependency**: the functional equation
branch (DLMF 25.4.2, NOT 25.4.1 per
`feedback_derive_dont_recall_coefficients` instance 5 / ADR-0026)
composes `sin`, `pow`, `gamma`. Each must already be five-mode
correct; Phase 1f's ordering puts zeta last for this reason.
CVZ acceleration (Borwein "Algorithm 1") for positive `Re(s)`
stays; the algorithm pin via reproducing the paper's `n = 1, 2`
worked examples (`2a₀/3`, `(16a₀ − 8a₁)/17`) per ADR-0026 is the
durable artifact.

## Multi-arg confirmation (slice p1.36)

### pow, fma

Already five-mode correct (pow via ADR-0022; fma by IEEE
construction). Slice p1.36 confirms five-mode coverage via the
existing `differential_pow` and `differential_fma` lanes (the
former already exercises `BIT_EXACT_ROUNDING_MODES`; the latter
is in the slice p1.23 scaffolding sweep).

### atan2

Migrates within family p1.25 (inverse trig). Two-argument shape:
the audit (next session) records whether atan2 composes through a
five-mode-correct `atan` plus a domain reduction (which would
inherit) or requires its own `ziv_round` envelope.

### beta

Migrates within family p1.29 (gamma family). Two-argument shape;
sign/pole math per ADR-0030 stays; case-4 algorithm O(1) form per
slice 8a.4c stays. The five-mode pass widens
`differential_beta`.

### agm

**AUDIT TBD** — populated in slice p1.22 next session. AGM
converges quadratically (Brent-Salamin). Expected Ziv strategy:
single-arg-eval closure that captures both arguments, or a
multi-arg sibling helper in `ziv.rs` (audit decides). The
quadratic convergence means at the 1024-bit cap precision Ziv adds
~1 iteration over the existing path.

## Parametric Bessel N-ladder recommendation

**AUDIT TBD** — populated in slice p1.22 next session. Today's
status table has only `Jn_5.toml` (n=5). The audit will recommend
the N-ladder for `Jn`, `Yn`, `In`, `Kn` parametric sweep coverage
at slice p1.35. Candidate ladders include: `n ∈ {2, 5, 10}` (small
+ representative + Miller-regime test), `n ∈ {2, 5, 10, 25, 100}`
(wider coverage including the Miller-backward-recurrence stress
regime). The choice depends on the per-N sweep wall-clock cost
(measured at slice p1.22 against a representative N).

## Deferred from v1.0 surface (one-line entries; no audit work)

Per ADR-0038's "no extension to non-frozen-surface kernels"
clause, the following functions are deferred to post-v1.0 and
inherit Phase 1f's scaffolding when each lands:

- `lambert_w` — branch structure + Halley iteration; multi-month
  primary-source derivation work. Deferred.
- `incomplete_gamma`, `incomplete_beta` — tier-2 specials with
  multi-regime evaluation. Deferred.
- `bessel_yn` at non-integer order — full ν ∈ ℝ surface.
  Deferred.
- Hypergeometric forms (`hyper_0F1`, `hyper_1F1`, `hyper_2F1`) —
  multi-parameter convergence regions. Deferred.

Each will get its own ziv_round wrap, cancellation analysis,
oracle coverage, and five-mode sweep when implemented; the
scaffolding from slice p1.23 makes this drop-in for future work.

## Audit cross-references

- ADR-0038 ratifies the strategic commitment this audit derives
  the per-family slices against.
- `~/.claude/plans/phase-1f-dynamic-fog.md` is the phase plan
  ADR-0038 ratifies.
- The per-family slices (p1.24 through p1.34) execute against the
  per-kernel entries above; each entry's "Migration commit shape"
  field records the commit count the slice ships for that kernel.
- The audit is built incrementally: this document's first revision
  ships the structural skeleton, the per-kernel template, the
  already-Ziv cohort sanity entries, and the elementary
  completions family. Subsequent slice p1.22 sessions populate
  one family at a time, with the **AUDIT TBD** markers replaced
  by derived content from kernel-source reading.
