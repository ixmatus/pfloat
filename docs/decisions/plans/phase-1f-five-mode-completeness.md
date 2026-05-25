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

**AUDIT TBD** — populated in slice p1.22 next session. Kernels:
`asin`, `acos`, `atan`, `atan2`. Expected cancellation regimes:
asin/acos near ±1 (the `sqrt(1 - x²)` term loses bits as |x| → 1);
atan2 quadrant boundaries.

## Family p1.26 (1f.4): Forward trig

**AUDIT TBD** — populated in slice p1.22 next session. Kernels:
`sin`, `cos`, `tan`. Expected shape: Payne-Hanek-style argument
reduction via `trig_reduce.rs` (4096-bit working-precision cap);
Maclaurin series on reduced argument. Cancellation regimes: tan
near odd multiples of π/2; large arguments where the reduction
table's 4096-bit cap is the precision floor.

## Family p1.27 (1f.5): Hyperbolic + inverse

**AUDIT TBD** — populated in slice p1.22 next session. Kernels:
`sinh`, `cosh`, `asinh`, `acosh`, `atanh`. Expected cancellation
regimes: sinh/cosh near 0 (Taylor expansion required, follow tanh
short-circuit precedent); asinh near 0 (`ln(x + √(x² + 1))`
cancellation when x small); acosh near 1 (`ln(x + √(x² - 1))`
cancellation when x → 1⁺); atanh near 0 (Taylor) and near ±1
(infinite gradient).

## Family p1.28 (1f.6): erf family (mode widening only)

`erf` is already Ziv-wrapped (covered in the already-Ziv cohort).
`erfc` is the remaining migration target.

### erfc

**AUDIT TBD** — populated in slice p1.22 next session. Expected
shape: short-circuit `erfc(x) = 1 - erf(x)` for small |x|;
continued-fraction or asymptotic for large |x|. Cancellation
regime: `1 - erf(x)` near `x = 0` requires the asymptotic branch;
no collapse-to-exact-zero expected with the existing branch logic.

## Family p1.29 (1f.7): Gamma family

**AUDIT TBD** — populated in slice p1.22 next session. Kernels:
`gamma`, `digamma`, `beta`. **Inter-family dependency**: gamma's
reflection branch composes through `sin`, which must already be
five-mode correct (family p1.26 lands first). beta sign/pole math
per ADR-0030 stands; the case-4 algorithm cost fix (ADR-0030
addendum, slice 8a.4c) stays. digamma special cases per existing
kernel.

## Family p1.30 (1f.8): Integrals

**AUDIT TBD** — populated in slice p1.22 next session. Kernels:
`Ei`, `Si`, `Ci`, `li`. Expected cancellation regimes: asymptotic
branches of each. `li(x) = Ei(ln(x))` composition stays. Slice
p1.6 li-at-+0 exact-bracket handling preserved.

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
