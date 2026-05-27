# Phase 1g audit: verification architecture closure

This is the in-tree audit document ratifying ADR-0039's strategic
commitment. It walks the four-bead structure (pf-kk16, pf-yupm,
pf-tqzz, pf-hdh8), derives each kernel's exact-value subset and
ZIV_ERROR_GUARD bound from the kernel source (not recalled), and
records the per-kernel implementation derived from the analysis.

The audit is the load-bearing artifact for Phase 1g. Per ADR-0039's
scope-narrowing commitment (Kani discharge at fixed t ∈ {24, 53,
113}, structural-analogy claim for arbitrary t) and CLAUDE.md's
derive-don't-recall discipline, each per-kernel entry derives its
content from the kernel source citation, not from a similar
kernel's pattern. The 8a case-4 O(m) DoS precedent and the
pf-jn1y misdiagnosis precedent stay live: recalled error bounds
and recalled exact-value subsets ship as latent has-errors.

This document is built incrementally across the four Phase 1g
slices (p1g.1 through p1g.4). The structural skeleton lands in
slice p1g.0 (phase entry). Subsequent slices populate one section
at a time per their bead.

## Context

Phase 1f closed at merge `b08a831` (2026-05-26), making every
v1.0-surface kernel correctly rounded under all five IEEE 754-2019
rounding modes. ADR-0038 is the load-bearing decision.

The remaining gap, named in `DESIGN.md` "Caveats and open questions"
§1, is that `ZIV_ERROR_GUARD = 24` at `src/math/ziv.rs:59` is an
assumed internal-error bound, not a proven one. Phase 1g closes
that gap as a v1.0 blocker.

The kickoff prompt at
`docs/decisions/plans/verification-architecture-kickoff-prompt.md`
records the prior session's discussion that led to ADR-0039.

## Phase shape

Long-arc phase branch (`phase-1g-verification-closure`), Phase 1f
precedent. Per-slice sub-branches off the phase branch, unsigned
commits, fast-forward merge back at slice closure. Slice boundaries
are pause-to-debrief checkpoints. One signed merge from the phase
branch into `main` at phase closure (p1g.5). Prompt before the
YubiKey touch.

## Slice index

- **p1g.0** — phase entry: bead surgery, this doc skeleton, ADR-0039
  proposed.
- **p1g.1** — `pf-kk16` exact-value pre-Ziv dispatch audit. § Exact-value
  pre-Ziv dispatch audit below.
- **p1g.2** — `pf-yupm` per-function `ZIV_ERROR_GUARD` calibration.
  § Per-function `ZIV_ERROR_GUARD` calibration below.
- **p1g.3** — `pf-tqzz` Arb cross-check assertion. § Arb cross-check
  protocol extension below.
- **p1g.4** — `pf-hdh8` Kani Ziv interval-test soundness at
  t ∈ {24, 53, 113}. § Kani soundness theorem below.
- **p1g.5** — phase closure: ADR-0039 accepted, ADR-0038 amended,
  DESIGN.md Caveats §1 narrowed, README Verification posture
  tightened, signed merge.

## Exact-value pre-Ziv dispatch audit (p1g.1, pf-kk16)

The audit derives the defect class precondition explicitly:
the true mathematical value `f(x)` must be **exactly representable
at target precision** (integer, exact zero, or m / 2^k for non-
negative integers m, k). Rationals with odd denominators (e.g.
`-1/12 = ζ(-1)`) are binary-irrational and never lie on a target-
precision rounding boundary; the Ziv interval test certifies them
correctly without dispatch. The recall-able candidate list from
the kickoff prompt and `feedback_exact_value_defeats_ziv` was
re-derived against this precondition at p1g.1; several entries
that looked like candidates do not in fact satisfy it.

### Already-shipped (audit verifies coverage; no new code)

| Kernel | Source file | Exact subset | Dispatch shape |
|--------|-------------|--------------|----------------|
| gamma | src/math/gamma.rs:188-210 | positive integers n | `try_gamma_pos_integer_exact`: iterate (n-1)! at target_precision under NE, bail on INEXACT (slice p1.29) |
| Forward trig (sin/cos) at 0 | src/math/sin.rs, cos.rs | x = ±0 | Class::Zero arm returns ±0 or 1 |
| Hyperbolic (sinh/tanh/asinh/atanh) at 0 | src/math/{sinh,tanh,asinh,atanh}.rs | x = ±0 | Class::Zero arm |
| expm1(0) | src/math/expm1.rs | x = ±0 | Class::Zero arm |
| log1p(0) | src/math/log1p.rs | x = ±0 | Class::Zero arm |
| asin / acos / atan / atan2 boundary cases | src/math/{asin,acos,atan,atan2}.rs | x = ±1, etc. | `pi_at_round` / `pi_over_2_at_round` mode-aware constants (slice p1.25) |
| erf / erfc at boundaries | src/math/{erf,erfc}.rs | x = ±0, ±∞ | Class::Zero / Class::Infinity arms |
| Bessel J_n(0), I_n(0) | src/math/{bessel_j,bessel_i}.rs:184-199 | x = 0 | Class::Zero arm: J_0(0)=I_0(0)=1; J_n(0)=I_n(0)=0 (n ≠ 0) |
| zeta(0)=−1/2, zeta(−2n)=0, zeta(+∞)=1 | src/math/zeta.rs:153-191 | exact special cases | Inline pre-Ziv dispatches per DLMF 25.6 |
| Si(0)=0, Li(0)=0, Si(±∞)=±π/2 | src/math/{si,li}.rs | exact boundary cases | Class::Zero / Class::Infinity arms (+ p1.30 mode-aware π/2 for Si(±∞)) |
| `pow` 12 special cases | src/math/pow.rs | integer-y fast path + special cases | `try_pow_*` dispatches (slice 7c) |
| agm(x, 0), agm(0, x), agm(±∞, _) | src/math/agm.rs | zero / infinity operands | Inline dispatches |

### New dispatches shipped at p1g.1

Each of these is a transferable instance of the gamma(7) defect
class: the kernel's composition returns `T(x) + epsilon` at the
exact-value input, and directed modes tip rounding 1 ULP off the
exact value. Pre-dispatch returns the exact value mode-independent.

| Kernel | Source file | Exact subset | Dispatch shape | Commit |
|--------|-------------|--------------|----------------|--------|
| ln | src/math/ln.rs (in `ln_kernel`) | x = 1 | `try_from_i64_exact(1, target_precision)` equality test; return `try_new_zero(Sign::Positive, target_precision)` | p1g.1 |
| lgamma | src/math/lgamma.rs (`try_lgamma_small_pos_int_exact`) | x ∈ {1, 2} | Free function: equality test against 1 and 2 (Γ(1)=Γ(2)=1, so ln Γ=0); return exact +0 | p1g.1 |
| agm | src/math/agm.rs (in `agm_kernel`) | a == b (fixed point) | Equality test `a.partial_cmp(b) == Equal`; return `a.round_to_precision(target_precision, mode)` | p1g.1 |
| log2 | src/math/log2.rs (`power_of_two_exponent`) | x = 2^k for integer k | Free function: mantissa-bit-pattern check (raw top limb = `1u64 << 63`, all other limbs zero); return `try_from_i64_exact(k, target_precision)`; fall through if k doesn't fit | p1g.1 |

Each new dispatch has a directed-mode pinning test
(`*_under_every_directed_mode`) in the kernel's existing
`#[cfg(test)] mod tests` that exercises NE, TP, TM, TZ, RNA at
p ∈ {24, 53, 113}. The gamma(7) → 720.00000000000011 failure shape
is the template: assert the dispatch returns the exact value
under every directed mode.

### Rule-outs (candidates that look in-scope but are NOT in the defect class)

| Candidate | Rationale for rule-out |
|-----------|------------------------|
| ζ at negative odd integers (ζ(-1) = -1/12, ζ(-3) = 1/120, …) | Rational with denominators 12, 120, 252, … containing odd prime factors (Bernoulli-number denominators per the von Staudt-Clausen theorem). Binary-irrational, never on a target-precision rounding boundary. Ziv handles correctly. |
| log10(10^k) for integer k | `10 = 1010₂` is binary-irrational, so `10^k` is itself not exactly representable as a `BigFloat` input at any finite precision. Even when log10(10^k) = k is exact, the input `10^k` is approximate, so the composition's epsilon is at the precision of the input rather than the precision of the kernel; Ziv handles correctly within the 64-bit guard. |
| digamma(n) at positive integer n | `ψ(n) = -γ + Σ_{k=1}^{n-1} 1/k` involves the irrational Euler-Mascheroni constant and a harmonic sum with odd denominators. Binary-irrational. |
| beta(m, n) at positive integers m, n | `Β(m,n) = (m-1)!(n-1)! / (m+n-1)!`. The denominator `(m+n-1)!` factorial contains odd prime factors for m+n ≥ 3, so the rational is binary-irrational. (For m=n=1 the result is exactly 1, but that case is already covered by the `Class::Zero`/boundary dispatch in beta.rs.) |
| `pow` integer-y fast path | Already covered by slice 7c's `try_pow_*` dispatches (square-and-multiply at target precision, INEXACT-bail). |

### Discipline note

The audit pass is mechanical against the defect-class precondition.
A future addition to the v1.0 surface that introduces a new kernel
must re-run this audit against its `eval(w)` and special-case
table; the `KNOWN_CALIBRATED_KERNELS` table for pf-yupm at p1g.2
becomes a structural place to enforce that the pf-kk16 audit ran
as part of the calibration sign-off.

## Per-function `ZIV_ERROR_GUARD` calibration (p1g.2, pf-yupm)

The `ziv_round` driver carried a single global `ZIV_ERROR_GUARD =
24` at `src/math/ziv.rs:59` pre-Phase-1g, justified empirically by
the doc comment at `src/math/ziv.rs:50-58`. Phase 1g moves the
bound into a per-kernel calibrated value (`pub(crate) const`
per kernel in `src/math/ziv_calibration.rs`), forces every
`ziv_round` call site to pass an explicitly-named constant, and
prepares the structural ground for the pf-tqzz active sweep guard
(slice p1g.3) which will assert each kernel's calibrated bound
against the rigorous Arb midpoint on every f32 input.

### Driver signature

`ziv_round`'s `pub(crate)` signature gains a required
`error_guard: u32` fourth parameter. No implicit default: every
call site cites a per-kernel constant. The `KNOWN_CALIBRATED_KERNELS`
acceptance criterion 4 is met structurally by the compiler
(adding a new kernel that calls `ziv_round` without naming its
constant fails to type-check) and by the
`every_per_kernel_bound_fits_under_base_guard_margin` unit test
in `src/math/ziv_calibration.rs` (which enumerates the constants
and rules out anything ≥ `ZIV_BASE_GUARD - 16 = 48`).

### Per-kernel calibration table

All kernels at p1g.2 land at `DEFAULT_ERROR_GUARD = 24` by
algebraic analysis. The shared analysis template (per
`src/math/ziv.rs:50-58` pre-Phase-1g) is:

> Count floating-point operations on the `eval(w)` path; each
> NE-rounded op contributes ≤ 1 ULP of accumulated error; sum
> ≤ op_count ULPs (linear sum is pessimistic); upper bound in
> bits = ceil(log₂(op_count + safety)). Round up to the next
> power of two; almost every transcendental kernel lands well
> under 2²⁴ ULP at all working precisions this driver runs.

Per-kernel evidence cited in the kernel-side `ziv_calibration.rs`
doc comment. Empirical confirmation across the f32 grid is
delivered by pf-tqzz at p1g.3; any kernel whose calibrated bound
is violated at any swept input surfaces as a fatal report, the
constant widens to the smallest power-of-two passing the sweep,
and the provenance flips from `algebraic` to
`empirical (sweep at <commit>)`.

| Kernel | Source file | Constant | Bits | Provenance | Citation |
|--------|-------------|----------|------|------------|----------|
| exp | src/math/exp.rs:132 | `EXP_ERROR_GUARD` | 24 | algebraic | exp series ~4w iterations, sum ≤ 2¹⁴ ULP at 1024-bit cap |
| exp2 | src/math/exp2.rs:119 | `EXP2_ERROR_GUARD` | 24 | algebraic | composition `exp(x·ln(2))` |
| exp10 | src/math/exp10.rs:120 | `EXP10_ERROR_GUARD` | 24 | algebraic | composition `exp(x·ln(10))` |
| expm1 | src/math/expm1.rs:149 | `EXPM1_ERROR_GUARD` | 24 | algebraic | cancellation boost inside the eval closure (slice p1.24) |
| ln | src/math/ln.rs:167 | `LN_ERROR_GUARD` | 24 | algebraic | atanh series ~w/3 iterations |
| log1p | src/math/log1p.rs:169 | `LOG1P_ERROR_GUARD` | 24 | algebraic | atanh series with tiny-x boost (slice p1.24) |
| sin | src/math/sin.rs:135 | `SIN_ERROR_GUARD` | 24 | algebraic | Payne-Hanek reduction + quadrant Taylor ~w/2 |
| cos | src/math/cos.rs:121 | `COS_ERROR_GUARD` | 24 | algebraic | shared range reduction |
| tan | src/math/tan.rs:121 | `TAN_ERROR_GUARD` | 24 | algebraic | sin/cos composition |
| asin | src/math/asin.rs:142 | `ASIN_ERROR_GUARD` | 24 | algebraic | `2·atan(|x|/(1+sqrt(1-x²)))` (slice p1.25) |
| acos | src/math/acos.rs:157 | `ACOS_ERROR_GUARD` | 24 | algebraic | `π - 2·atan(sqrt((1+x)/(1-x)))` (slice p1.25) |
| atan | src/math/atan.rs:115 | `ATAN_ERROR_GUARD` | 24 | algebraic | unsigned composition on |x| (slice p1.25) |
| atan2 | src/math/atan2.rs:244 | `ATAN2_ERROR_GUARD` | 24 | algebraic | quadrant-shifted `atan(y/x)` (slice p1.25) |
| sinh | src/math/sinh.rs:106 | `SINH_ERROR_GUARD` | 24 | algebraic | `(expm1(x)−expm1(−x))/2` (slice p1.27) |
| cosh | src/math/cosh.rs:107 | `COSH_ERROR_GUARD` | 24 | algebraic | `(exp(x)+exp(−x))/2` |
| tanh | src/math/tanh.rs:122 | `TANH_ERROR_GUARD` | 24 | algebraic | composition through `tanh_at_w` (slice p1.27) |
| asinh | src/math/asinh.rs:111 | `ASINH_ERROR_GUARD` | 24 | algebraic | `log1p(|x| + x²/(1+sqrt(1+x²)))` |
| acosh | src/math/acosh.rs:143 | `ACOSH_ERROR_GUARD` | 24 | algebraic | `log1p((x−1) + sqrt((x−1)(x+1)))` |
| atanh | src/math/atanh.rs:131 | `ATANH_ERROR_GUARD` | 24 | algebraic | `(log1p(x) − log1p(−x))/2` |
| pow (`exp·ln`) | src/math/pow.rs:301 | `POW_ERROR_GUARD` | 24 | algebraic | `ln + mul + exp` composition; product bound ≪ 2²⁴ |
| pow (integer-y) | src/math/pow.rs:343 | `POW_INT_ERROR_GUARD` | 24 | algebraic | ~log₂(|n|) multiplications; n ≤ 2³¹ keeps sum ≤ 2⁵ ULP |
| gamma | src/math/gamma.rs:155 | `GAMMA_ERROR_GUARD` | 24 | algebraic | `sign(x)·exp(lgamma(x))`; integer-fast-path dispatches exactly (pf-kk16) |
| lgamma | src/math/lgamma.rs:155 | `LGAMMA_ERROR_GUARD` | 24 | algebraic | Spouge + reflection; ~30 ops at z_min≈20. Empirical confirmation pending pf-tqzz |
| digamma | src/math/digamma.rs:132 | `DIGAMMA_ERROR_GUARD` | 24 | algebraic | composition through `digamma_at_w` (slice p1.29) |
| beta | src/math/beta.rs:220, 304 | `BETA_ERROR_GUARD` | 24 | algebraic | `exp(lgamma(x)+lgamma(y)−lgamma(x+y))` composition |
| erf | src/math/erf.rs:133 | `ERF_ERROR_GUARD` | 24 | algebraic | asymptotic / Maclaurin dispatched at working precision (slice p1.4) |
| erfc | src/math/erfc.rs:142 | `ERFC_ERROR_GUARD` | 24 | algebraic | `1 − erf(...)` or direct asymptotic (slice p1.28) |
| Si | src/math/si.rs:136 | `SI_ERROR_GUARD` | 24 | algebraic | Maclaurin or asymptotic (slice p1.30) |
| Ci | src/math/ci.rs:145 | `CI_ERROR_GUARD` | 24 | algebraic | Maclaurin or asymptotic (slice p1.30) |
| Li | src/math/li.rs:143 | `LI_ERROR_GUARD` | 24 | algebraic | series summation (slice p1.30) |
| Ei | src/math/ei.rs:137 | `EI_ERROR_GUARD` | 24 | algebraic | series or asymptotic (slice p1.30) |
| Airy (Ai/Bi/Ai′/Bi′) | src/math/airy.rs:254 | `AIRY_ERROR_GUARD` | 24 | algebraic | shared eval body (slice p1.31) |
| Bessel J_n | src/math/bessel_j.rs:210 | `BESSEL_J_ERROR_GUARD` | 24 | algebraic | Maclaurin / Miller / asymptotic dispatched (slice p1.32; oracle bumps to p=320) |
| Bessel Y_n | src/math/bessel_y.rs:256 | `BESSEL_Y_ERROR_GUARD` | 24 | algebraic | reflection through J's; the oscillatory regime near J zeros may surface as the first empirical-widening candidate at p1g.3 |
| Bessel I_n | src/math/bessel_i.rs:234 | `BESSEL_I_ERROR_GUARD` | 24 | algebraic | Maclaurin / Miller / asymptotic (slice p1.33) |
| Bessel K_n | src/math/bessel_k.rs:251 | `BESSEL_K_ERROR_GUARD` | 24 | algebraic | reflection through I's (slice p1.33) |
| zeta | src/math/zeta.rs:204 | `ZETA_ERROR_GUARD` | 24 | algebraic | Borwein for s>0; FE composing `gamma·sin·pow·zeta_borwein` for s<0 (slice p1.34). Deepest composition on the surface; empirical confirmation pending pf-tqzz |
| agm | src/math/agm.rs:184 | `AGM_ERROR_GUARD` | 24 | algebraic | Gauss AGM iteration; quadratic convergence, ~log w ops |

39 call sites across 37 kernel modules; 38 distinct per-kernel
constants (pow has two paths — `exp·ln` and integer-y — each with
its own constant; beta has two `ziv_round` calls but both inside
`beta_kernel` so they share `BETA_ERROR_GUARD`). The
`every_per_kernel_bound_fits_under_base_guard_margin` and
`calibration_table_enumerates_expected_kernel_count` tests in
`src/math/ziv_calibration.rs` enforce the count drift guard.

### Risk-mitigation note

Any kernel whose calibrated bound exceeds `ZIV_BASE_GUARD - 16 =
48` (currently `ZIV_BASE_GUARD = 64` at `src/math/ziv.rs:37`)
triggers a paired bump to `ZIV_BASE_GUARD` (likely to 96),
recorded in this doc as a paired decision. At Phase 1g landing no
kernel hits this threshold.

## Arb cross-check protocol extension (p1g.3, pf-tqzz)

### Status at p1g.3 parts 1 & 2 landing

**Landed:**

- **Driver foundation (p1g.3 part 1, commit 35e3ee3).** New
  `ziv_round_capturing` in `src/math/ziv.rs` returns
  `ZivTrace = (BigFloat, Status, u32, BigFloat)` exposing the
  converged working precision and the eval(w) intermediate. The
  pre-existing `ziv_round` becomes a thin wrapper that destructures
  with `_` for the trailing pair; the 39 existing call sites stay
  unchanged.

- **Worker protocol extension (p1g.3 part 2).** New `MIDPOINT` verb
  in `scripts/arb_oracle_worker.py` reads
  `MIDPOINT <fn_id> <order_or_dash> <input_hex> <oracle_prec>`,
  computes the function at `ctx.prec = oracle_prec` via `python-
  flint` ball arithmetic, and returns the ball midpoint via
  `arf.man_exp()` as a lossless triple
  `OK <sign> <mantissa_hex> <exponent>` (or `INC` for non-finite,
  `ERR <msg>` on failure). The verb dispatches on the first
  request token; the original 4-token implicit "CERTIFY" form
  stays backward-compatible.

- **Rust parser and oracle method (p1g.3 part 2).** New
  `ArbOracle::midpoint(f, input, oracle_prec) -> Result<Float,
  MidpointError>` in `tests/oracle/arb.rs` parses the wire format
  back into a `rug::Float` at `oracle_prec`: hex-encoded
  absolute mantissa via `rug::Integer::parse_radix`, signed lift
  via `Float::with_val(prec, &signed)`, scaled by `<< exp` or `>>
  -exp` for the binary exponent. Mode-independent (the midpoint
  request has no mode parameter; the mode applies downstream when
  the cross-check assertion rounds the gap).

- **End-to-end smoke (p1g.3 part 2).** New
  `tests/oracle_arb_midpoint_smoke.rs` exercises the wire format
  on three Arb-primary kernels: `Si(0) = 0` (zero-encoded wire
  form), `Si(1) ≈ 0.946083070367183` (NIST DLMF 6.7.1 reference,
  matches within 1e-13 f64 tolerance at oracle_prec=128), and
  `K_0(1) ≈ 0.421...` (NIST DLMF 10.32.9 sanity, finite-positive
  range check). All three pass under the Arb venv at
  `${HOME}/.cache/pfloat-arb-oracle/venv`.

### Remaining for full pf-tqzz acceptance (follow-up sub-slice)

The wire format is validated; the cross-check sweep harness
remains. The deferred work:

1. **Per-kernel `<fn>_round_capturing` wrappers**, one per
   five-mode-correct `FnId`. Each kernel function gets a thin
   `#[cfg(any(test, feature = "ziv-instrumented"))] pub fn
   <fn>_round_capturing(...) -> ZivTrace` that calls
   `ziv_round_capturing` with the same per-kernel `error_guard`
   the production path uses. Mechanical 47-kernel pass.

2. **MPFR-side midpoint** for the 35 MPFR-primary kernels (the
   `MIDPOINT` verb only knows the 12 Arb-primary `FnId`s today).
   Add a `MpfrOracle::midpoint` method that calls
   `mpfr_<fn>(..., RoundNearest)` at `oracle_prec` and returns
   the result as `rug::Float`. MPFR's correct-rounding at
   `oracle_prec >= working_prec + 64` is itself the rigorous
   midpoint (no ball-radius adjustment needed because MPFR's
   directed rounding at `oracle_prec` already brackets the true
   value within sub-ULP).

3. **`tests/oracle/cross_check.rs` sweep harness**. For each
   `(kernel, input, mode)` triple in the 65536 × 5 × 47 sweep:
   call `<fn>_round_capturing` to obtain `(_, _, working,
   eval_w)`; call `oracle.midpoint(f, input, working + 64)` to
   obtain `arb_mid` (or MPFR-mid); compute
   `error = |eval_w - arb_mid|` and `bound = 2^(error_guard -
   working) * |arb_mid|`; assert `error <= bound`. Fail-fast
   structured report on any violation.

4. **Cargo.toml `ziv-instrumented` feature** to gate the per-
   kernel capturing wrappers behind a release-time-only flag.

Estimated effort: 2-4 hours for items 1-3 (mechanical), bounded
by per-kernel `pub fn` boilerplate generation. The 3.1M Arb
midpoint calls (~one per kernel-input pair, mode-independent)
fits the per-release runtime budget; per-push CI unchanged.

### Protocol reference (frozen)

| Direction | Wire format |
|-----------|-------------|
| Request   | `MIDPOINT <fn_id> <order_or_dash> <input_hex> <oracle_prec>` |
| Response (success) | `OK <sign> <mantissa_hex> <exponent>` — value = `sign * mantissa * 2^exponent` |
| Response (zero) | `OK + 0 0` — exact zero |
| Response (non-finite) | `INC` — NaN or unbounded ball |
| Response (error) | `ERR <message>` — worker-side failure |

The `<sign>` token is `+` or `-`; `<mantissa_hex>` is the absolute
integer mantissa as lowercase hex with no `0x` prefix; `<exponent>`
is signed decimal. The triple is a faithful representation of the
Arb ball midpoint at `oracle_prec` precision.

## Kani soundness theorem (p1g.4, pf-hdh8)

To be populated at slice p1g.4.

Theorem (Phase 1g form, at fixed t):

```
For all BigFloat y at precision w, BigFloat h at precision w with h ≥ 0,
for all RoundingMode m, for fixed t ∈ {24, 53, 113}:

  if round_to_precision(y − h, t, m) == round_to_precision(y + h, t, m)
  then for all BigFloat y' with |y' − y| ≤ h:
       round_to_precision(y', t, m) == round_to_precision(y, t, m)
```

Bounded encoding: `BoundedBigFloat<80>` (LIMBS=80 covers 5120
working bits = max(target=4096) + ZIV_GUARD_CAP=1024). Fixed-size
`[u64; LIMBS]` mantissa unrolls under CBMC where `Vec<u64>` does
not (ADR-0012 lesson).

Conversion-shim soundness: `BoundedBigFloat<80>` ↔ `BigFloat`
shims Kani-checked.

Per-precision proof output (Kani CBMC output, time-to-solve,
counterexample status) lands at slice p1g.4.

## Closure prose (p1g.5)

To be populated at slice p1g.5.

- ADR-0039 status `proposed` → `accepted`; final form names the
  three Kani-discharged target precisions and the structural-
  analogy claim for arbitrary t.
- ADR-0038 Consequences amendment (1 paragraph).
- DESIGN.md Caveats §1 narrowing (the empirical-slack framing at
  `src/math/ziv.rs:50-58` becomes a cross-reference to
  `ziv_calibration.rs`).
- README Verification posture (lines 189–194) tightened to
  reference ADR-0039 + sweep cross-check + Kani-discharged
  soundness at IEEE binary32/64/128. Protected disclosure block
  at "How pfloat is developed" stays bit-identical
  (`feedback_disclosure_update_under_explicit_permission`;
  verify `git apply --check docs/disclosure-correction-v1.0.diff`).
- Close pf-kk16, pf-yupm, pf-tqzz, pf-hdh8.
- Prompt before YubiKey; signed merge
  `phase-1g-verification-closure` → `main`.
