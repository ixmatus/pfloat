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

To be populated at slice p1g.2.

Per-kernel table (skeleton):

| Kernel | Source file | Calibrated bound (bits) | Provenance | Citation |
|--------|-------------|-------------------------|------------|----------|
| exp | src/math/exp.rs:132 | 24 (default) | algebraic | exp series at ~4w iterations, each ≤ 1 ULP NE-rounding error, sum well under 2^24 ULP per the ziv.rs:51-58 analysis template |
| ln | src/math/ln.rs:148 | 24 (default) | algebraic | atanh series at ~w iterations |
| lgamma | src/math/lgamma.rs:141 | TODO p1g.2 | TODO | Stirling + reflection + recurrence depth |
| bessel_y | src/math/bessel_y.rs:256 | TODO p1g.2 | TODO (empirical likely; oscillatory regime near zeros) | widened sweep at <commit> |
| ... | ... | ... | ... | ... |

Full enumeration across the 44 `ziv_round` call sites lands at
slice p1g.2.

## Arb cross-check protocol extension (p1g.3, pf-tqzz)

To be populated at slice p1g.3.

Protocol extension (skeleton):

- **New worker verb:** `MIDPOINT <fn_id> <order_or_dash> <input_hex> <working_prec> <oracle_prec>`.
- **Worker computation:** `python-flint` ball arithmetic at `oracle_prec ≥ working_prec + 64`; return ball midpoint as a lossless BigFloat triple.
- **Response shape:** `OK <sign_hex> <exp_hex> <mantissa_hex>` (lossless), `INC`, or `ERR <msg>`.
- **Spike kernel:** `exp` (simplest, well-characterized).
- **Sweep cost:** ~3.1M new Arb calls per release (mode-independent midpoint, one per `(kernel, input)`).

Spike output and per-kernel cross-check results land at slice p1g.3.

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
