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

To be populated at slice p1g.1.

Per-kernel table (skeleton):

| Kernel | Source file | Exact-value subset | Subset predicate | Construction at target_precision | Mode-independence argument | Status |
|--------|-------------|--------------------|-------------------|---------------------------------|----------------------------|--------|
| gamma | src/math/gamma.rs:155 | positive integers n | `is_integer_test(x) && x.sign() == Positive && !x.is_zero()` | Iterate (n-1)! at target_precision via `BigFloat::mul` under NE; bail out on INEXACT | The exact factorial value is representable at target_precision iff iterated mul stays INEXACT-clean; mode irrelevant once the value is exact | shipped p1.29 |
| zeta | src/math/zeta.rs | negative odd integers (DLMF 25.6.3) | `is_negative_odd_integer(x)` (already negative-even has-exact path; pair with the negative-odd-rational path) | TODO p1g.1 — build numerator and denominator as exact i64 (e.g., n=1 → -1/12; n=3 → 1/120; n=5 → -1/252), divide at target_precision with INEXACT check | rational at any precision admitting the denominator | TODO p1g.1 |
| lgamma | src/math/lgamma.rs | x = 1 or x = 2 | `x == 1 || x == 2` (both lgamma values are exactly 0) | `BigFloat::try_new_zero(Sign::Positive, target_precision)` | lgamma(1) = lgamma(2) = 0 exactly under every mode | TODO p1g.1 |
| ... | ... | ... | ... | ... | ... | ... |

Full enumeration across the 47 v1.0-surface kernels lands at slice
p1g.1.

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
