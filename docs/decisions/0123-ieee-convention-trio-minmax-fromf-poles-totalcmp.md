# ADR-0123: IEEE convention trio — 754-2019 minimumNumber/maximumNumber, lossless from_f32/f64, log-pole ±0, and LSB-anchored payload order

- **Status**: accepted
- **Date**: 2026-07-03

## Context

Four convention defects from the 2026-06 review (epic pf-8iji), each a
place where pfloat disagreed with a standard (or with itself). The
conventions were decided 2026-06-10 with verbatim primary-source
citations (`~/.claude/plans/pfloat-convention-citations-2026-06-10.md`).

1. **pf-l5en — min/max were 754-2008 minNum, not 754-2019 minimumNumber.**
   `min(+0, −0)` was order-dependent, and `min(sNaN, 1)` returned a quiet
   NaN + INVALID. 754-2019 §9.6 **withdrew** minNum/maxNum precisely
   because the sNaN-beats-number rule is non-associative.
2. **pf-jd4s — from_f32/f64 quieted signaling NaNs and dropped payloads**,
   contradicting their "lossless" docstring.
3. **pf-k8ax — the ±0 pole convention was inconsistent.** `ln(±0)` was
   correct (−∞ pole both signs) but `y0(−0)`/`k0(−0)` returned NaN +
   INVALID (grouping −0 with the negative axis), and the Ci doc claimed
   −0 → NaN.
4. **pf-qe5c — total_cmp compared NaN payloads MSB-aligned** though they
   are LSB-anchored, so equal payloads at different precisions compared
   unequal.

## Decision

1. **pf-l5en — implement 754-2019 §9.6 minimumNumber/maximumNumber**
   (C23 N3088 7.12.12.8/9, F.10.9.5): a NaN argument is missing data —
   return the NUMBER when the other is a number; a signaling-NaN argument
   raises INVALID (even though the number is returned); a both-NaN pair
   returns a quiet NaN (INVALID iff either is signaling). Signed zeros are
   ordered: `−0 < +0` (N3088 7.12.12.4/5, carried to the `_num` variants
   by "differ only in their treatment of NaN arguments"). This matches
   C23 fminimum_num/fmaximum_num, RISC-V FMIN/FMAX, and LLVM
   minimumnum/maximumnum. A NaN-propagating minimum/maximum pair (§9.6's
   other two operations) is a separate additive v1.x decision.

2. **pf-jd4s — from_f32/f64 are a lossless representation embed** in the
   IEEE §5.5.1 quiet-copy class (like copysign/abs, which "raise no
   floating-point exception, even if x is a signaling NaN"): preserve the
   signaling bit AND the payload verbatim (f32: 22 payload bits, f64: 51,
   LSB-anchored), raise nothing. Justified by the no-Status signature (the
   function cannot honestly signal) and by BigFloat representing
   sNaN+payload exactly. This is a **documented divergence** from IEC 60559
   convertFormat (which DOES signal on sNaN), not a claim IEEE blesses
   quiet conversion; the arithmetic direction `to_f32/f64_round` keeps
   §7.2 semantics (sNaN → INVALID + quiet), giving a coherent asymmetry:
   exact embed in, arithmetic conversion out.

3. **pf-k8ax — ±0 belongs to the POLE for every log/power-singularity
   kernel**, for BOTH zero signs: `f(±0) = pole + DIV_BY_ZERO`; strictly
   negative `x` → NaN + INVALID. Grounds: C11 F.10.3.7 (`log(±0) = −∞ +
   divideByZero` for both zero signs) and IEEE 754-2019 §9.2; `−0.0 < 0.0`
   is false in C, so "domain error for x < 0" cannot capture `−0`, and
   POSIX y0's "if x is 0.0, pole error" therefore covers `−0`; Ci ~ γ +
   ln x and Y0/K0 ~ ln x share ln's singularity class. Fixes: `y0(−0)` and
   `k0(−0)` → the pole (was NaN + INVALID); the Ci/Y/K module docs. Pole
   sign by order parity: `Y_n(±0) = −∞` unless `n` is negative AND odd,
   then `+∞` (`Y_{−n} = (−1)^n Y_n`, DLMF 10.4.1; POSIX yn verbatim) — this
   also fixed a latent unconditional-−∞ bug for negative-odd orders;
   `K_n(±0) = +∞` for every order (`K_{−n} = K_n`, positive).

4. **pf-qe5c — compare NaN payloads LSB-anchored** (`limb_cmp_lsb_aligned`,
   indexing from limb 0 with missing high limbs zero), keeping the
   MSB-aligned `limb_cmp_aligned` for left-aligned mantissas. Equal
   payloads at different precisions now compare equal.

## Consequences

- pfloat's comparison, conversion, and pole conventions match the cited
  standards and are internally consistent. Existing in-module tests that
  encoded the withdrawn minNum, the quieting conversion, and the −0→NaN
  poles were updated to the new semantics (they asserted the defects).
- Full lib suite green (879); `differential_yn`/`differential_ik`
  unchanged.

### Inversion (failure paragraphs considered)

- *"min/max: re-cite the doc to match the code."* Rejected: the code was
  754-2008 minNum, which 754-2019 withdrew for non-associativity; the doc
  already cited §9.6, so the code was wrong, and §9.6 is what modern
  hardware/compilers implement.
- *"from_f32: a secret thread-local raise on sNaN."* Rejected: worst
  option — impossible under no_std, and dishonest for a no-Status
  constructor; the lossless embed is the coherent choice.
- *"−0 is off the positive axis, so NaN is right."* Refuted by the
  C-semantics argument (`−0.0 < 0.0` is false) and by ln's own convention,
  which pfloat already followed; grouping −0 with the pole is what the
  committee and POSIX do.

## References

- pf-l5en, pf-jd4s, pf-k8ax, pf-qe5c (epic pf-8iji); the 2026-06-10
  convention decisions + verbatim citations (C23 N3088 F.10.9.5 /
  7.12.12.4-9, C11 N1570 F.10.3.7, POSIX Issue 8 y0/yn, RISC-V F ext).
- `src/cmp.rs` (`min`/`max`/`nan_number_reduce`/`min_signed_zero`/
  `max_signed_zero`/`limb_cmp_lsb_aligned`), `src/convert.rs` (`from_ieee`),
  `src/math/bessel_y.rs`, `src/math/bessel_k.rs`, `src/math/ci.rs`.
- Tests: `cmp::tests::{min_max_signed_zero_is_deterministic,
  min_max_signaling_nan_returns_the_number, total_cmp_equal_payload_across_precisions,
  min_max_nan_handling}`, `convert::tests::from_f32_f64_preserve_signaling_ness_losslessly`,
  `bessel_y::tests::{y0_negative_zero_is_the_pole_not_nan, y_negative_zero_is_the_pole}`,
  `bessel_k::tests::k_negative_zero_is_the_pole`.
