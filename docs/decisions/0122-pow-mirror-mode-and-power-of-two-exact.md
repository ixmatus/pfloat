# ADR-0122: pow — mirror the mode under odd-parity negation, and dispatch the power-of-two exact set

- **Status**: accepted
- **Date**: 2026-07-03

## Context

Two `pow` contract defects from the 2026-06 review (epic pf-8iji).

1. **pf-l38k — odd-parity negation drops the directed mode.** For a
   negative base and integer exponent, `pow` computes `|x|^y` under
   `mode` and negates it for odd parity. Directed rounding is relative to
   the true (negative) value, so rounding `|x|^y` under `mode` and then
   flipping the sign lands on the wrong side: `(-3)^7 = -2187` at target 8
   under TowardPositive rounded `2187` UP to `2192` then negated to
   `-2192`, where `-2176` is due (TowardPositive on a negative value
   rounds toward less-negative). The pf-egxm/pf-yprp negate-after-directed-
   round family.

2. **pf-vzim — an exact result over-reports INEXACT.** `pow(4, 1.5) =
   2^(2·1.5) = 2^3 = 8` is exactly representable, but the general path
   `exp(y·ln x)` computes `8 ± ε` and the Ziv driver reports INEXACT — the
   exact-value-defeats-Ziv family (ADR-0039/0063), here for an algebraic
   (not transcendental) exact case.

## Decision

1. pf-l38k: in the negative-base odd-parity branch, round `|x|^y` under
   `mirror_mode_for_negation(mode)` before negating (the existing helper
   used by `signed_constant_at_round`). Even parity is unchanged.

2. pf-vzim: before the `exp·ln` Ziv path, dispatch the **power-of-two
   base** exact subset (`try_pow_pow2_base_exact`): if `x = 2^e` (a single
   significant bit — `round_to_precision(x, 1)` is exact) and `e·y` is an
   exact integer in the representable range (`integer_exponent(e_bf · y)`
   at a width sized to hold the exact dyadic product), then `x^y =
   2^(e·y)` is a single-bit value, exactly representable — return it with
   OK. A representable-range overflow falls through to the Ziv path (which
   produces the honest ±∞/±0 + OVERFLOW/UNDERFLOW). Non-power-of-two bases
   and non-integer `e·y` fall through unchanged, so `4^1.6` (irrational)
   stays correctly-rounded INEXACT.

   The **general** perfect-power exact set — a base `x = m^j` with rational
   `y` (`9^0.5 = 3`, `27^(1/3) = 3`) — needs integer-root reasoning over
   the base and is out of scope for this pass; it is the pf-vzim Fable
   escalation the review flagged, and is left filed. Only the dyadic
   power-of-two subset (which covers the reported `4^1.5`) is dispatched.

## Consequences

- Negative-base powers round correctly under every mode; power-of-two
  exact powers return OK. The common path (positive base, general
  exponent; non-power-of-two base) is unchanged — the exact dispatch is
  cheap (one round-to-1-bit, one multiply, one integer check) and only
  fires for a power-of-two base with integer `e·y`.
- Verified against MPFR: `differential_pow` unchanged.

### Inversion (failure paragraphs considered)

- *"Trust the Ziv INEXACT flag for pf-vzim."* The flag is honest about the
  *computation* (the `exp·ln` chain rounded) but wrong about the *IEEE
  result* (`8` is exact); the exact-value cases need a pre-Ziv dispatch,
  the same shape as `gamma`'s integer walk (ADR-0039).
- *"Handle the whole perfect-power exact set now."* Deferred: detecting `x
  = m^j` for a general integer `m` and reconciling with a rational `y`
  needs integer-root machinery and careful exact-set reasoning (the review's
  Fable escalation); the power-of-two subset is the sound, self-contained
  slice that covers the reported case, and over-firing it would wrongly
  flag an inexact result OK — so it is gated to provable exactness.
- *"Mirror is a no-op for TowardZero/nearest."* Correct, and intended:
  `mirror_mode_for_negation` is identity for TZ/NE/NA and swaps TP↔TN, so
  only the sign-asymmetric directed modes change — exactly the ones the
  bug mis-rounded.

## References

- pf-l38k, pf-vzim (epic pf-8iji); reproducers P1, Q6 in the review
  harness (`(-3)^7 @8 TP → -2192`; `4^1.5 → 8` INEXACT).
- `mirror_mode_for_negation` / `signed_constant_at_round` (src/math/mod.rs),
  ADR-0039/0063 (the exact-value-defeats-Ziv family), the pf-egxm/pf-yprp
  negate-after-directed-round siblings.
- `src/math/pow.rs` (`pow_kernel`, `pow_positive`, `try_pow_pow2_base_exact`);
  `src/math/pow.rs` tests `pow_odd_negative_base_mirrors_directed_mode`,
  `pow_power_of_two_base_exact_is_ok`.
