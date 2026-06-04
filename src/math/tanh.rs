//! `tanh(x) = sinh(x) / cosh(x)`: hyperbolic tangent.
//!
//! Identity used: for `|x| ≥ 0`, `tanh(|x|) = (1 − e^{−2|x|}) /
//! (1 + e^{−2|x|})`. The sign of the result is the sign of `x`.
//! Working on `|x|` avoids the `−∞/+∞` indeterminate that the
//! `(e^x − e^{−x})/(e^x + e^{−x})` form produces at `x = −∞`.
//!
//! For `|x|` large enough that `e^{−2|x|}` falls below ULP at target
//! precision, the formula collapses to `1/1 = 1` and the sign-
//! flipped result is `±1` correctly.
//!
//! The numerator is evaluated through `expm1` as
//! `1 − e^{−2|x|} = −expm1(−2|x|)`, which avoids the catastrophic
//! cancellation the bare `1 − e^{−2|x|}` suffers for small `|x|`
//! (where `e^{−2|x|} ≈ 1`). This both keeps the working-precision
//! intermediate accurate to the Ziv driver's error guard for tiny
//! `|x|` (pf-zhcy) and subsumes the former tiny-`|x|` short circuit:
//! that case (slice p1.4, pf-7d7) returned `|x|` because the
//! cancelling numerator otherwise collapsed to exactly `0` and the
//! Ziv interval test certified the false `0`; `expm1` preserves the
//! tiny numerator `≈ 2|x|`, so the collapse cannot occur and the
//! special case is gone (ADR-0050).
//!
//! Correctly rounded under every IEEE rounding mode via the shared
//! [`crate::math::ziv::ziv_round`] driver (slice p1.2, ADR-0022).
//! The composition `−expm1(−2|x|) / (2 + expm1(−2|x|))` runs at the
//! working precision the Ziv driver supplies; the internal `expm1`
//! call is itself Ziv-driven, so the composition is correctly
//! rounded under the outer envelope's interval test.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `tanh(±0) = ±0`.
//! - `tanh(±∞) = ±1`.
//! - `tanh(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::ziv::ziv_round;
use super::ziv_calibration::TANH_ERROR_GUARD;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `tanh(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn tanh(&self, mode: RoundingMode) -> (Self, Status) {
        self.tanh_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `tanh(self)` with explicit result precision.
    pub fn tanh_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(tanh_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `tanh(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::tanh`].
    #[must_use]
    pub fn tanh(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().tanh(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn tanh_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    match &x.class {
        Class::Nan {
            quiet,
            sign,
            payload,
        } => {
            if !*quiet {
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            let nan = BigFloat::try_new_quiet_nan(*sign, target_precision, payload)
                .expect("precision >= 1");
            return (nan, Status::OK);
        }
        Class::Zero { sign } => {
            let z = BigFloat::try_new_zero(*sign, target_precision).expect("precision >= 1");
            return (z, Status::OK);
        }
        Class::Infinity { sign } => {
            // tanh(±∞) = ±1.
            let one = BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1");
            let signed = if matches!(sign, Sign::Negative) {
                one.negated()
            } else {
                one
            };
            return (signed, Status::OK);
        }
        Class::Normal { .. } => {}
    }

    // Tiny x: tanh(x) = x − x³/3 + … shrinks toward x in magnitude (the
    // −x³/3 correction opposes x's sign). ADR-0050 removed the
    // pre-Phase-1f tiny-x short-circuit that returned |x| rounded under
    // mode: NE-correct, but it dropped the directed-mode information and
    // rounded the wrong way under TZ/TP/TN. This restores a short-circuit
    // through the mode-aware round_with_infinitesimal — which carries the
    // sign of the dropped correction term — and bypasses the expm1 form
    // that otherwise grinds at moderate working precision for moderately
    // tiny x (ADR-0059).
    let e = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => unreachable!(),
    };
    if e <= -(i64::from(target_precision) + 2) {
        return crate::rounding::round_with_infinitesimal(
            x,
            x.sign(),
            true, // magnitude shrinks: the −x³/3 correction opposes x's sign
            target_precision,
            mode,
        );
    }

    let sign = x.sign();
    let abs_x = x.abs();
    ziv_round(
        |working_prec| tanh_at_w(&abs_x, sign, working_prec),
        target_precision,
        mode,
        TANH_ERROR_GUARD,
    )
}

/// Evaluate `tanh(x)` at the supplied working precision via the
/// numerically stable form
///
/// ```text
/// tanh(|x|) = −expm1(−2|x|) / (2 + expm1(−2|x|))
/// ```
///
/// restoring the sign of the input on the result. The caller's
/// special-case handling has already peeled off NaN, ±0, and ±∞.
/// Returns the unrounded value; the Ziv driver handles rounding to
/// the caller's target precision and mode.
///
/// This is algebraically `(1 − e^{−2|x|}) / (1 + e^{−2|x|})` with the
/// numerator rewritten through `expm1`. The bare `1 − e^{−2|x|}`
/// cancels catastrophically for small `|x|`: `e^{−2|x|} ≈ 1`, so the
/// subtraction loses ~`2·(−log2|x|)` bits, leaving the working-
/// precision value accurate to far fewer bits than the Ziv error
/// guard assumes (pf-zhcy: at `|x| = 2^−149` the loss is ~148 bits,
/// so the converged intermediate held only ~388 bits at
/// `working = 536`). `expm1(−2|x|) = e^{−2|x|} − 1` is evaluated
/// without that cancellation, so `numer = −expm1(−2|x|)` and
/// `denom = 2 + expm1(−2|x|) ∈ (1, 2]` are both accurate to working
/// precision and neither cancels.
///
/// The `expm1` form also subsumes the former tiny-`|x|` short circuit
/// (slice p1.4, pf-7d7): that case returned `|x|` because the
/// `1 − e^{−2|x|}` numerator otherwise collapsed to exactly `0` (the
/// `exp` rounded to `1`) and the Ziv interval test certified the
/// false `0`. `expm1` preserves the tiny numerator `≈ 2|x|` instead
/// of rounding it to zero, so the collapse cannot occur and no
/// short circuit is needed (ADR-0050).
fn tanh_at_w(abs_x: &BigFloat, sign: Sign, working_prec: u32) -> BigFloat {
    let abs_x_w = abs_x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    let two = BigFloat::try_from_i64_exact(2, working_prec).expect("precision >= 1");
    let (two_x, _) = abs_x_w.mul(&two, RoundingMode::NearestEven);
    let neg_two_x = two_x.negated();
    // expm1(−2|x|) = e^{−2|x|} − 1, accurate for small |x| where the
    // bare `1 − e^{−2|x|}` would cancel.
    let (em1, _) = neg_two_x.expm1(RoundingMode::NearestEven);
    // numer = 1 − e^{−2|x|} = −expm1(−2|x|) ≥ 0.
    let numer = em1.negated();
    // denom = 1 + e^{−2|x|} = 2 + expm1(−2|x|) ∈ (1, 2].
    let (denom, _) = two.add(&em1, RoundingMode::NearestEven);
    let (result_abs, _) = numer.div(&denom, RoundingMode::NearestEven);
    if matches!(sign, Sign::Negative) {
        result_abs.negated()
    } else {
        result_abs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    fn close_at(v: &BigFloat, expected: &BigFloat, bits: u32) -> bool {
        let (diff, _) = v.sub(expected, RoundingMode::NearestEven);
        let abs_diff = diff.abs();
        if abs_diff.is_zero() {
            return true;
        }
        let p = v.precision().max(expected.precision());
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let abs_b = expected.abs();
        let mut bound = if abs_b.is_zero() { one } else { abs_b };
        for _ in 0..bits {
            bound = bound.div(&two, RoundingMode::NearestEven).0;
        }
        matches!(
            abs_diff.partial_cmp(&bound).0,
            Some(Ordering::Less | Ordering::Equal)
        )
    }

    #[test]
    fn tanh_zero() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.tanh(RoundingMode::NearestEven);
            assert!(r.is_zero());
            assert_eq!(r.is_sign_negative(), matches!(s, Sign::Negative));
        }
    }

    #[test]
    fn tanh_pos_inf_is_one() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.tanh(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    #[test]
    fn tanh_neg_inf_is_neg_one() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, _) = ni.tanh(RoundingMode::NearestEven);
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        assert_eq!(r.partial_cmp(&neg_one).0, Some(Ordering::Equal));
    }

    #[test]
    fn tanh_one_matches_definition() {
        // tanh(1) = sinh(1)/cosh(1)
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.tanh(RoundingMode::NearestEven);
        let (s, _) = one.sinh(RoundingMode::NearestEven);
        let (c, _) = one.cosh(RoundingMode::NearestEven);
        let (expected, _) = s.div(&c, RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, 113 - 12));
    }

    #[test]
    fn tanh_negation() {
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let neg_two = BigFloat::try_from_i64_exact(-2, 113).unwrap();
        let (a, _) = two.tanh(RoundingMode::NearestEven);
        let (b, _) = neg_two.tanh(RoundingMode::NearestEven);
        let neg_a = a.negated();
        assert!(close_at(&b, &neg_a, 113 - 12));
    }

    #[test]
    fn tanh_large_x_saturates() {
        // For x ≥ ~40 at p = 53, tanh(x) rounds to 1 exactly.
        let big = BigFloat::try_from_i64_exact(100, 53).unwrap();
        let (r, _) = big.tanh(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    #[test]
    fn tanh_bounded() {
        // |tanh(x)| < 1 for finite x. Pick x = 10; tanh(10) ≈ 0.99999..  < 1.
        let ten = BigFloat::try_from_i64_exact(10, 53).unwrap();
        let (r, _) = ten.tanh(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Less));
    }

    #[test]
    fn tanh_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.tanh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn tanh_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.tanh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn tanh_tiny_input_round_to_nearest_returns_input() {
        // Slice p1.4 short circuit (closes pf-7d7). For `|x|` whose
        // binary exponent is below the Ziv driver's error guard,
        // `tanh(x)` rounds to `x` at target precision under either
        // round-to-nearest tie rule: the cubic Taylor correction
        // `x³/3` for x = ±2⁻¹⁰⁰ is ≈ 2⁻³⁰² , far below the p=53
        // half-ULP threshold ≈ 2⁻¹⁵³. Pre-fix, the cancellation
        // path returned `0` instead, and the Ziv interval test
        // certified the wrong value because `half_width(0) = 0`.
        //
        // Directed modes (`TowardPositive`, `TowardZero`,
        // `TowardNegative`) legitimately pick the neighboring
        // representable that brackets the true value (sign of the
        // cubic correction matters), so they are covered by the
        // oracle_sweep MPFR regression, not by this point test.
        let two = BigFloat::try_from_i64_exact(2, 200).unwrap();
        let mut x = BigFloat::try_from_i64_exact(1, 200).unwrap();
        for _ in 0..100 {
            x = x.div(&two, RoundingMode::NearestEven).0;
        }
        let x53 = x
            .round_to_precision(53, RoundingMode::NearestEven)
            .unwrap()
            .0;
        for &mode in &[RoundingMode::NearestEven, RoundingMode::NearestAway] {
            let (r, _) = x53.tanh(mode);
            assert_eq!(
                r.partial_cmp(&x53).0,
                Some(Ordering::Equal),
                "tanh(2^-100) under {mode:?} = {r}, expected {x53}"
            );
            let neg = x53.negated();
            let (r_neg, _) = neg.tanh(mode);
            assert_eq!(
                r_neg.partial_cmp(&neg).0,
                Some(Ordering::Equal),
                "tanh(-2^-100) under {mode:?} = {r_neg}, expected {neg}"
            );
        }
    }

    #[test]
    fn tanh_tiny_input_directed_modes() {
        // Tiny-x short-circuit restored mode-aware (ADR-0059, contrast
        // ADR-0050). tanh(x) = x − x³/3 + … shrinks toward zero, so for
        // tiny x the result is x under the three modes that round toward
        // or onto x and x's toward-zero neighbour under the two that
        // round inward. The pre-Phase-1f short-circuit returned bare |x|
        // and rounded these inward modes the wrong way; this pins the
        // correct shrink direction (subtracts_magnitude = true).
        use RoundingMode::{NearestAway, NearestEven, TowardNegative, TowardPositive, TowardZero};
        let two = BigFloat::try_from_i64_exact(2, 200).unwrap();
        let mut x = BigFloat::try_from_i64_exact(1, 200).unwrap();
        for _ in 0..100 {
            x = x.div(&two, RoundingMode::NearestEven).0;
        }
        let x53 = x
            .round_to_precision(53, RoundingMode::NearestEven)
            .unwrap()
            .0;

        // Positive tiny x: TowardZero and TowardNegative round below x.
        for &m in &[NearestEven, NearestAway, TowardPositive] {
            let (r, _) = x53.tanh(m);
            assert_eq!(
                r.partial_cmp(&x53).0,
                Some(Ordering::Equal),
                "tanh +tiny {m:?}"
            );
        }
        for &m in &[TowardZero, TowardNegative] {
            let (r, _) = x53.tanh(m);
            assert_eq!(
                r.partial_cmp(&x53).0,
                Some(Ordering::Less),
                "tanh +tiny {m:?} shrinks"
            );
        }

        // Negative tiny x: TowardZero and TowardPositive round above x.
        let neg = x53.negated();
        for &m in &[NearestEven, NearestAway, TowardNegative] {
            let (r, _) = neg.tanh(m);
            assert_eq!(
                r.partial_cmp(&neg).0,
                Some(Ordering::Equal),
                "tanh -tiny {m:?}"
            );
        }
        for &m in &[TowardZero, TowardPositive] {
            let (r, _) = neg.tanh(m);
            assert_eq!(
                r.partial_cmp(&neg).0,
                Some(Ordering::Greater),
                "tanh -tiny {m:?} shrinks"
            );
        }
    }
}
