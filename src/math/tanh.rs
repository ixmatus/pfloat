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
//! Correctly rounded under every IEEE rounding mode via the shared
//! [`crate::math::ziv::ziv_round`] driver (slice p1.2, ADR-0022).
//! The composition `(1 − e^{−2|x|}) / (1 + e^{−2|x|})` runs at the
//! working precision the Ziv driver supplies; the internal `exp`
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

    let sign = x.sign();
    let abs_x = x.abs();
    ziv_round(
        |working_prec| tanh_at_w(&abs_x, sign, working_prec),
        target_precision,
        mode,
    )
}

/// Evaluate `tanh(x)` at the supplied working precision via
/// `tanh(|x|) = (1 − e^{−2|x|}) / (1 + e^{−2|x|})`, restoring the
/// sign of the input on the result. The caller's special-case
/// handling has already peeled off NaN, ±0, and ±∞. Returns the
/// unrounded value; the Ziv driver handles rounding to the
/// caller's target precision and mode.
fn tanh_at_w(abs_x: &BigFloat, sign: Sign, working_prec: u32) -> BigFloat {
    let abs_x_w = abs_x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let two = BigFloat::try_from_i64_exact(2, working_prec).expect("precision >= 1");
    let (two_x, _) = abs_x_w.mul(&two, RoundingMode::NearestEven);
    let neg_two_x = two_x.negated();
    let (exp_neg, _) = neg_two_x.exp(RoundingMode::NearestEven);
    let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let (numer, _) = one.sub(&exp_neg, RoundingMode::NearestEven);
    let (denom, _) = one.add(&exp_neg, RoundingMode::NearestEven);
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
}
