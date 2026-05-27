//! `cosh(x) = (e^x + e^−x) / 2`: hyperbolic cosine.
//!
//! Direct evaluation: `(exp(x) + exp(−x)) / 2`. Both summands are
//! non-negative, so the addition has no cancellation regardless of
//! the sign or magnitude of `x`. For `|x| → ∞`, one summand
//! dominates and the result is `+∞`.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `cosh(±0) = 1`.
//! - `cosh(±∞) = +∞`.
//! - `cosh(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::ziv::ziv_round;
use super::ziv_calibration::COSH_ERROR_GUARD;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `cosh(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn cosh(&self, mode: RoundingMode) -> (Self, Status) {
        self.cosh_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `cosh(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.27, ADR-0038).
    pub fn cosh_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(cosh_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `cosh(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::cosh`].
    #[must_use]
    pub fn cosh(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().cosh(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn cosh_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
        Class::Zero { .. } => {
            return (
                BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Infinity { .. } => {
            return (
                BigFloat::try_new_infinity(Sign::Positive, target_precision)
                    .expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Normal { .. } => {}
    }

    // Ziv-driven correct rounding under every IEEE mode. The eval
    // closure runs the cancellation-free `(exp(x) + exp(-x))/2`
    // composition at working precision `w`. exp is Ziv-driven
    // (slice p1.2) and both summands are positive, so no
    // cancellation regime.
    let (result, status) = ziv_round(
        |w| {
            let x_w = x
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let neg_x = x_w.negated();
            let (e_pos, _) = x_w.exp(RoundingMode::NearestEven);
            let (e_neg, _) = neg_x.exp(RoundingMode::NearestEven);
            let (sum, _) = e_pos.add(&e_neg, RoundingMode::NearestEven);
            let two = BigFloat::try_from_i64_exact(2, w).expect("precision >= 1");
            sum.div(&two, RoundingMode::NearestEven).0
        },
        target_precision,
        mode,
        COSH_ERROR_GUARD,
    );
    auto_raise(status);
    (result, status)
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
    fn cosh_zero_is_one() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.cosh(RoundingMode::NearestEven);
            let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
            assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
        }
    }

    #[test]
    fn cosh_inf_is_pos_inf() {
        for s in [Sign::Positive, Sign::Negative] {
            let i = BigFloat::try_new_infinity(s, 53).unwrap();
            let (r, _) = i.cosh(RoundingMode::NearestEven);
            assert!(r.is_infinite() && r.is_sign_positive());
        }
    }

    #[test]
    fn cosh_is_even() {
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let neg_two = BigFloat::try_from_i64_exact(-2, 113).unwrap();
        let (a, _) = two.cosh(RoundingMode::NearestEven);
        let (b, _) = neg_two.cosh(RoundingMode::NearestEven);
        assert!(close_at(&a, &b, 113 - 12));
    }

    #[test]
    fn cosh_one_matches_definition() {
        // cosh(1) = (e + 1/e)/2 ≈ 1.5430806348152437
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.cosh(RoundingMode::NearestEven);
        let (e, _) = one.exp(RoundingMode::NearestEven);
        let neg_one = BigFloat::try_from_i64_exact(-1, 113).unwrap();
        let (e_neg, _) = neg_one.exp(RoundingMode::NearestEven);
        let (sum, _) = e.add(&e_neg, RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (expected, _) = sum.div(&two, RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, 113 - 12));
    }

    #[test]
    fn cosh_pythagorean_identity() {
        // cosh²(x) − sinh²(x) = 1
        let x = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (c, _) = x.cosh(RoundingMode::NearestEven);
        let (s, _) = x.sinh(RoundingMode::NearestEven);
        let (c2, _) = c.mul(&c, RoundingMode::NearestEven);
        let (s2, _) = s.mul(&s, RoundingMode::NearestEven);
        let (diff, _) = c2.sub(&s2, RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        assert!(close_at(&diff, &one, 113 - 16));
    }

    #[test]
    fn cosh_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.cosh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn cosh_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.cosh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
