//! `exp2(x) = 2^x`: binary exponential.
//!
//! Composition: `2^x = exp(x · ln(2))`. The kernel computes the
//! product at working precision (with `ln(2)` from the shared
//! 1024-bit constant), then calls `exp`. All special cases flow
//! through composition.
//!
//! Correctly rounded under every IEEE 754-2019 rounding mode via
//! the shared [`crate::math::ziv::ziv_round`] driver (slice p1.24,
//! ADR-0038). The `exp · ln(2)` composition has no cancellation
//! regime; the Ziv envelope's working-precision growth certifies
//! the rounding-mode interval test on the final round.
//!
//! Special cases per IEEE 754-2019 §9.2 reduce to:
//!
//! - `exp2(±0) = 1`.
//! - `exp2(+∞) = +∞`, `exp2(−∞) = +0`.
//! - `exp2(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::ln_2_at;
use super::ziv::ziv_round;

impl BigFloat {
    /// `2^self` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn exp2(&self, mode: RoundingMode) -> (Self, Status) {
        self.exp2_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `exp2(self)` with explicit result precision.
    pub fn exp2_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(exp2_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `exp2(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::exp2`].
    #[must_use]
    pub fn exp2(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().exp2(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn exp2_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
        Class::Infinity {
            sign: Sign::Positive,
        } => {
            return (
                BigFloat::try_new_infinity(Sign::Positive, target_precision)
                    .expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Infinity {
            sign: Sign::Negative,
        } => {
            return (
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Normal { .. } => {}
    }

    // Ziv-driven correct rounding under every IEEE mode. The
    // composition `exp(x · ln(2))` has no cancellation regime; the
    // Ziv driver's working-precision growth handles the rounding-
    // mode interval test at the final round to target.
    let (result, status) = ziv_round(
        |w| {
            let x_w = x
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let ln_2 = ln_2_at(w);
            let (product, _) = x_w.mul(&ln_2, RoundingMode::NearestEven);
            let (e_val, _) = product.exp(RoundingMode::NearestEven);
            e_val
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0
        },
        target_precision,
        mode,
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
    fn exp2_zero_is_one() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, _) = z.exp2(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    #[test]
    fn exp2_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.exp2(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn exp2_neg_inf_is_zero() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, _) = ni.exp2(RoundingMode::NearestEven);
        assert!(r.is_zero() && r.is_sign_positive());
    }

    #[test]
    fn exp2_one_is_two() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.exp2(RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        assert!(close_at(&r, &two, 113 - 8));
    }

    #[test]
    fn exp2_ten_is_1024() {
        let ten = BigFloat::try_from_i64_exact(10, 113).unwrap();
        let (r, _) = ten.exp2(RoundingMode::NearestEven);
        let k = BigFloat::try_from_i64_exact(1024, 113).unwrap();
        assert!(close_at(&r, &k, 113 - 12));
    }

    #[test]
    fn exp2_negative_ten() {
        // 2^-10 = 1/1024
        let neg_ten = BigFloat::try_from_i64_exact(-10, 113).unwrap();
        let (r, _) = neg_ten.exp2(RoundingMode::NearestEven);
        let recip = BigFloat::parse_str("0.0009765625", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        assert!(close_at(&r, &recip, 113 - 12));
    }

    #[test]
    fn exp2_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.exp2(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn exp2_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.exp2(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixed_exp2() {
        let one = FixedFloat::<53>::try_from_i64_exact(1).unwrap();
        let (r, _) = one.exp2(RoundingMode::NearestEven);
        let two = FixedFloat::<53>::try_from_i64_exact(2).unwrap();
        assert_eq!(r.partial_cmp(&two).0, Some(Ordering::Equal));
    }
}
