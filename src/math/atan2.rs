//! `atan2(y, x)`: two-argument arc tangent, returning the polar
//! angle of the point `(x, y)` in `(−π, π]`.
//!
//! Dispatch (IEEE 754-2019 §9.2.1):
//!
//! | y \ x      | `+0`     | `−0`     | `x > 0` finite | `x < 0` finite | `+∞`   | `−∞`   |
//! |------------|----------|----------|----------------|----------------|--------|--------|
//! | `+0`       | `+0`     | `+π`     | `+0`           | `+π`           | `+0`   | `+π`   |
//! | `−0`       | `−0`     | `−π`     | `−0`           | `−π`           | `−0`   | `−π`   |
//! | `y > 0`    | `+π/2`   | `+π/2`   | `atan(y/x)`    | `π − atan(...)`| `+0`   | `+π`   |
//! | `y < 0`    | `−π/2`   | `−π/2`   | `atan(y/x)`    | `−π + atan(...)`| `−0`  | `−π`   |
//! | `+∞`       | `+π/2`   | `+π/2`   | `+π/2`         | `+π/2`         | `+π/4` | `+3π/4`|
//! | `−∞`       | `−π/2`   | `−π/2`   | `−π/2`         | `−π/2`         | `−π/4` | `−3π/4`|
//!
//! `atan2(NaN, _) = atan2(_, NaN) = NaN`; signaling NaN raises
//! `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::{pi_at, pi_over_2_at};

impl BigFloat {
    /// `atan2(self, x)` returns the polar angle of `(x, self)`.
    /// Result rounded under `mode` to a precision of
    /// `max(self.precision, x.precision)`.
    #[must_use]
    pub fn atan2(&self, x: &Self, mode: RoundingMode) -> (Self, Status) {
        let target = self.precision.max(x.precision);
        self.atan2_round(x, target, mode)
            .expect("max of two valid precisions is valid")
    }

    /// `atan2(self, x)` with explicit result precision.
    pub fn atan2_round(
        &self,
        x: &Self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(atan2_kernel(self, x, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `atan2(self, x)` for `FixedFloat`. Delegates to
    /// [`BigFloat::atan2`].
    #[must_use]
    pub fn atan2(&self, x: &Self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().atan2(&x.to_big(), mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn atan2_kernel(
    y: &BigFloat,
    x: &BigFloat,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    // sNaN first: raise INVALID, return qNaN.
    if y.is_signaling_nan() || x.is_signaling_nan() {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }

    // qNaN propagation.
    if y.is_nan() {
        let nan =
            BigFloat::try_new_quiet_nan(y.sign(), target_precision, &[]).expect("precision >= 1");
        return (nan, Status::OK);
    }
    if x.is_nan() {
        let nan =
            BigFloat::try_new_quiet_nan(x.sign(), target_precision, &[]).expect("precision >= 1");
        return (nan, Status::OK);
    }

    let y_sign = y.sign();
    let x_sign = x.sign();

    // Both infinite: atan2(±∞, ±∞) ∈ {±π/4, ±3π/4}.
    if y.is_infinite() && x.is_infinite() {
        let pi = pi_at(target_precision);
        let one = BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1");
        let four = BigFloat::try_from_i64_exact(4, target_precision).expect("precision >= 1");
        let (pi_4, _) = pi.div(&four, RoundingMode::NearestEven);
        let result_abs = if matches!(x_sign, Sign::Positive) {
            pi_4
        } else {
            // 3π/4
            let three = BigFloat::try_from_i64_exact(3, target_precision).expect("precision >= 1");
            let (three_pi_4, _) = three.mul(&pi_4, RoundingMode::NearestEven);
            let _ = one;
            three_pi_4
        };
        let signed = if matches!(y_sign, Sign::Negative) {
            result_abs.negated()
        } else {
            result_abs
        };
        return (signed, Status::OK);
    }

    // y is ±∞, x finite: ±π/2.
    if y.is_infinite() {
        let pi_2 = pi_over_2_at(target_precision);
        let signed = if matches!(y_sign, Sign::Negative) {
            pi_2.negated()
        } else {
            pi_2
        };
        return (signed, Status::OK);
    }

    // x is ±∞, y finite (and not infinite from above).
    if x.is_infinite() {
        let result = if matches!(x_sign, Sign::Positive) {
            // +∞: angle = ±0 with sign of y.
            BigFloat::try_new_zero(y_sign, target_precision).expect("precision >= 1")
        } else {
            // −∞: angle = ±π with sign of y.
            let pi = pi_at(target_precision);
            if matches!(y_sign, Sign::Negative) {
                pi.negated()
            } else {
                pi
            }
        };
        return (result, Status::OK);
    }

    // y = ±0, x finite (and not infinite from above).
    if y.is_zero() {
        let result = match (y_sign, &x.class) {
            (
                Sign::Positive,
                Class::Zero {
                    sign: Sign::Positive,
                },
            ) => BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1"),
            (
                Sign::Negative,
                Class::Zero {
                    sign: Sign::Positive,
                },
            ) => BigFloat::try_new_zero(Sign::Negative, target_precision).expect("precision >= 1"),
            (
                Sign::Positive,
                Class::Zero {
                    sign: Sign::Negative,
                },
            ) => pi_at(target_precision),
            (
                Sign::Negative,
                Class::Zero {
                    sign: Sign::Negative,
                },
            ) => pi_at(target_precision).negated(),
            // x finite normal (sign-only dispatch).
            (s, _) if matches!(x_sign, Sign::Positive) => {
                BigFloat::try_new_zero(s, target_precision).expect("precision >= 1")
            }
            (s, _) => {
                let pi = pi_at(target_precision);
                if matches!(s, Sign::Negative) {
                    pi.negated()
                } else {
                    pi
                }
            }
        };
        return (result, Status::OK);
    }

    // x = ±0, y finite non-zero: ±π/2 with sign of y.
    if x.is_zero() {
        let pi_2 = pi_over_2_at(target_precision);
        let signed = if matches!(y_sign, Sign::Negative) {
            pi_2.negated()
        } else {
            pi_2
        };
        return (signed, Status::OK);
    }

    // Both finite and nonzero. Compute atan(y/x) and adjust by π
    // for the second/third quadrant.
    let working_prec = target_precision.saturating_add(64).min(1024);
    let y_w = y
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let x_w = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let (ratio, _) = y_w.div(&x_w, RoundingMode::NearestEven);
    let (at, _) = ratio.atan(RoundingMode::NearestEven);

    let result = if matches!(x_sign, Sign::Positive) {
        at
    } else {
        // x < 0: shift by ±π depending on sign of y.
        let pi = pi_at(working_prec);
        if matches!(y_sign, Sign::Positive) {
            let (shifted, _) = at.add(&pi, RoundingMode::NearestEven);
            shifted
        } else {
            let (shifted, _) = at.sub(&pi, RoundingMode::NearestEven);
            shifted
        }
    };

    let (rounded, status) = result
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    auto_raise(status);
    (rounded, status)
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
    fn atan2_quadrant_one() {
        // atan2(1, 1) = π/4
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.atan2(&one, RoundingMode::NearestEven);
        let pi_2 = super::super::pi_over_2_at(113);
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (pi_4, _) = pi_2.div(&two, RoundingMode::NearestEven);
        assert!(close_at(&r, &pi_4, 113 - 12));
    }

    #[test]
    fn atan2_quadrant_two() {
        // atan2(1, -1) = 3π/4
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let neg_one = BigFloat::try_from_i64_exact(-1, 113).unwrap();
        let (r, _) = one.atan2(&neg_one, RoundingMode::NearestEven);
        let pi = super::super::pi_at(113);
        let three = BigFloat::try_from_i64_exact(3, 113).unwrap();
        let four = BigFloat::try_from_i64_exact(4, 113).unwrap();
        let (three_pi, _) = three.mul(&pi, RoundingMode::NearestEven);
        let (three_pi_4, _) = three_pi.div(&four, RoundingMode::NearestEven);
        assert!(close_at(&r, &three_pi_4, 113 - 12));
    }

    #[test]
    fn atan2_quadrant_three() {
        // atan2(-1, -1) = -3π/4
        let neg_one = BigFloat::try_from_i64_exact(-1, 113).unwrap();
        let (r, _) = neg_one.atan2(&neg_one, RoundingMode::NearestEven);
        let pi = super::super::pi_at(113);
        let three = BigFloat::try_from_i64_exact(3, 113).unwrap();
        let four = BigFloat::try_from_i64_exact(4, 113).unwrap();
        let (three_pi, _) = three.mul(&pi, RoundingMode::NearestEven);
        let (three_pi_4, _) = three_pi.div(&four, RoundingMode::NearestEven);
        let expected = three_pi_4.negated();
        assert!(close_at(&r, &expected, 113 - 12));
    }

    #[test]
    fn atan2_quadrant_four() {
        // atan2(-1, 1) = -π/4
        let neg_one = BigFloat::try_from_i64_exact(-1, 113).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = neg_one.atan2(&one, RoundingMode::NearestEven);
        let pi_2 = super::super::pi_over_2_at(113);
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (pi_4, _) = pi_2.div(&two, RoundingMode::NearestEven);
        let expected = pi_4.negated();
        assert!(close_at(&r, &expected, 113 - 12));
    }

    #[test]
    fn atan2_pos_y_zero_x_is_pi_over_2() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let z = BigFloat::try_new_zero(Sign::Positive, 113).unwrap();
        let (r, _) = one.atan2(&z, RoundingMode::NearestEven);
        let pi_2 = super::super::pi_over_2_at(113);
        assert!(close_at(&r, &pi_2, 100));
    }

    #[test]
    fn atan2_zero_zero() {
        // atan2(+0, +0) = +0
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, _) = z.atan2(&z, RoundingMode::NearestEven);
        assert!(r.is_zero() && r.is_sign_positive());

        // atan2(-0, +0) = -0
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, _) = nz.atan2(&z, RoundingMode::NearestEven);
        assert!(r.is_zero() && r.is_sign_negative());

        // atan2(+0, -0) = +π
        let (r, _) = z.atan2(&nz, RoundingMode::NearestEven);
        let pi = super::super::pi_at(53);
        assert_eq!(r.partial_cmp(&pi).0, Some(Ordering::Equal));

        // atan2(-0, -0) = -π
        let (r, _) = nz.atan2(&nz, RoundingMode::NearestEven);
        let neg_pi = super::super::pi_at(53).negated();
        assert_eq!(r.partial_cmp(&neg_pi).0, Some(Ordering::Equal));
    }

    #[test]
    fn atan2_pos_inf_pos_inf_is_pi_over_4() {
        let pi_inf = BigFloat::try_new_infinity(Sign::Positive, 113).unwrap();
        let (r, _) = pi_inf.atan2(&pi_inf, RoundingMode::NearestEven);
        let pi = super::super::pi_at(113);
        let four = BigFloat::try_from_i64_exact(4, 113).unwrap();
        let (pi_4, _) = pi.div(&four, RoundingMode::NearestEven);
        assert!(close_at(&r, &pi_4, 113 - 12));
    }

    #[test]
    fn atan2_pos_y_neg_inf_x_is_pi() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let neg_inf = BigFloat::try_new_infinity(Sign::Negative, 113).unwrap();
        let (r, _) = one.atan2(&neg_inf, RoundingMode::NearestEven);
        let pi = super::super::pi_at(113);
        assert!(close_at(&r, &pi, 100));
    }

    #[test]
    fn atan2_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, _) = q.atan2(&one, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        let (r, _) = one.atan2(&q, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn atan2_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, status) = sn.atan2(&one, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
