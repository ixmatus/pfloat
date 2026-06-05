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
//!
//! Correctly rounded under every IEEE 754-2019 rounding mode via
//! the shared [`crate::math::ziv::ziv_round`] driver on the
//! finite/nonzero finite path (slice p1.25, ADR-0038). Every
//! special case that returns an irrational constant (`±π`,
//! `±π/2`, `±π/4`, `±3π/4`) rounds via the mode-aware helpers
//! [`super::pi_at_round`] / [`super::pi_over_2_at_round`] (or the
//! local boost-then-round pattern for the `π/4` / `3π/4` cases).

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::ziv::ziv_round;
use super::ziv_calibration::ATAN2_ERROR_GUARD;
use super::{pi_at, pi_at_round, pi_over_2_at_round};

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

    // Both infinite: atan2(±∞, ±∞) ∈ {±π/4, ±3π/4}. Compute at
    // boosted precision and round to target under mode (slice
    // p1.25 — directed-mode awareness for the irrational-constant
    // returns).
    if y.is_infinite() && x.is_infinite() {
        let boost = target_precision.saturating_add(128);
        let pi = pi_at(boost);
        let four = BigFloat::try_from_i64_exact(4, boost).expect("precision >= 1");
        let (pi_4, _) = pi.div(&four, RoundingMode::NearestEven);
        let abs_unrounded = if matches!(x_sign, Sign::Positive) {
            pi_4
        } else {
            // 3π/4
            let three = BigFloat::try_from_i64_exact(3, boost).expect("precision >= 1");
            three.mul(&pi_4, RoundingMode::NearestEven).0
        };
        let signed_unrounded = if matches!(y_sign, Sign::Negative) {
            abs_unrounded.negated()
        } else {
            abs_unrounded
        };
        let (rounded, status) = signed_unrounded
            .round_to_precision(target_precision, mode)
            .expect("target precision >= 1");
        auto_raise(status);
        return (rounded, status);
    }

    // y is ±∞, x finite: ±π/2.
    if y.is_infinite() {
        let (pi_2, status) = pi_over_2_at_round(target_precision, mode);
        let signed = if matches!(y_sign, Sign::Negative) {
            pi_2.negated()
        } else {
            pi_2
        };
        auto_raise(status);
        return (signed, status);
    }

    // x is ±∞, y finite (and not infinite from above).
    if x.is_infinite() {
        let result = if matches!(x_sign, Sign::Positive) {
            // +∞: angle = ±0 with sign of y. Exact; no mode needed.
            (
                BigFloat::try_new_zero(y_sign, target_precision).expect("precision >= 1"),
                Status::OK,
            )
        } else {
            // −∞: angle = ±π with sign of y. Mode-aware.
            let (pi, status) = pi_at_round(target_precision, mode);
            if matches!(y_sign, Sign::Negative) {
                (pi.negated(), status)
            } else {
                (pi, status)
            }
        };
        auto_raise(result.1);
        return result;
    }

    // y = ±0, x finite (and not infinite from above).
    if y.is_zero() {
        let (result, status) = match (y_sign, &x.class) {
            (
                Sign::Positive,
                Class::Zero {
                    sign: Sign::Positive,
                },
            ) => (
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1"),
                Status::OK,
            ),
            (
                Sign::Negative,
                Class::Zero {
                    sign: Sign::Positive,
                },
            ) => (
                BigFloat::try_new_zero(Sign::Negative, target_precision).expect("precision >= 1"),
                Status::OK,
            ),
            (
                Sign::Positive,
                Class::Zero {
                    sign: Sign::Negative,
                },
            ) => pi_at_round(target_precision, mode),
            (
                Sign::Negative,
                Class::Zero {
                    sign: Sign::Negative,
                },
            ) => {
                let (pi, status) = pi_at_round(target_precision, mode);
                (pi.negated(), status)
            }
            // x finite normal (sign-only dispatch).
            (s, _) if matches!(x_sign, Sign::Positive) => (
                BigFloat::try_new_zero(s, target_precision).expect("precision >= 1"),
                Status::OK,
            ),
            (s, _) => {
                let (pi, status) = pi_at_round(target_precision, mode);
                if matches!(s, Sign::Negative) {
                    (pi.negated(), status)
                } else {
                    (pi, status)
                }
            }
        };
        auto_raise(status);
        return (result, status);
    }

    // x = ±0, y finite non-zero: ±π/2 with sign of y.
    if x.is_zero() {
        let (pi_2, status) = pi_over_2_at_round(target_precision, mode);
        let signed = if matches!(y_sign, Sign::Negative) {
            pi_2.negated()
        } else {
            pi_2
        };
        auto_raise(status);
        return (signed, status);
    }

    // Both finite and nonzero. Ziv-driven correct rounding under
    // every IEEE mode. The eval closure captures y and x and runs
    // the existing finite-case composition (`atan(y/x)` + quadrant
    // shift) at working precision `w` under NE. The quadrant shift
    // for x < 0 uses pi_at(w) so the working-precision π scales
    // with the Ziv guard.
    let (result, status) = ziv_round(
        |w| {
            let y_w = y
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let x_w = x
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let (ratio, _) = y_w.div(&x_w, RoundingMode::NearestEven);
            let (at, _) = ratio.atan(RoundingMode::NearestEven);

            if matches!(x_sign, Sign::Positive) {
                at
            } else {
                let pi = pi_at(w);
                if matches!(y_sign, Sign::Positive) {
                    let (shifted, _) = at.add(&pi, RoundingMode::NearestEven);
                    shifted
                } else {
                    let (shifted, _) = at.sub(&pi, RoundingMode::NearestEven);
                    shifted
                }
            }
        },
        target_precision,
        mode,
        ATAN2_ERROR_GUARD,
    );
    // atan2(y, x) for finite nonzero y and x off the axes is
    // transcendental (Lindemann–Weierstrass), hence irrational, hence
    // INEXACT even where it rounds onto a grid value (pf-uqd1, ADR-0063).
    // The only exact result, atan2(+0, x > 0) = +0, and the axis cases
    // (±π/2, ±π, ±π/4, ±3π/4 via pi_*_at_round, already INEXACT) are
    // dispatched above.
    let status = if matches!(result.class, Class::Normal { .. }) {
        status | Status::INEXACT
    } else {
        status
    };
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
