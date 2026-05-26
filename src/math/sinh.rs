//! `sinh(x) = (e^x − e^−x) / 2`: hyperbolic sine.
//!
//! Composition: `sinh(x) = (expm1(x) − expm1(−x)) / 2`. The
//! `expm1`-based form avoids the cancellation of `exp(x) − exp(−x)`
//! for `|x| < 1`. `expm1` for large `|x|` returns the same magnitude
//! as `exp` (the `−1` is below ULP), so the subtraction at large
//! `|x|` is dominated by the larger term.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `sinh(±0) = ±0`.
//! - `sinh(±∞) = ±∞`.
//! - `sinh(NaN) = NaN`; `sNaN` raises `INVALID`.

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
    /// `sinh(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn sinh(&self, mode: RoundingMode) -> (Self, Status) {
        self.sinh_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `sinh(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.27, ADR-0038).
    pub fn sinh_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(sinh_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `sinh(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::sinh`].
    #[must_use]
    pub fn sinh(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().sinh(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn sinh_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            return (
                BigFloat::try_new_infinity(*sign, target_precision).expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Normal { .. } => {}
    }

    // Ziv-driven correct rounding under every IEEE mode. The eval
    // closure runs the existing `(expm1(x) - expm1(-x))/2`
    // composition at working precision `w` under NE; the outer
    // envelope certifies the rounding-mode interval test on the
    // final round. expm1 is itself Ziv-driven (slice p1.24) and
    // its cancellation boost handles the small-|x| regime.
    let (result, status) = ziv_round(
        |w| {
            let x_w = x
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let neg_x = x_w.negated();
            let (em1_pos, _) = x_w.expm1(RoundingMode::NearestEven);
            let (em1_neg, _) = neg_x.expm1(RoundingMode::NearestEven);
            let (diff, _) = em1_pos.sub(&em1_neg, RoundingMode::NearestEven);
            let two = BigFloat::try_from_i64_exact(2, w).expect("precision >= 1");
            diff.div(&two, RoundingMode::NearestEven).0
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
    fn sinh_zero() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.sinh(RoundingMode::NearestEven);
            assert!(r.is_zero());
            assert_eq!(r.is_sign_negative(), matches!(s, Sign::Negative));
        }
    }

    #[test]
    fn sinh_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.sinh(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn sinh_neg_inf() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, _) = ni.sinh(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative());
    }

    #[test]
    fn sinh_one_matches_definition() {
        // sinh(1) = (e − 1/e)/2 ≈ 1.1752011936438014
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.sinh(RoundingMode::NearestEven);
        let (e, _) = one.exp(RoundingMode::NearestEven);
        let neg_one = BigFloat::try_from_i64_exact(-1, 113).unwrap();
        let (e_neg, _) = neg_one.exp(RoundingMode::NearestEven);
        let (diff, _) = e.sub(&e_neg, RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (expected, _) = diff.div(&two, RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, 113 - 12));
    }

    #[test]
    fn sinh_negation() {
        // sinh(-x) = -sinh(x)
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let neg_two = BigFloat::try_from_i64_exact(-2, 113).unwrap();
        let (a, _) = two.sinh(RoundingMode::NearestEven);
        let (b, _) = neg_two.sinh(RoundingMode::NearestEven);
        let neg_a = a.negated();
        assert!(close_at(&b, &neg_a, 113 - 12));
    }

    #[test]
    fn sinh_small_x() {
        // sinh(small) ≈ small.
        let p = 113u32;
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let mut x = one;
        for _ in 0..30 {
            x = x.div(&two, RoundingMode::NearestEven).0;
        }
        let (r, _) = x.sinh(RoundingMode::NearestEven);
        // sinh(x) - x = x³/6 + ..., relative error from "sinh(x) = x" is x²/6.
        // For x = 2^-30, that's 2^-60 / 6 ≈ 2^-63. So they agree to ~60 bits.
        assert!(close_at(&r, &x, 58));
    }

    #[test]
    fn sinh_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.sinh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn sinh_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.sinh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
