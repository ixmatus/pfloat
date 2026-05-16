//! `asin(x)`: inverse sine, defined for `|x| ≤ 1`.
//!
//! Identity used: `asin(x) = 2 · atan(x / (1 + sqrt(1 − x²)))`. The
//! divisor is bounded below by `1 + 0 = 1` for all `|x| ≤ 1`, so
//! the argument to `atan` never blows up — even at `x = ±1`, where
//! the divisor is `1` and the argument is `±1`, giving
//! `2 · atan(±1) = ±π/2` as required.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `asin(±0) = ±0`.
//! - `asin(±1) = ±π/2`.
//! - `asin(x) = qNaN + INVALID` for `|x| > 1`, including `±∞`.
//! - `asin(NaN) = NaN`; `sNaN` raises `INVALID`.

use core::cmp::Ordering;

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::atan::atan_finite_unsigned;
use super::pi_over_2_at;

impl BigFloat {
    /// `asin(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn asin(&self, mode: RoundingMode) -> (Self, Status) {
        self.asin_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `asin(self)` with explicit result precision.
    pub fn asin_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(asin_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `asin(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::asin`].
    #[must_use]
    pub fn asin(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().asin(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn asin_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
        Class::Infinity { .. } => {
            // |∞| > 1: domain error.
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        Class::Normal { .. } => {}
    }

    // Domain dispatch on |x| vs 1.
    let sign = x.sign();
    let abs_x = x.abs();
    let one_at_input = BigFloat::try_from_i64_exact(1, x.precision).expect("precision >= 1");
    match abs_x.partial_cmp(&one_at_input).0 {
        Some(Ordering::Greater) => {
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        Some(Ordering::Equal) => {
            // asin(±1) = ±π/2.
            let pi_2 = pi_over_2_at(target_precision);
            let signed = if matches!(sign, Sign::Negative) {
                pi_2.negated()
            } else {
                pi_2
            };
            return (signed, Status::OK);
        }
        _ => {}
    }

    let working_prec = target_precision.saturating_add(64);
    let abs_x_w = abs_x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let two = BigFloat::try_from_i64_exact(2, working_prec).expect("precision >= 1");

    let (x_sq, _) = abs_x_w.mul(&abs_x_w, RoundingMode::NearestEven);
    let (one_minus_sq, _) = one.sub(&x_sq, RoundingMode::NearestEven);
    let (s, _) = one_minus_sq.sqrt(RoundingMode::NearestEven);
    let (denom, _) = one.add(&s, RoundingMode::NearestEven);
    let (y, _) = abs_x_w.div(&denom, RoundingMode::NearestEven);
    let atan_y = atan_finite_unsigned(&y, working_prec);
    let (twice, _) = two.mul(&atan_y, RoundingMode::NearestEven);

    let signed = if matches!(sign, Sign::Negative) {
        twice.negated()
    } else {
        twice
    };
    let (rounded, status) = signed
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    auto_raise(status);
    (rounded, status)
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn asin_zero() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.asin(RoundingMode::NearestEven);
            assert!(r.is_zero());
            assert_eq!(r.is_sign_negative(), matches!(s, Sign::Negative));
        }
    }

    #[test]
    fn asin_one_is_pi_over_2() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.asin(RoundingMode::NearestEven);
        let pi_2 = super::super::pi_over_2_at(113);
        assert!(close_at(&r, &pi_2, 100));
    }

    #[test]
    fn asin_neg_one_is_neg_pi_over_2() {
        let neg_one = BigFloat::try_from_i64_exact(-1, 113).unwrap();
        let (r, _) = neg_one.asin(RoundingMode::NearestEven);
        let pi_2 = super::super::pi_over_2_at(113);
        let neg = pi_2.negated();
        assert!(close_at(&r, &neg, 100));
    }

    #[test]
    fn asin_above_one_is_invalid() {
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, status) = two.asin(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn asin_below_neg_one_is_invalid() {
        let neg_two = BigFloat::try_from_i64_exact(-2, 53).unwrap();
        let (r, status) = neg_two.asin(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn asin_pos_inf_is_invalid() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, status) = pi.asin(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn asin_half_is_pi_over_6() {
        // asin(0.5) = π/6
        let half = BigFloat::parse_str("0.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = half.asin(RoundingMode::NearestEven);
        let pi_2 = super::super::pi_over_2_at(113);
        let three = BigFloat::try_from_i64_exact(3, 113).unwrap();
        let (pi_6, _) = pi_2.div(&three, RoundingMode::NearestEven);
        assert!(close_at(&r, &pi_6, 113 - 12));
    }

    #[test]
    fn asin_sin_round_trip() {
        // asin(sin(x)) = x for x ∈ [−π/2, π/2].
        let p = 113u32;
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (s, _) = one.sin(RoundingMode::NearestEven);
        let (back, _) = s.asin(RoundingMode::NearestEven);
        assert!(close_at(&back, &one, p - 16));
    }

    #[test]
    fn asin_is_odd() {
        let p = 113u32;
        let half = BigFloat::parse_str("0.5", p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let neg_half = half.negated();
        let (a, _) = half.asin(RoundingMode::NearestEven);
        let (b, _) = neg_half.asin(RoundingMode::NearestEven);
        let neg_a = a.negated();
        assert!(close_at(&b, &neg_a, p - 12));
    }

    #[test]
    fn asin_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.asin(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn asin_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.asin(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
