//! `asinh(x) = ln(x + sqrt(x² + 1))`: inverse hyperbolic sine.
//!
//! Identity used: `asinh(x) = sign(x) · log1p(|x| + |x|² /
//! (sqrt(|x|² + 1) + 1))`. The formula has no cancellation for
//! `|x| ≥ 0`: every term in the `log1p` argument is non-negative,
//! and `sqrt(|x|² + 1) + 1 ≥ 2 > 0` for any `x`. Naively evaluating
//! `ln(x + sqrt(x² + 1))` for large negative `x` cancels leading
//! bits because `x + sqrt(x² + 1) → 0+`; computing on `|x|` and
//! applying the sign avoids the issue.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `asinh(±0) = ±0`.
//! - `asinh(±∞) = ±∞`.
//! - `asinh(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `asinh(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn asinh(&self, mode: RoundingMode) -> (Self, Status) {
        self.asinh_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `asinh(self)` with explicit result precision.
    pub fn asinh_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(asinh_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `asinh(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::asinh`].
    #[must_use]
    pub fn asinh(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().asinh(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn asinh_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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

    let working_prec = target_precision.saturating_add(64).min(1024);
    let sign = x.sign();
    let abs_x = x.abs();
    let x_w = abs_x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let (x_sq, _) = x_w.mul(&x_w, RoundingMode::NearestEven);
    let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let (x_sq_plus_one, _) = x_sq.add(&one, RoundingMode::NearestEven);
    let (s, _) = x_sq_plus_one.sqrt(RoundingMode::NearestEven);
    let (s_plus_one, _) = s.add(&one, RoundingMode::NearestEven);
    let (correction, _) = x_sq.div(&s_plus_one, RoundingMode::NearestEven);
    let (arg, _) = x_w.add(&correction, RoundingMode::NearestEven);
    let (lp, _) = arg.log1p(RoundingMode::NearestEven);
    let signed = if matches!(sign, Sign::Negative) {
        lp.negated()
    } else {
        lp
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
    fn asinh_zero() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.asinh(RoundingMode::NearestEven);
            assert!(r.is_zero());
            assert_eq!(r.is_sign_negative(), matches!(s, Sign::Negative));
        }
    }

    #[test]
    fn asinh_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.asinh(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn asinh_neg_inf() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, _) = ni.asinh(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative());
    }

    #[test]
    fn asinh_sinh_round_trip() {
        // asinh(sinh(x)) = x for moderate x.
        let p = 113u32;
        for n in &[-3i64, -1, 1, 3, 5] {
            let x = BigFloat::try_from_i64_exact(*n, p).unwrap();
            let (s, _) = x.sinh(RoundingMode::NearestEven);
            let (back, _) = s.asinh(RoundingMode::NearestEven);
            assert!(
                close_at(&back, &x, p - 16),
                "asinh(sinh({n})) = {back}, expected {x}"
            );
        }
    }

    #[test]
    fn asinh_one() {
        // asinh(1) = ln(1 + sqrt(2)) ≈ 0.8813735870195429
        let p = 113u32;
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (r, _) = one.asinh(RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let (sqrt2, _) = two.sqrt(RoundingMode::NearestEven);
        let (arg, _) = one.add(&sqrt2, RoundingMode::NearestEven);
        let (expected, _) = arg.ln(RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, p - 12));
    }

    #[test]
    fn asinh_negation() {
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let neg_two = BigFloat::try_from_i64_exact(-2, 113).unwrap();
        let (a, _) = two.asinh(RoundingMode::NearestEven);
        let (b, _) = neg_two.asinh(RoundingMode::NearestEven);
        let neg_a = a.negated();
        assert!(close_at(&b, &neg_a, 113 - 8));
    }

    #[test]
    fn asinh_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.asinh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn asinh_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.asinh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
