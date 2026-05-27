//! `atanh(x) = (1/2) · ln((1 + x)/(1 − x))`: inverse hyperbolic
//! tangent, defined for `|x| < 1`.
//!
//! Identity used: `atanh(x) = (log1p(x) − log1p(−x)) / 2`. Both
//! `log1p` calls handle their small-argument cancellation regimes
//! internally; for `x → 0` the subtraction `log1p(x) − log1p(−x)`
//! reduces to `2x` in real math without leading-bit cancellation
//! in the `BigFloat` representation.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `atanh(±0) = ±0`.
//! - `atanh(±1) = ±∞ + DIV_BY_ZERO`.
//! - `atanh(x) = qNaN + INVALID` for `|x| > 1`, including `±∞`.
//! - `atanh(NaN) = NaN`; `sNaN` raises `INVALID`.

use core::cmp::Ordering;

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::ziv::ziv_round;
use super::ziv_calibration::ATANH_ERROR_GUARD;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `atanh(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn atanh(&self, mode: RoundingMode) -> (Self, Status) {
        self.atanh_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `atanh(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.27, ADR-0038).
    pub fn atanh_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(atanh_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `atanh(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::atanh`].
    #[must_use]
    pub fn atanh(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().atanh(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn atanh_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
        Some(Ordering::Equal) => {
            // atanh(±1) = ±∞ + DIV_BY_ZERO.
            let inf = BigFloat::try_new_infinity(sign, target_precision).expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            return (inf, Status::DIV_BY_ZERO);
        }
        Some(Ordering::Greater) => {
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        _ => {}
    }

    // Ziv-driven correct rounding under every IEEE mode. The eval
    // closure runs the cancellation-resistant identity
    // `atanh(x) = (log1p(x) − log1p(-x)) / 2` at working precision
    // `w`. log1p is Ziv-driven (slice p1.24) and each log1p
    // handles its own small-argument cancellation; the outer
    // subtraction adds the leading terms (no cancellation).
    let (result, status) = ziv_round(
        |w| {
            let x_w = x
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let neg_x = x_w.negated();
            let (lp_pos, _) = x_w.log1p(RoundingMode::NearestEven);
            let (lp_neg, _) = neg_x.log1p(RoundingMode::NearestEven);
            let (diff, _) = lp_pos.sub(&lp_neg, RoundingMode::NearestEven);
            let two = BigFloat::try_from_i64_exact(2, w).expect("precision >= 1");
            diff.div(&two, RoundingMode::NearestEven).0
        },
        target_precision,
        mode,
        ATANH_ERROR_GUARD,
    );
    auto_raise(status);
    (result, status)
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
    fn atanh_zero() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.atanh(RoundingMode::NearestEven);
            assert!(r.is_zero());
            assert_eq!(r.is_sign_negative(), matches!(s, Sign::Negative));
        }
    }

    #[test]
    fn atanh_one_is_pos_inf_div() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, status) = one.atanh(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
        assert!(status.div_by_zero());
    }

    #[test]
    fn atanh_neg_one_is_neg_inf_div() {
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        let (r, status) = neg_one.atanh(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative());
        assert!(status.div_by_zero());
    }

    #[test]
    fn atanh_above_one_is_invalid() {
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, status) = two.atanh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn atanh_pos_inf_is_invalid() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, status) = pi.atanh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn atanh_half() {
        // atanh(0.5) = (1/2)·ln(3) ≈ 0.5493061443340549
        let p = 113u32;
        let half = BigFloat::parse_str("0.5", p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = half.atanh(RoundingMode::NearestEven);
        let three = BigFloat::try_from_i64_exact(3, p).unwrap();
        let (ln3, _) = three.ln(RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let (expected, _) = ln3.div(&two, RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, p - 12));
    }

    #[test]
    fn atanh_tanh_round_trip() {
        // atanh(tanh(x)) = x for moderate x.
        let p = 113u32;
        for n in &[-3i64, -1, 1, 2, 4] {
            let x = BigFloat::try_from_i64_exact(*n, p).unwrap();
            let (t, _) = x.tanh(RoundingMode::NearestEven);
            let (back, _) = t.atanh(RoundingMode::NearestEven);
            assert!(
                close_at(&back, &x, p - 16),
                "atanh(tanh({n})) = {back}, expected {x}"
            );
        }
    }

    #[test]
    fn atanh_negation() {
        let p = 113u32;
        let half = BigFloat::parse_str("0.5", p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let neg_half = half.negated();
        let (a, _) = half.atanh(RoundingMode::NearestEven);
        let (b, _) = neg_half.atanh(RoundingMode::NearestEven);
        let neg_a = a.negated();
        assert!(close_at(&b, &neg_a, p - 12));
    }

    #[test]
    fn atanh_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.atanh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn atanh_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.atanh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
