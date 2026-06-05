//! `acosh(x) = ln(x + sqrt(x² − 1))`: inverse hyperbolic cosine,
//! defined for `x ≥ 1`.
//!
//! Identity used: `acosh(x) = log1p((x − 1) + sqrt((x − 1)(x + 1)))`.
//! The `log1p` form avoids the cancellation near `x = 1` (the
//! naive `ln(x + sqrt(x² − 1))` collapses to `ln(1 + tiny)` and
//! loses leading bits). For large `x`, the formula reduces to
//! `log1p(2x − 1) ≈ ln(2x)` as expected.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `acosh(1) = +0`.
//! - `acosh(x) = qNaN + INVALID` for `x < 1` (including `−∞`).
//! - `acosh(+∞) = +∞`.
//! - `acosh(NaN) = NaN`; `sNaN` raises `INVALID`.

use core::cmp::Ordering;

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::ziv::ziv_round;
use super::ziv_calibration::ACOSH_ERROR_GUARD;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `acosh(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn acosh(&self, mode: RoundingMode) -> (Self, Status) {
        self.acosh_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `acosh(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.27, ADR-0038).
    pub fn acosh_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(acosh_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `acosh(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::acosh`].
    #[must_use]
    pub fn acosh(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().acosh(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn acosh_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        Class::Zero { .. } => {
            // 0 < 1: domain error.
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        Class::Normal { .. } => {}
    }

    // Domain: x ≥ 1.
    let one_at_input = BigFloat::try_from_i64_exact(1, x.precision).expect("precision >= 1");
    match x.partial_cmp(&one_at_input).0 {
        Some(Ordering::Less) => {
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        Some(Ordering::Equal) => {
            return (
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1"),
                Status::OK,
            );
        }
        _ => {}
    }

    // Ziv-driven correct rounding under every IEEE mode. The eval
    // closure runs the existing log1p-based identity
    // `acosh(x) = log1p((x − 1) + sqrt((x − 1)(x + 1)))` at
    // working precision `w` under NE; the outer envelope certifies
    // the rounding-mode interval test. log1p is Ziv-driven (slice
    // p1.24) and the identity avoids the near-1 cancellation of
    // the naive `ln(x + sqrt(x² − 1))` form.
    let (result, status) = ziv_round(
        |w| {
            let x_w = x
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let one = BigFloat::try_from_i64_exact(1, w).expect("precision >= 1");
            let (x_minus_one, _) = x_w.sub(&one, RoundingMode::NearestEven);
            let (x_plus_one, _) = x_w.add(&one, RoundingMode::NearestEven);
            let (prod, _) = x_minus_one.mul(&x_plus_one, RoundingMode::NearestEven);
            let (s, _) = prod.sqrt(RoundingMode::NearestEven);
            let (arg, _) = x_minus_one.add(&s, RoundingMode::NearestEven);
            arg.log1p(RoundingMode::NearestEven).0
        },
        target_precision,
        mode,
        ACOSH_ERROR_GUARD,
    );
    // acosh(x) for finite normal x > 1 is transcendental (Lindemann–
    // Weierstrass), hence irrational, hence INEXACT even where it rounds
    // onto a grid value (pf-uqd1, ADR-0063). acosh(1) = 0 is the only
    // exact input and is dispatched above.
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
    fn acosh_one_is_zero() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, _) = one.acosh(RoundingMode::NearestEven);
        assert!(r.is_zero());
    }

    #[test]
    fn acosh_below_one_is_invalid() {
        let half = BigFloat::parse_str("0.5", 53, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, status) = half.acosh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn acosh_zero_is_invalid() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, status) = z.acosh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn acosh_negative_is_invalid() {
        let neg = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        let (r, status) = neg.acosh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn acosh_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.acosh(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn acosh_neg_inf_is_invalid() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, status) = ni.acosh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn acosh_cosh_round_trip() {
        // acosh(cosh(x)) = |x| for x with cosh(x) > 1.
        let p = 113u32;
        for n in &[1i64, 2, 3, 5, 10] {
            let x = BigFloat::try_from_i64_exact(*n, p).unwrap();
            let (c, _) = x.cosh(RoundingMode::NearestEven);
            let (back, _) = c.acosh(RoundingMode::NearestEven);
            assert!(
                close_at(&back, &x, p - 16),
                "acosh(cosh({n})) = {back}, expected {x}"
            );
        }
    }

    #[test]
    fn acosh_two() {
        // acosh(2) = ln(2 + sqrt(3)) ≈ 1.3169578969248166
        let p = 113u32;
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let (r, _) = two.acosh(RoundingMode::NearestEven);
        let three = BigFloat::try_from_i64_exact(3, p).unwrap();
        let (sqrt3, _) = three.sqrt(RoundingMode::NearestEven);
        let (arg, _) = two.add(&sqrt3, RoundingMode::NearestEven);
        let (expected, _) = arg.ln(RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, p - 12));
    }

    #[test]
    fn acosh_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.acosh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn acosh_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.acosh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
