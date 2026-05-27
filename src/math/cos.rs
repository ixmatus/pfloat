//! `cos(x)`: trigonometric cosine.
//!
//! Argument reduction (see [`super::trig_reduce`]) gives a quadrant
//! `q` and reduced argument `r ∈ [−π/4, π/4]`. Then:
//!
//! - `q = 0`: `cos(x) = cos(r)`
//! - `q = 1`: `cos(x) = −sin(r)`
//! - `q = 2`: `cos(x) = −cos(r)`
//! - `q = 3`: `cos(x) = +sin(r)`
//!
//! Taylor series for `sin(r)` and `cos(r)` are shared with
//! [`super::sin`].
//!
//! Correctly rounded under every IEEE 754-2019 rounding mode via
//! the shared [`crate::math::ziv::ziv_round`] driver (slice p1.26,
//! ADR-0038). Same shape as `sin`: range-cap pre-check, then the
//! composition runs inside the eval closure at each Ziv working
//! precision.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `cos(±0) = 1`.
//! - `cos(±∞) = qNaN + INVALID`.
//! - `cos(NaN) = NaN`; `sNaN` raises `INVALID`.
//! - `|x|` past the reduction table budget: `qNaN + INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::sin::{cos_taylor, sin_taylor};
use super::trig_reduce::{reduce, Reduction};
use super::ziv::{ziv_round, ZIV_BASE_GUARD};
use super::ziv_calibration::COS_ERROR_GUARD;

impl BigFloat {
    /// `cos(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn cos(&self, mode: RoundingMode) -> (Self, Status) {
        self.cos_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `cos(self)` with explicit result precision.
    pub fn cos_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(cos_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `cos(self)` for `FixedFloat`. Delegates to [`BigFloat::cos`].
    #[must_use]
    pub fn cos(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().cos(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn cos_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        Class::Normal { .. } => {}
    }

    // Range-cap check at the Ziv first-iteration working precision
    // (sin family precedent; see `sin.rs` for the full pf-1axr
    // rationale on why this is `target + ZIV_BASE_GUARD` rather than
    // `target + ZIV_GUARD_CAP`).
    let ziv_first_working = target_precision.saturating_add(ZIV_BASE_GUARD);
    if reduce(x, ziv_first_working).is_none() {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }

    let (result, status) = ziv_round(
        |w| match reduce(x, w) {
            Some(Reduction { quadrant, r }) => match quadrant {
                0 => cos_taylor(&r, w),
                1 => sin_taylor(&r, w).negated(),
                2 => cos_taylor(&r, w).negated(),
                _ => sin_taylor(&r, w),
            },
            None => BigFloat::try_new_quiet_nan(Sign::Positive, w, &[])
                .expect("precision >= 1"),
        },
        target_precision,
        mode,
        COS_ERROR_GUARD,
    );
    // Post-Ziv NaN-to-INVALID surfacing (pf-1axr, sin.rs precedent).
    if matches!(result.class, Class::Nan { .. }) && !status.invalid() {
        let merged = status.merge(Status::INVALID);
        auto_raise(Status::INVALID);
        return (result, merged);
    }
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
    fn cos_zero_is_one() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.cos(RoundingMode::NearestEven);
            let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
            assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
        }
    }

    #[test]
    fn cos_pi_is_neg_one() {
        let pi = super::super::pi_at(113);
        let (r, _) = pi.cos(RoundingMode::NearestEven);
        let neg_one = BigFloat::try_from_i64_exact(-1, 113).unwrap();
        assert!(close_at(&r, &neg_one, 113 - 12));
    }

    #[test]
    fn cos_pi_over_2_is_zero() {
        let pi_2 = super::super::pi_over_2_at(113);
        let (r, _) = pi_2.cos(RoundingMode::NearestEven);
        let zero = BigFloat::try_new_zero(Sign::Positive, 113).unwrap();
        assert!(close_at(&r, &zero, 100));
    }

    #[test]
    fn cos_one() {
        // cos(1) ≈ 0.5403023058681398
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.cos(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "0.54030230586813971740093660744297661",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 100));
    }

    #[test]
    fn cos_is_even() {
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let neg_two = BigFloat::try_from_i64_exact(-2, 113).unwrap();
        let (a, _) = two.cos(RoundingMode::NearestEven);
        let (b, _) = neg_two.cos(RoundingMode::NearestEven);
        assert!(close_at(&a, &b, 113 - 12));
    }

    #[test]
    fn sin_cos_pythagorean() {
        // sin²(x) + cos²(x) = 1
        let x = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (s, _) = x.sin(RoundingMode::NearestEven);
        let (c, _) = x.cos(RoundingMode::NearestEven);
        let (s2, _) = s.mul(&s, RoundingMode::NearestEven);
        let (c2, _) = c.mul(&c, RoundingMode::NearestEven);
        let (sum, _) = s2.add(&c2, RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        assert!(close_at(&sum, &one, 100));
    }

    #[test]
    fn cos_pos_inf_is_invalid() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, status) = pi.cos(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn cos_neg_inf_is_invalid() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, status) = ni.cos(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn cos_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.cos(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn cos_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.cos(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
