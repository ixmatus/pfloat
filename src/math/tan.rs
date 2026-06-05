//! `tan(x) = sin(x) / cos(x)`: trigonometric tangent.
//!
//! Implemented as the ratio of `sin` and `cos` after a single
//! shared Payne-Hanek reduction. Quadrants 0 and 2 use `sin(r) /
//! cos(r)`; quadrants 1 and 3 swap the roles and negate, yielding
//! `−cos(r) / sin(r)`.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `tan(±0) = ±0`.
//! - `tan(±∞) = qNaN + INVALID`.
//! - `tan(NaN) = NaN`; `sNaN` raises `INVALID`.
//! - `|x|` past the reduction table budget: `qNaN + INVALID`.
//!
//! Note: `tan(x)` diverges at the odd multiples of `π/2`. Because
//! the reduction stores `r` at finite working precision, the
//! division `sin(r) / cos(r)` near these poles produces a very
//! large but finite result rather than `±∞`. This matches MPFR's
//! behavior and IEEE 754-2019's expectation for "transcendental"
//! tangent.
//!
//! Correctly rounded under every IEEE 754-2019 rounding mode via
//! the shared [`crate::math::ziv::ziv_round`] driver (slice p1.26,
//! ADR-0038). The ratio composition runs inside the eval closure
//! at each Ziv working precision; the range-cap NaN check
//! pre-empts Ziv at the maximum working precision the driver
//! could request.

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
use super::ziv_calibration::TAN_ERROR_GUARD;

impl BigFloat {
    /// `tan(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn tan(&self, mode: RoundingMode) -> (Self, Status) {
        self.tan_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `tan(self)` with explicit result precision.
    pub fn tan_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(tan_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `tan(self)` for `FixedFloat`. Delegates to [`BigFloat::tan`].
    #[must_use]
    pub fn tan(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().tan(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn tan_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        Class::Normal { .. } => {}
    }

    // Range-cap check at the Ziv first-iteration working precision
    // (pf-1axr; see `sin.rs` for the full rationale).
    let ziv_first_working = target_precision.saturating_add(ZIV_BASE_GUARD);
    if reduce(x, ziv_first_working).is_none() {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }

    let (result, status) = ziv_round(
        |w| match reduce(x, w) {
            Some(Reduction { quadrant, r }) => {
                let s = sin_taylor(&r, w);
                let c = cos_taylor(&r, w);
                match quadrant {
                    0 | 2 => s.div(&c, RoundingMode::NearestEven).0,
                    _ => {
                        // quadrants 1, 3: tan(x) = −cos(r) / sin(r).
                        let neg_c = c.negated();
                        neg_c.div(&s, RoundingMode::NearestEven).0
                    }
                }
            }
            None => BigFloat::try_new_quiet_nan(Sign::Positive, w, &[]).expect("precision >= 1"),
        },
        target_precision,
        mode,
        TAN_ERROR_GUARD,
    );
    // Post-Ziv NaN-to-INVALID surfacing (pf-1axr, sin.rs precedent).
    if matches!(result.class, Class::Nan { .. }) && !status.invalid() {
        let merged = status.merge(Status::INVALID);
        auto_raise(Status::INVALID);
        return (result, merged);
    }
    // tan(x) for finite normal x is transcendental (Lindemann–
    // Weierstrass), hence irrational, hence INEXACT even where it rounds
    // onto a grid value (pf-uqd1, ADR-0063). tan(±0) = ±0 is the only
    // exact input and is special-cased above.
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
    fn tan_zero() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.tan(RoundingMode::NearestEven);
            assert!(r.is_zero());
            assert_eq!(r.is_sign_negative(), matches!(s, Sign::Negative));
        }
    }

    #[test]
    fn tan_pi_is_zero() {
        let pi = super::super::pi_at(113);
        let (r, _) = pi.tan(RoundingMode::NearestEven);
        let zero = BigFloat::try_new_zero(Sign::Positive, 113).unwrap();
        assert!(close_at(&r, &zero, 100));
    }

    #[test]
    fn tan_one_matches_definition() {
        // tan(1) = sin(1)/cos(1)
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.tan(RoundingMode::NearestEven);
        let (s, _) = one.sin(RoundingMode::NearestEven);
        let (c, _) = one.cos(RoundingMode::NearestEven);
        let (expected, _) = s.div(&c, RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, 113 - 16));
    }

    #[test]
    fn tan_pi_over_4_is_one() {
        let pi_4 = super::super::pi_over_2_at(113);
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (pi_4, _) = pi_4.div(&two, RoundingMode::NearestEven);
        let (r, _) = pi_4.tan(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        assert!(close_at(&r, &one, 113 - 16));
    }

    #[test]
    fn tan_is_odd() {
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let neg_two = BigFloat::try_from_i64_exact(-2, 113).unwrap();
        let (a, _) = two.tan(RoundingMode::NearestEven);
        let (b, _) = neg_two.tan(RoundingMode::NearestEven);
        let neg_a = a.negated();
        assert!(close_at(&b, &neg_a, 113 - 16));
    }

    #[test]
    fn tan_pos_inf_is_invalid() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, status) = pi.tan(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn tan_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.tan(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn tan_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.tan(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
