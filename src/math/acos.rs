//! `acos(x)`: inverse cosine, defined for `|x| ≤ 1`, returning a
//! value in `[0, π]`.
//!
//! Identity used:
//!
//! - For `x ≥ 0`: `acos(x) = 2 · atan(sqrt((1 − x)/(1 + x)))`.
//! - For `x < 0`: `acos(x) = π − 2 · atan(sqrt((1 + x)/(1 − x)))`.
//!
//! The two-branch form keeps the `atan` argument bounded in `[0, ∞)`
//! and avoids the catastrophic cancellation that `π/2 − asin(x)`
//! would suffer near `x = 1` (where `asin → π/2`).
//!
//! Correctly rounded under every IEEE 754-2019 rounding mode via
//! the shared [`crate::math::ziv::ziv_round`] driver (slice p1.25,
//! ADR-0038). The `x = 0` and `x = -1` special cases return
//! `π/2` and `π` respectively, rounded under the caller's mode via
//! [`super::pi_over_2_at_round`] and [`super::pi_at_round`].
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `acos(+1) = +0`.
//! - `acos(−1) = π`.
//! - `acos(x) = qNaN + INVALID` for `|x| > 1`, including `±∞`.
//! - `acos(NaN) = NaN`; `sNaN` raises `INVALID`.

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
use super::ziv::ziv_round;
use super::{pi_at, pi_at_round, pi_over_2_at_round};

impl BigFloat {
    /// `acos(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn acos(&self, mode: RoundingMode) -> (Self, Status) {
        self.acos_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `acos(self)` with explicit result precision.
    pub fn acos_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(acos_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `acos(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::acos`].
    #[must_use]
    pub fn acos(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().acos(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn acos_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            // acos(0) = π/2. Round under the caller's mode (slice
            // p1.25 surfaced the NE-only `pi_over_2_at` return as
            // 1-ULP wrong under TZ/TN; the helper boosts precision
            // and rounds to target under mode).
            let (rounded, status) = pi_over_2_at_round(target_precision, mode);
            crate::status::auto_raise(status);
            return (rounded, status);
        }
        Class::Infinity { .. } => {
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        Class::Normal { .. } => {}
    }

    let one_at_input = BigFloat::try_from_i64_exact(1, x.precision).expect("precision >= 1");
    let neg_one_at_input = BigFloat::try_from_i64_exact(-1, x.precision).expect("precision >= 1");
    match x.partial_cmp(&one_at_input).0 {
        Some(Ordering::Equal) => {
            return (
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1"),
                Status::OK,
            );
        }
        Some(Ordering::Greater) => {
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        _ => {}
    }
    match x.partial_cmp(&neg_one_at_input).0 {
        Some(Ordering::Equal) => {
            // acos(-1) = π. Round under the caller's mode.
            let (rounded, status) = pi_at_round(target_precision, mode);
            crate::status::auto_raise(status);
            return (rounded, status);
        }
        Some(Ordering::Less) => {
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        _ => {}
    }

    // Ziv-driven correct rounding under every IEEE mode. The
    // eval closure carries the two-branch composition: for x ≥ 0
    // use `2·atan(sqrt((1−x)/(1+x)))`; for x < 0 use
    // `π − 2·atan(sqrt((1+x)/(1−x)))`. Both forms avoid the
    // catastrophic cancellation `π/2 − asin(x)` would suffer near
    // x = ±1. The branch dispatch is on `x.sign()` (Class-level,
    // independent of working precision) so it stays inside eval.
    let is_negative = matches!(x.sign(), Sign::Negative);
    let (result, status) = ziv_round(
        |w| {
            let x_w = x
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let one = BigFloat::try_from_i64_exact(1, w).expect("precision >= 1");
            let two = BigFloat::try_from_i64_exact(2, w).expect("precision >= 1");

            if is_negative {
                let (one_plus_x, _) = one.add(&x_w, RoundingMode::NearestEven);
                let (one_minus_x, _) = one.sub(&x_w, RoundingMode::NearestEven);
                let (ratio, _) = one_plus_x.div(&one_minus_x, RoundingMode::NearestEven);
                let (arg, _) = ratio.sqrt(RoundingMode::NearestEven);
                let at = atan_finite_unsigned(&arg, w);
                let (twice, _) = two.mul(&at, RoundingMode::NearestEven);
                let pi = pi_at(w);
                let (diff, _) = pi.sub(&twice, RoundingMode::NearestEven);
                diff
            } else {
                let (one_minus_x, _) = one.sub(&x_w, RoundingMode::NearestEven);
                let (one_plus_x, _) = one.add(&x_w, RoundingMode::NearestEven);
                let (ratio, _) = one_minus_x.div(&one_plus_x, RoundingMode::NearestEven);
                let (arg, _) = ratio.sqrt(RoundingMode::NearestEven);
                let at = atan_finite_unsigned(&arg, w);
                let (twice, _) = two.mul(&at, RoundingMode::NearestEven);
                twice
            }
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
    fn acos_zero_is_pi_over_2() {
        let z = BigFloat::try_new_zero(Sign::Positive, 113).unwrap();
        let (r, _) = z.acos(RoundingMode::NearestEven);
        let pi_2 = super::super::pi_over_2_at(113);
        assert!(close_at(&r, &pi_2, 100));
    }

    #[test]
    fn acos_one_is_zero() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.acos(RoundingMode::NearestEven);
        assert!(r.is_zero());
    }

    #[test]
    fn acos_neg_one_is_pi() {
        let neg_one = BigFloat::try_from_i64_exact(-1, 113).unwrap();
        let (r, _) = neg_one.acos(RoundingMode::NearestEven);
        let pi = super::super::pi_at(113);
        assert!(close_at(&r, &pi, 100));
    }

    #[test]
    fn acos_above_one_is_invalid() {
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, status) = two.acos(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn acos_below_neg_one_is_invalid() {
        let neg_two = BigFloat::try_from_i64_exact(-2, 53).unwrap();
        let (r, status) = neg_two.acos(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn acos_half_is_pi_over_3() {
        let half = BigFloat::parse_str("0.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = half.acos(RoundingMode::NearestEven);
        let pi = super::super::pi_at(113);
        let three = BigFloat::try_from_i64_exact(3, 113).unwrap();
        let (pi_3, _) = pi.div(&three, RoundingMode::NearestEven);
        assert!(close_at(&r, &pi_3, 113 - 12));
    }

    #[test]
    fn acos_cos_round_trip() {
        // acos(cos(x)) = x for x ∈ [0, π].
        let p = 113u32;
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (c, _) = one.cos(RoundingMode::NearestEven);
        let (back, _) = c.acos(RoundingMode::NearestEven);
        assert!(close_at(&back, &one, p - 16));
    }

    #[test]
    fn acos_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.acos(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn acos_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.acos(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
