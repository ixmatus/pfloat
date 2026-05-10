//! `beta(a, b) = Γ(a) · Γ(b) / Γ(a + b)`: the Euler beta function.
//!
//! Implementation: route through `lgamma`:
//!
//! ```text
//! ln β(a, b) = lgamma(a) + lgamma(b) − lgamma(a + b)
//! β(a, b)    = sign · exp(ln β(a, b))
//! ```
//!
//! The combined sign is the product of the three Γ signs. For
//! positive `a, b`, all three are positive and the sign is `+`.
//!
//! Slice 4c restricts the kernel to `a, b > 0` for full numerical
//! support. Non-positive integer inputs hit Γ poles and the
//! division produces an indeterminate form; the kernel returns
//! `qNaN + INVALID` in those cases rather than attempting the
//! delicate cancellation analysis. Negative non-integer inputs
//! that produce a well-defined result are also coerced to
//! `qNaN + INVALID` for now; a follow-up slice can extend the
//! domain by tracking signs explicitly.
//!
//! Special cases:
//!
//! - `beta(NaN, _) = beta(_, NaN) = NaN`; `sNaN` raises `INVALID`.
//! - `beta(a, b)` with `a ≤ 0` or `b ≤ 0`: `qNaN + INVALID`.
//! - `beta(+∞, b)` finite positive `b`: `+0`.
//! - `beta(a, +∞)` finite positive `a`: `+0`.

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
    /// `beta(self, other)` rounded under `mode` to
    /// `max(self.precision, other.precision)`.
    #[must_use]
    pub fn beta(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let target = self.precision.max(other.precision);
        self.beta_round(other, target, mode)
            .expect("max of two valid precisions is valid")
    }

    /// `beta(self, other)` with explicit result precision.
    pub fn beta_round(
        &self,
        other: &Self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(beta_kernel(self, other, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `beta(self, other)` for `FixedFloat`. Delegates to
    /// [`BigFloat::beta`].
    #[must_use]
    pub fn beta(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().beta(&other.to_big(), mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn beta_kernel(
    a: &BigFloat,
    b: &BigFloat,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    if a.is_signaling_nan() || b.is_signaling_nan() {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }
    if a.is_nan() {
        let nan =
            BigFloat::try_new_quiet_nan(a.sign(), target_precision, &[]).expect("precision >= 1");
        return (nan, Status::OK);
    }
    if b.is_nan() {
        let nan =
            BigFloat::try_new_quiet_nan(b.sign(), target_precision, &[]).expect("precision >= 1");
        return (nan, Status::OK);
    }

    // Infinity: beta(+∞, finite_positive) = +0.
    if (a.is_infinite() && matches!(a.sign(), Sign::Positive) && is_finite_positive(b))
        || (b.is_infinite() && matches!(b.sign(), Sign::Positive) && is_finite_positive(a))
    {
        return (
            BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1"),
            Status::OK,
        );
    }

    // Domain check: both a, b must be finite and strictly positive.
    if !is_finite_positive(a) || !is_finite_positive(b) {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }

    let working_prec = target_precision
        .saturating_add(64)
        .min(target_precision.saturating_add(512));
    let (lg_a, _) = a
        .lgamma_round(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1");
    let (lg_b, _) = b
        .lgamma_round(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1");
    let (sum, _) = a.add(b, RoundingMode::NearestEven);
    let (lg_sum, _) = sum
        .lgamma_round(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1");
    let (lg_a_plus_b, _) = lg_a.add(&lg_b, RoundingMode::NearestEven);
    let (lb, _) = lg_a_plus_b.sub(&lg_sum, RoundingMode::NearestEven);
    let (result, _) = lb.exp(RoundingMode::NearestEven);
    let (rounded, status) = result
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    auto_raise(status);
    (rounded, status)
}

fn is_finite_positive(x: &BigFloat) -> bool {
    matches!(
        &x.class,
        Class::Normal {
            sign: Sign::Positive,
            ..
        }
    )
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
    fn beta_2_3_is_one_twelfth() {
        // β(2, 3) = 1/12.
        let a = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let b = BigFloat::try_from_i64_exact(3, 113).unwrap();
        let (r, _) = a.beta(&b, RoundingMode::NearestEven);
        let twelve = BigFloat::try_from_i64_exact(12, 113).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (expected, _) = one.div(&twelve, RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn beta_half_half_is_pi() {
        // β(1/2, 1/2) = π.
        let half = BigFloat::parse_str("0.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = half.beta(&half, RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "3.1415926535897932384626433832795028841971693993751",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn beta_3_5() {
        // β(3, 5) = 1/105 ≈ 0.00952380952380952.
        let a = BigFloat::try_from_i64_exact(3, 113).unwrap();
        let b = BigFloat::try_from_i64_exact(5, 113).unwrap();
        let (r, _) = a.beta(&b, RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "0.0095238095238095238095238095238095238095238095238095",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn beta_is_symmetric() {
        // β(a, b) = β(b, a).
        let a = BigFloat::parse_str("2.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let b = BigFloat::parse_str("3.7", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (ab, _) = a.beta(&b, RoundingMode::NearestEven);
        let (ba, _) = b.beta(&a, RoundingMode::NearestEven);
        assert!(close_at(&ab, &ba, 80));
    }

    #[test]
    fn beta_negative_is_invalid() {
        let a = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        let b = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, status) = a.beta(&b, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn beta_zero_is_invalid() {
        let a = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let b = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, status) = a.beta(&b, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn beta_pos_inf_is_zero() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let b = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, _) = pi.beta(&b, RoundingMode::NearestEven);
        assert!(r.is_zero() && r.is_sign_positive());
    }

    #[test]
    fn beta_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, _) = q.beta(&one, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        let (r, _) = one.beta(&q, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn beta_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, status) = sn.beta(&one, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
