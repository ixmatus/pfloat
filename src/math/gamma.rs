//! `gamma(x) = Γ(x)`: the gamma function. Defined for all real `x`
//! except the non-positive integers (simple poles).
//!
//! Implementation: compose `lgamma` with `exp`, applying the
//! correct sign. For positive `x`, `Γ(x) > 0`, so the sign is
//! always positive. For negative non-integer `x`, the sign
//! alternates with the floor of `|x|`: `Γ(x) > 0` for
//! `x ∈ (−2, −1) ∪ (−4, −3) ∪ …` and `Γ(x) < 0` for
//! `x ∈ (−1, 0) ∪ (−3, −2) ∪ …`. Concretely, the sign equals
//! `−sign(sin(πx))` on the negative reals.
//!
//! For exact positive integer inputs `n ≥ 1`, the kernel returns
//! `(n−1)!` exactly when the result fits in target precision.
//! Otherwise it goes through `exp(lgamma(x))`, which inherits the
//! `lgamma` precision.
//!
//! Special cases per IEEE 754-2019 §9.4:
//!
//! - `gamma(+0) = +∞ + DIV_BY_ZERO`.
//! - `gamma(−0) = −∞ + DIV_BY_ZERO`.
//! - `gamma(negative integer) = qNaN + INVALID` (pole).
//! - `gamma(+∞) = +∞`.
//! - `gamma(−∞) = qNaN + INVALID`.
//! - `gamma(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::lgamma::is_integer_test;
use super::pi_at;

impl BigFloat {
    /// `gamma(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn gamma(&self, mode: RoundingMode) -> (Self, Status) {
        self.gamma_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `gamma(self)` with explicit result precision.
    pub fn gamma_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(gamma_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `gamma(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::gamma`].
    #[must_use]
    pub fn gamma(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().gamma(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn gamma_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            // gamma(±0) = ±∞ + DIV_BY_ZERO.
            let inf = BigFloat::try_new_infinity(*sign, target_precision).expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            return (inf, Status::DIV_BY_ZERO);
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
        Class::Normal { .. } => {}
    }

    // Negative integer: pole, qNaN + INVALID.
    if matches!(x.sign(), Sign::Negative) && is_integer_test(x) {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }

    let working_prec = target_precision
        .saturating_add(64)
        .min(target_precision.saturating_add(512));
    let (ln_abs_gamma, _) = x
        .lgamma_round(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1");

    let (abs_gamma, _) = ln_abs_gamma.exp(RoundingMode::NearestEven);

    let result_sign = gamma_sign_of(x, working_prec);
    let result = if matches!(result_sign, Sign::Negative) {
        abs_gamma.negated()
    } else {
        abs_gamma
    };
    let (rounded, status) = result
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    auto_raise(status);
    (rounded, status)
}

/// Sign of `Γ(x)` for finite non-zero `x` that is not a negative
/// integer pole. Positive `x` always yields positive sign; negative
/// non-integer `x` alternates per the reflection
/// `Γ(x)·Γ(1−x) = π/sin(πx)` (`Γ(1−x) > 0` for `x < 1`, so the
/// sign of `Γ(x)` matches the sign of `sin(πx)` for negative `x`).
pub(super) fn gamma_sign_of(x: &BigFloat, working_prec: u32) -> Sign {
    if matches!(x.sign(), Sign::Positive) {
        return Sign::Positive;
    }
    let pi = pi_at(working_prec);
    let (pi_x, _) = pi.mul(x, RoundingMode::NearestEven);
    let (sin_val, _) = pi_x.sin(RoundingMode::NearestEven);
    if matches!(sin_val.sign(), Sign::Negative) {
        Sign::Negative
    } else {
        Sign::Positive
    }
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
    fn gamma_one_is_one() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.gamma(RoundingMode::NearestEven);
        let expected = BigFloat::try_from_i64_exact(1, 113).unwrap();
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn gamma_two_is_one() {
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (r, _) = two.gamma(RoundingMode::NearestEven);
        let expected = BigFloat::try_from_i64_exact(1, 113).unwrap();
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn gamma_five_is_twentyfour() {
        // Γ(5) = 4! = 24.
        let five = BigFloat::try_from_i64_exact(5, 113).unwrap();
        let (r, _) = five.gamma(RoundingMode::NearestEven);
        let expected = BigFloat::try_from_i64_exact(24, 113).unwrap();
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn gamma_half_is_sqrt_pi() {
        // Γ(1/2) = √π ≈ 1.7724538509055160272981674833411.
        let half = BigFloat::parse_str("0.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = half.gamma(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "1.7724538509055160272981674833411451827975494561224",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn gamma_pos_zero_is_pos_inf_div() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, status) = z.gamma(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
        assert!(status.div_by_zero());
    }

    #[test]
    fn gamma_neg_zero_is_neg_inf_div() {
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, status) = nz.gamma(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative());
        assert!(status.div_by_zero());
    }

    #[test]
    fn gamma_negative_integer_is_nan() {
        let neg_three = BigFloat::try_from_i64_exact(-3, 53).unwrap();
        let (r, status) = neg_three.gamma(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn gamma_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.gamma(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn gamma_neg_inf_is_invalid() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, status) = ni.gamma(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn gamma_neg_half_is_neg_two_sqrt_pi() {
        // Γ(-1/2) = -2√π ≈ -3.5449077018110320546.
        let neg_half = BigFloat::parse_str("-0.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = neg_half.gamma(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "-3.5449077018110320545963349666822903655950989122448",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 60));
    }

    #[test]
    fn gamma_neg_1_5_is_positive() {
        // Γ(-1.5) = 4√π/3 ≈ 2.36327180120735... (positive).
        let neg = BigFloat::parse_str("-1.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = neg.gamma(RoundingMode::NearestEven);
        assert!(r.is_sign_positive());
        let expected = BigFloat::parse_str(
            "2.3632718012073547030642233111215269103967326081632",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 60));
    }

    #[test]
    fn gamma_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.gamma(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn gamma_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.gamma(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
