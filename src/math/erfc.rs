//! `erfc(x) = 1 − erf(x) = (2/√π) · ∫ₓ^∞ e^(−t²) dt`: the
//! complementary error function.
//!
//! For `x ≥ 0` two regimes split the work:
//!
//! - **Small `x`:** evaluate `1 − erf_maclaurin(x)` at working
//!   precision boosted to absorb the `1 − tiny` cancellation as
//!   `erf(x) → 1`.
//!
//! - **Large `x`:** evaluate [`erfc_asymptotic`] (the
//!   `(e^(−x²)/(x √π)) · Σ …` series, truncated at the smallest
//!   term).
//!
//! For `x < 0` the symmetry `erfc(−x) = 2 − erfc(x)` brings the
//! argument back to the positive regime. The shared
//! [`super::erf::erf_maclaurin`] and `erfc_asymptotic` helpers
//! avoid `erf` ↔ `erfc` mutual recursion in the kernels.
//!
//! Special cases per IEEE 754-2019 §9.4.3:
//!
//! - `erfc(±0) = 1`.
//! - `erfc(+∞) = +0`, `erfc(−∞) = 2`.
//! - `erfc(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::erf::{asymptotic_threshold_exponent, erf_maclaurin};
use super::two_over_sqrt_pi_at;

impl BigFloat {
    /// `erfc(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn erfc(&self, mode: RoundingMode) -> (Self, Status) {
        self.erfc_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `erfc(self)` with explicit result precision.
    pub fn erfc_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(erfc_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `erfc(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::erfc`].
    #[must_use]
    pub fn erfc(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().erfc(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn erfc_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
        Class::Infinity {
            sign: Sign::Positive,
        } => {
            return (
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Infinity {
            sign: Sign::Negative,
        } => {
            return (
                BigFloat::try_from_i64_exact(2, target_precision).expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Normal { .. } => {}
    }

    // For x < 0, reflect via 2 − erfc(|x|). The subtraction loses
    // at most one bit because erfc(|x|) ≤ 2.
    if matches!(x.sign(), Sign::Negative) {
        let abs_x = x.abs();
        let working_prec = target_precision.saturating_add(8);
        let (erfc_abs, _) = abs_x
            .erfc_round(working_prec, RoundingMode::NearestEven)
            .expect("precision >= 1");
        let two = BigFloat::try_from_i64_exact(2, working_prec).expect("precision >= 1");
        let (result, _) = two.sub(&erfc_abs, RoundingMode::NearestEven);
        let (rounded, status) = result
            .round_to_precision(target_precision, mode)
            .expect("precision >= 1");
        auto_raise(status);
        return (rounded, status);
    }

    let e_x = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };

    let threshold = asymptotic_threshold_exponent(target_precision);
    let result = if e_x >= threshold {
        // Asymptotic produces target precision directly.
        let working_prec = target_precision
            .saturating_add(64)
            .min(target_precision.saturating_add(512));
        erfc_asymptotic(x, working_prec)
    } else {
        // Small-x path: 1 − erf_maclaurin(x) at boosted working
        // precision so the cancellation preserves target bits.
        let working_prec = target_precision
            .saturating_add(64)
            .saturating_add(128)
            .min(target_precision.saturating_add(512));
        let erf_val = erf_maclaurin(x, working_prec);
        let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
        let (diff, _) = one.sub(&erf_val, RoundingMode::NearestEven);
        diff
    };

    let (rounded, status) = result
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    auto_raise(status);
    (rounded, status)
}

/// `erfc(x)` for positive `x` past the asymptotic threshold, via
/// the divergent expansion
///
/// ```text
/// erfc(x) = (e^(−x²) / (x √π)) · Σ_{k=0}^N (−1)^k · (2k−1)!! / (2x²)^k.
/// ```
///
/// Truncation at the term of smallest magnitude gives an absolute
/// error bounded by that term, so for `|x|` large enough the series
/// dwarfs target precision. Caller guarantees the threshold check.
pub(super) fn erfc_asymptotic(x: &BigFloat, working_prec: u32) -> BigFloat {
    let x_w = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    let (x_sq, _) = x_w.mul(&x_w, RoundingMode::NearestEven);
    let neg_x_sq = x_sq.negated();
    let (exp_neg_x_sq, _) = neg_x_sq.exp(RoundingMode::NearestEven);

    // Leading envelope: e^(−x²) / (x · √π) = e^(−x²) · (2/√π) / (2x).
    let coef_two_sqrt_pi = two_over_sqrt_pi_at(working_prec);
    let two = BigFloat::try_from_i64_exact(2, working_prec).expect("precision >= 1");
    let (two_x, _) = two.mul(&x_w, RoundingMode::NearestEven);
    let (leading, _) = coef_two_sqrt_pi.div(&two_x, RoundingMode::NearestEven);
    let (envelope, _) = exp_neg_x_sq.mul(&leading, RoundingMode::NearestEven);

    // Series in t = 1 / (2x²): S = Σ (−1)^k · (2k − 1)!! · t^k.
    // Recurrence: t_k = t_{k−1} · (2k − 1) · t (magnitude); sign
    // alternates.
    let (two_x_sq, _) = two.mul(&x_sq, RoundingMode::NearestEven);
    let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let (t, _) = one.div(&two_x_sq, RoundingMode::NearestEven);

    let mut term = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let mut sum = term.clone();
    let mut prev_exp: i64 = 1;
    let mut sign_negative = true;

    let max_iter = working_prec.saturating_mul(2).max(256);
    for k in 1u32..=max_iter {
        let factor_int = i64::from(2 * k - 1);
        let factor_int_bf =
            BigFloat::try_from_i64_exact(factor_int, working_prec).expect("precision >= 1");
        let (step, _) = factor_int_bf.mul(&t, RoundingMode::NearestEven);
        let (next_term, _) = term.mul(&step, RoundingMode::NearestEven);
        let next_exp = match &next_term.class {
            Class::Normal { exponent, .. } => *exponent,
            _ => -i64::from(working_prec) - 1,
        };
        // Stop before the divergent tail starts growing again.
        if next_exp >= prev_exp {
            break;
        }
        prev_exp = next_exp;
        let signed = if sign_negative {
            next_term.negated()
        } else {
            next_term.clone()
        };
        sign_negative = !sign_negative;
        term = next_term;
        let (next_sum, _) = sum.add(&signed, RoundingMode::NearestEven);
        sum = next_sum;
        if next_exp < -i64::from(working_prec) - 4 {
            break;
        }
    }

    let (result, _) = envelope.mul(&sum, RoundingMode::NearestEven);
    result
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
    fn erfc_zero_is_one() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.erfc(RoundingMode::NearestEven);
            let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
            assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
        }
    }

    #[test]
    fn erfc_pos_inf_is_zero() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.erfc(RoundingMode::NearestEven);
        assert!(r.is_zero() && r.is_sign_positive());
    }

    #[test]
    fn erfc_neg_inf_is_two() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, _) = ni.erfc(RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert_eq!(r.partial_cmp(&two).0, Some(Ordering::Equal));
    }

    #[test]
    fn erfc_one() {
        // erfc(1) ≈ 0.15729920705028513
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.erfc(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "0.15729920705028513065877936491739",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 100));
    }

    #[test]
    fn erfc_two() {
        // erfc(2) ≈ 0.00467773498104726583793
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (r, _) = two.erfc(RoundingMode::NearestEven);
        let expected =
            BigFloat::parse_str("0.00467773498104726583793", 113, RoundingMode::NearestEven)
                .unwrap()
                .0;
        assert!(close_at(&r, &expected, 60));
    }

    #[test]
    fn erfc_six_is_tiny() {
        let six = BigFloat::try_from_i64_exact(6, 113).unwrap();
        let (r, _) = six.erfc(RoundingMode::NearestEven);
        let upper = BigFloat::parse_str("1e-15", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let lower = BigFloat::try_new_zero(Sign::Positive, 113).unwrap();
        assert_eq!(r.partial_cmp(&upper).0, Some(Ordering::Less));
        assert_eq!(r.partial_cmp(&lower).0, Some(Ordering::Greater));
    }

    #[test]
    fn erfc_negative_is_two_minus_positive() {
        let neg_one = BigFloat::try_from_i64_exact(-1, 113).unwrap();
        let (r, _) = neg_one.erfc(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "1.84270079294971486934122063508260925",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 100));
    }

    #[test]
    fn erf_erfc_sum_is_one() {
        let x = BigFloat::parse_str("0.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (a, _) = x.erf(RoundingMode::NearestEven);
        let (b, _) = x.erfc(RoundingMode::NearestEven);
        let (sum, _) = a.add(&b, RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        assert!(close_at(&sum, &one, 100));
    }

    #[test]
    fn erfc_large_x() {
        // erfc(20) ≈ 5.39e-176. Use the asymptotic and verify the
        // result is positive and finite with the expected order of
        // magnitude.
        let twenty = BigFloat::try_from_i64_exact(20, 113).unwrap();
        let (r, _) = twenty.erfc(RoundingMode::NearestEven);
        let upper = BigFloat::parse_str("1e-170", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let lower = BigFloat::parse_str("1e-200", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        assert_eq!(r.partial_cmp(&upper).0, Some(Ordering::Less));
        assert_eq!(r.partial_cmp(&lower).0, Some(Ordering::Greater));
    }

    #[test]
    fn erfc_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.erfc(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn erfc_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.erfc(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
