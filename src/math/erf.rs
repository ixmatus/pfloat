//! `erf(x) = (2/√π) · ∫₀ˣ e^(−t²) dt`: the error function.
//!
//! Algorithm: a single shared Maclaurin evaluator
//! [`erf_maclaurin`] computes
//!
//! ```text
//! erf(x) = (2/√π) · (x − x³/3 + x⁵/(2!·5) − x⁷/(3!·7) + …)
//! ```
//!
//! at working precision boosted by approximately `x² · log₂ e` so
//! the peak term (near `n ≈ x²`) does not exhaust the precision
//! budget. For `|x|` past a target-dependent threshold the kernel
//! switches to `1 − erfc_asymptotic(|x|)` to avoid the slow tail
//! of the Maclaurin.
//!
//! `erf_maclaurin` and [`super::erfc::erfc_asymptotic`] are both
//! `pub(super)` so the two kernels can share code without
//! recursing into each other's `kernel` functions.
//!
//! Special cases per IEEE 754-2019 §9.4.3:
//!
//! - `erf(±0) = ±0`.
//! - `erf(±∞) = ±1`.
//! - `erf(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::erfc::erfc_asymptotic;
use super::two_over_sqrt_pi_at;

impl BigFloat {
    /// `erf(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn erf(&self, mode: RoundingMode) -> (Self, Status) {
        self.erf_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `erf(self)` with explicit result precision.
    pub fn erf_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(erf_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `erf(self)` for `FixedFloat`. Delegates to [`BigFloat::erf`].
    #[must_use]
    pub fn erf(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().erf(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn erf_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            // erf(±∞) = ±1.
            let one = BigFloat::try_from_i64_exact(
                if matches!(sign, Sign::Negative) {
                    -1
                } else {
                    1
                },
                target_precision,
            )
            .expect("precision >= 1");
            return (one, Status::OK);
        }
        Class::Normal { .. } => {}
    }

    let sign = x.sign();
    let abs_x = x.abs();
    let e_x = match &abs_x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };

    let threshold = asymptotic_threshold_exponent(target_precision);
    let result_abs = if e_x >= threshold {
        // erf(|x|) = 1 − erfc_asymptotic(|x|) at working precision
        // generous enough to absorb the 1 − tiny cancellation.
        let working_prec = target_precision
            .saturating_add(64)
            .min(target_precision.saturating_add(512));
        let erfc_val = erfc_asymptotic(&abs_x, working_prec);
        let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
        let (diff, _) = one.sub(&erfc_val, RoundingMode::NearestEven);
        diff
    } else {
        erf_maclaurin(&abs_x, target_precision)
    };
    let result = if matches!(sign, Sign::Negative) {
        result_abs.negated()
    } else {
        result_abs
    };
    let (rounded, status) = result
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    auto_raise(status);
    (rounded, status)
}

/// Returns the smallest binary exponent at which the `erfc`
/// asymptotic gives at least `target_precision + 32` bits of
/// accuracy. For `|x|` with exponent `e_x ≥ threshold`, the
/// asymptotic dominates; below it, the Maclaurin is the only
/// path that converges.
///
/// Derivation: the asymptotic's smallest-term truncation gives
/// roughly `2^(−x² · log₂ e)` relative error. Setting
/// `x² · log₂ e ≥ p + 32` and using `|x|² ≥ 4^e_x` yields
/// `4^e_x ≥ (p + 32) / log₂ e`. Solving in integer arithmetic
/// (so the helper stays no_std-clean): use `log₂ e ≈ 23/16` and
/// pick the smallest `e` with `4^e ≥ ⌈(p + 32) · 16 / 23⌉`.
pub(super) fn asymptotic_threshold_exponent(target_precision: u32) -> i64 {
    let bits_needed: u64 = u64::from(target_precision) + 32;
    // need = ⌈bits_needed · 16 / 23⌉
    let need: u64 = (bits_needed * 16).div_ceil(23);
    let mut e: i64 = 2;
    let mut pow_4: u64 = 16; // 4^2
    while pow_4 < need && e < 60 {
        e += 1;
        pow_4 = pow_4.saturating_mul(4);
    }
    e
}

/// Maclaurin series for `erf(|x|)`. Always returns a non-negative
/// result. Callers apply the sign of the original `x`. The
/// working precision is boosted by approximately `x² · log₂ e`
/// bits so the peak term at `n ≈ x²` does not exhaust the budget.
pub(super) fn erf_maclaurin(abs_x: &BigFloat, target_precision: u32) -> BigFloat {
    if abs_x.is_zero() {
        return BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1");
    }

    let e_x = match &abs_x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };

    // Maclaurin's peak term sits near `n ≈ x²`. To absorb it we
    // need extra bits proportional to `x² · log₂ e`. Estimate
    // `x² ≤ 4 · 4^e_x` (the upper bound of `|x| ∈ [2^e, 2^(e+1))`)
    // and use the rational approximation `log₂ e ≈ 23/16` to stay
    // no_std-clean.
    let extra = if e_x <= 0 {
        0
    } else {
        let shift = (2 * (e_x + 1)).min(20) as u32;
        let mag: u64 = 1u64 << shift;
        (mag.saturating_mul(23) / 16).min(4096) as u32
    };
    let working_prec = target_precision
        .saturating_add(64)
        .saturating_add(extra)
        .min(target_precision.saturating_add(4096));

    let x_w = abs_x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    let (x_sq, _) = x_w.mul(&x_w, RoundingMode::NearestEven);
    let mut x_power = x_w.clone(); // x^(2n+1) starting at n=0 with x^1
    let mut factorial = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    // First term (n=0): x / 1.
    let mut sum = x_w.clone();

    let max_iter = working_prec.saturating_mul(4).max(2048);
    let mut sign_negative_term = true;
    for n in 1u32..=max_iter {
        let (next_power, _) = x_power.mul(&x_sq, RoundingMode::NearestEven);
        x_power = next_power;
        let n_big =
            BigFloat::try_from_i64_exact(i64::from(n), working_prec).expect("precision >= 1");
        let (next_fact, _) = factorial.mul(&n_big, RoundingMode::NearestEven);
        factorial = next_fact;
        let two_n_plus_one = BigFloat::try_from_i64_exact(i64::from(2 * n + 1), working_prec)
            .expect("precision >= 1");
        let (denom, _) = factorial.mul(&two_n_plus_one, RoundingMode::NearestEven);
        let (mut term, _) = x_power.div(&denom, RoundingMode::NearestEven);
        if sign_negative_term {
            term = term.negated();
        }
        sign_negative_term = !sign_negative_term;
        let (next_sum, _) = sum.add(&term, RoundingMode::NearestEven);
        sum = next_sum;

        if let Class::Normal { exponent, .. } = &term.class {
            if *exponent < -i64::from(working_prec) - 4 {
                break;
            }
        } else {
            break;
        }
    }

    let coef = two_over_sqrt_pi_at(working_prec);
    let (result, _) = coef.mul(&sum, RoundingMode::NearestEven);
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
    fn erf_zero() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.erf(RoundingMode::NearestEven);
            assert!(r.is_zero());
            assert_eq!(r.is_sign_negative(), matches!(s, Sign::Negative));
        }
    }

    #[test]
    fn erf_pos_inf_is_one() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.erf(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    #[test]
    fn erf_neg_inf_is_neg_one() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, _) = ni.erf(RoundingMode::NearestEven);
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        assert_eq!(r.partial_cmp(&neg_one).0, Some(Ordering::Equal));
    }

    #[test]
    fn erf_is_odd() {
        let half = BigFloat::parse_str("0.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let neg_half = half.negated();
        let (a, _) = half.erf(RoundingMode::NearestEven);
        let (b, _) = neg_half.erf(RoundingMode::NearestEven);
        let neg_a = a.negated();
        assert!(close_at(&b, &neg_a, 113 - 12));
    }

    #[test]
    fn erf_one() {
        // erf(1) ≈ 0.84270079294971486934122063508260925
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.erf(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "0.84270079294971486934122063508260925",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 100));
    }

    #[test]
    fn erf_two() {
        // erf(2) ≈ 0.99532226501895273416206925636725
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (r, _) = two.erf(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "0.99532226501895273416206925636725",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 100));
    }

    #[test]
    fn erf_six_is_basically_one() {
        let six = BigFloat::try_from_i64_exact(6, 113).unwrap();
        let (r, _) = six.erf(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Less));
        let (diff, _) = one.sub(&r, RoundingMode::NearestEven);
        let bound = BigFloat::parse_str("0.0000000000000001", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        assert_eq!(diff.partial_cmp(&bound).0, Some(Ordering::Less));
    }

    #[test]
    fn erf_large_x_saturates_to_one() {
        // For |x| past the asymptotic threshold at p=113, erf rounds
        // to ±1 exactly.
        let big = BigFloat::try_from_i64_exact(30, 113).unwrap();
        let (r, _) = big.erf(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    #[test]
    fn erf_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.erf(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn erf_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.erf(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixed_erf() {
        let one = FixedFloat::<113>::try_from_i64_exact(1).unwrap();
        let (r, _) = one.erf(RoundingMode::NearestEven);
        let one_fixed = FixedFloat::<113>::try_from_i64_exact(1).unwrap();
        assert_eq!(r.partial_cmp(&one_fixed).0, Some(Ordering::Less));
    }
}
