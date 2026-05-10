//! `digamma(x) = ψ(x) = Γ'(x)/Γ(x) = d/dx ln Γ(x)`: the digamma
//! function (a.k.a. psi). Defined for all real `x` except the
//! non-positive integers (where `Γ` has simple poles, so
//! `lnΓ` has logarithmic singularities and `ψ` has simple poles).
//!
//! Algorithm: differentiate the lgamma kernel.
//!
//! - For `x ≤ 0` integer: `−∞ + DIV_BY_ZERO` (pole).
//! - For `x < 0` non-integer: reflection
//!   `ψ(1 − x) − ψ(x) = π · cot(πx)`, i.e.,
//!   `ψ(x) = ψ(1 − x) − π · cot(πx)`.
//! - For `0 < x < z_min`: shift via the recurrence
//!   `ψ(z+1) = ψ(z) + 1/z`, i.e.,
//!   `ψ(x) = ψ(x + n) − Σ_{k=0}^{n−1} 1/(x + k)`.
//! - For `x ≥ z_min`: direct Stirling-like asymptotic via
//!   [`super::gamma_stirling::stirling_digamma`].

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::gamma_stirling::stirling_digamma;
use super::pi_at;

impl BigFloat {
    /// `digamma(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn digamma(&self, mode: RoundingMode) -> (Self, Status) {
        self.digamma_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `digamma(self)` with explicit result precision.
    pub fn digamma_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(digamma_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `digamma(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::digamma`].
    #[must_use]
    pub fn digamma(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().digamma(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn digamma_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            // ψ(0) is −∞ (simple pole).
            let neg_inf = BigFloat::try_new_infinity(Sign::Negative, target_precision)
                .expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            return (neg_inf, Status::DIV_BY_ZERO);
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

    // Negative branch.
    if matches!(x.sign(), Sign::Negative) {
        if super::lgamma::is_integer_test(x) {
            let neg_inf = BigFloat::try_new_infinity(Sign::Negative, target_precision)
                .expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            return (neg_inf, Status::DIV_BY_ZERO);
        }
        // Reflection: ψ(x) = ψ(1 − x) − π·cot(πx).
        let working_prec = target_precision
            .saturating_add(64)
            .min(target_precision.saturating_add(512));
        let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
        let (y, _) = one.sub(x, RoundingMode::NearestEven);
        let pi = pi_at(working_prec);
        let (pi_x, _) = pi.mul(x, RoundingMode::NearestEven);
        let (sin_pi_x, _) = pi_x.sin(RoundingMode::NearestEven);
        let (cos_pi_x, _) = pi_x.cos(RoundingMode::NearestEven);
        let (cot_pi_x, _) = cos_pi_x.div(&sin_pi_x, RoundingMode::NearestEven);
        let (pi_cot, _) = pi.mul(&cot_pi_x, RoundingMode::NearestEven);
        let (psi_y, _) = y
            .digamma_round(working_prec, RoundingMode::NearestEven)
            .expect("precision >= 1");
        let (result, _) = psi_y.sub(&pi_cot, RoundingMode::NearestEven);
        let (rounded, status) = result
            .round_to_precision(target_precision, mode)
            .expect("precision >= 1");
        auto_raise(status);
        return (rounded, status);
    }

    // Positive branch: shift up if needed, then apply Stirling.
    let working_prec = target_precision
        .saturating_add(64)
        .min(target_precision.saturating_add(512));
    let x_w = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    let z_min = z_min_for_target(target_precision);
    let e_x = match &x_w.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    let approx_x: u32 = if e_x < 0 {
        0
    } else if e_x >= 30 {
        u32::MAX
    } else {
        1u32 << ((e_x + 1) as u32)
    };

    let result_full_prec = if approx_x >= z_min {
        stirling_digamma(&x_w, working_prec)
    } else {
        let shifts = z_min - approx_x;
        let n_big =
            BigFloat::try_from_i64_exact(i64::from(shifts), working_prec).expect("precision >= 1");
        let (shifted, _) = x_w.add(&n_big, RoundingMode::NearestEven);
        let psi_z = stirling_digamma(&shifted, working_prec);
        // Σ_{k=0}^{shifts-1} 1/(x + k).
        let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
        let mut sum_recip =
            BigFloat::try_new_zero(Sign::Positive, working_prec).expect("precision >= 1");
        for k in 0..shifts {
            let k_big =
                BigFloat::try_from_i64_exact(i64::from(k), working_prec).expect("precision >= 1");
            let (denom, _) = x_w.add(&k_big, RoundingMode::NearestEven);
            let (term, _) = one.div(&denom, RoundingMode::NearestEven);
            let (next_sum, _) = sum_recip.add(&term, RoundingMode::NearestEven);
            sum_recip = next_sum;
        }
        psi_z.sub(&sum_recip, RoundingMode::NearestEven).0
    };

    let (rounded, status) = result_full_prec
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    auto_raise(status);
    (rounded, status)
}

/// `z_min` for the digamma asymptotic. Same shape as the lgamma
/// helper but slightly tighter because digamma's tail starts at
/// `z^(−2)` (one power deeper than lgamma) and the truncated
/// remainder is bounded by `|c_17·(2·17−1)| · z^(−34)`.
fn z_min_for_target(target_precision: u32) -> u32 {
    // Match the lgamma sizing — the difference between a `z^(−33)`
    // and `z^(−34)` tail is one extra bit of headroom, well
    // inside our `+32` margin.
    let log_z_needed = (target_precision + 60).div_ceil(33);
    let shift = log_z_needed.min(28);
    let z_min = 1u32 << shift;
    z_min.max(25)
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
    fn digamma_one_is_neg_euler() {
        // ψ(1) = -γ (Euler-Mascheroni) ≈ -0.57721566490153286.
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.digamma(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "-0.57721566490153286060651209008240243104215933593992",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn digamma_two_is_one_minus_euler() {
        // ψ(2) = 1 - γ ≈ 0.42278433509846714.
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (r, _) = two.digamma(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "0.42278433509846713939348790991759756895784066406008",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn digamma_half() {
        // ψ(1/2) = -γ - 2·ln(2) ≈ -1.9635100260214234794.
        let half = BigFloat::parse_str("0.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = half.digamma(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "-1.9635100260214234794409763329987555671931596046604",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn digamma_large_argument() {
        // ψ(100) ≈ 4.6001618527380874.
        let hundred = BigFloat::try_from_i64_exact(100, 113).unwrap();
        let (r, _) = hundred.digamma(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "4.6001618527380874001986055855758507268668127907685",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn digamma_recurrence() {
        // ψ(x + 1) = ψ(x) + 1/x.
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(7, p).unwrap();
        let x_plus_1 = BigFloat::try_from_i64_exact(8, p).unwrap();
        let (psi_x, _) = x.digamma(RoundingMode::NearestEven);
        let (psi_x_plus_1, _) = x_plus_1.digamma(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (recip, _) = one.div(&x, RoundingMode::NearestEven);
        let (lhs, _) = psi_x.add(&recip, RoundingMode::NearestEven);
        assert!(close_at(&lhs, &psi_x_plus_1, p - 20));
    }

    #[test]
    fn digamma_zero_is_neg_inf_div() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, status) = z.digamma(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative());
        assert!(status.div_by_zero());
    }

    #[test]
    fn digamma_negative_integer_is_neg_inf_div() {
        let neg_three = BigFloat::try_from_i64_exact(-3, 53).unwrap();
        let (r, status) = neg_three.digamma(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative());
        assert!(status.div_by_zero());
    }

    #[test]
    fn digamma_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.digamma(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn digamma_neg_inf_is_invalid() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, status) = ni.digamma(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn digamma_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.digamma(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn digamma_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.digamma(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
