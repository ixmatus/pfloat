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

use super::gamma_stirling::stirling_digamma;
use super::pi_at;
use super::ziv::ziv_round;
use super::ziv_calibration::DIGAMMA_ERROR_GUARD;

impl BigFloat {
    /// `digamma(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn digamma(&self, mode: RoundingMode) -> (Self, Status) {
        self.digamma_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `digamma(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.29, ADR-0038).
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

    // Negative integer pole, before any working-precision work.
    if matches!(x.sign(), Sign::Negative) && super::lgamma::is_integer_test(x) {
        let neg_inf =
            BigFloat::try_new_infinity(Sign::Negative, target_precision).expect("precision >= 1");
        auto_raise(Status::DIV_BY_ZERO);
        return (neg_inf, Status::DIV_BY_ZERO);
    }

    // z_min pinned from target_precision so the shift count does not
    // flip across Ziv retries (mirrors the lgamma precedent). The
    // regime dispatch (direct Stirling vs shift-then-Stirling) reads
    // from this pinned value inside eval(w).
    let z_min = z_min_for_target(target_precision);
    let (result, status) = ziv_round(
        |w| digamma_at_w(x, z_min, w),
        target_precision,
        mode,
        DIGAMMA_ERROR_GUARD,
    );
    // Defensive INEXACT guard (pf-umlm, ADR-0066): ψ has no dyadic
    // outputs (every value is ψ(n) = −γ + H_{n−1} or a reflection of one,
    // irrational), so a finite-normal result is INEXACT. The ADR-0065
    // sweep showed this path already flags it everywhere, so the force is
    // a no-op hardening against regression; its worst-case soundness
    // rests on the irrationality of γ = −ψ(1), an open problem.
    let status = super::force_transcendental_inexact(&result, status);
    auto_raise(status);
    (result, status)
}

/// Evaluate `digamma(x)` at the supplied working precision under
/// `NearestEven`. The caller has peeled off NaN, ±0, ±∞, and the
/// negative-integer pole. The reflection branch (negative non-integer
/// `x`) composes `ψ(1 − x) − π·cot(πx)` with a recursive
/// `digamma_round` call that routes through the positive branch
/// (since `1 − x > 1` when `x < 0`), so it does not recurse
/// indefinitely. The positive branch dispatches on `z_min` (pinned
/// by the caller from `target_precision`) between direct Stirling
/// and shift-then-Stirling per the recurrence
/// `ψ(x+1) = ψ(x) + 1/x`.
fn digamma_at_w(x: &BigFloat, z_min: u32, working_prec: u32) -> BigFloat {
    if matches!(x.sign(), Sign::Negative) {
        // ψ has roots on the negative axis where the reflection
        // ψ(1 − x) − π·cot(πx) is a near-total cancellation of O(1)
        // terms; boost the working precision by the realised
        // cancellation so the Ziv half-width stays sound (review
        // 2026-05-29, root cause 2). Proximity to the negative-axis
        // POLES additionally collapses inside π·x before sin/cos see
        // it (the lgamma reflection analog, ADR-0098): pre-boost by
        // the exactly-computed depth to the nearest integer.
        let pole_boost = super::lgamma::pole_proximity_depth(x).saturating_add(8);
        let working_prec = working_prec.saturating_add(pole_boost);
        return super::ziv::cancellation_boosted(working_prec, |w| {
            // Reflection: ψ(x) = ψ(1 − x) − π·cot(πx).
            let one = BigFloat::try_from_i64_exact(1, w).expect("precision >= 1");
            let (y, _) = one.sub(x, RoundingMode::NearestEven);
            let pi = pi_at(w);
            let (pi_x, _) = pi.mul(x, RoundingMode::NearestEven);
            let (sin_pi_x, _) = pi_x.sin(RoundingMode::NearestEven);
            let (cos_pi_x, _) = pi_x.cos(RoundingMode::NearestEven);
            let (cot_pi_x, _) = cos_pi_x.div(&sin_pi_x, RoundingMode::NearestEven);
            let (pi_cot, _) = pi.mul(&cot_pi_x, RoundingMode::NearestEven);
            let (psi_y, _) = y
                .digamma_round(w, RoundingMode::NearestEven)
                .expect("precision >= 1");
            let (result, _) = psi_y.sub(&pi_cot, RoundingMode::NearestEven);
            let op_scale =
                super::ziv::value_exponent(&psi_y).max(super::ziv::value_exponent(&pi_cot));
            (result, op_scale)
        });
    }

    // Positive branch. Near ψ's positive root (1.46163…) the shift
    // composition ψ(z) − Σ 1/(x+k) is a near-total cancellation
    // whose depth is the input's proximity to the irrational root —
    // bounded by the input precision and invisible to the relative
    // half-width model (pf-wmv7, ADR-0097): the p100 rounding of the
    // root certified garbage at target 53. Mirror the negative
    // branch: inside the window [5/4, 7/4] boost by the realised
    // cancellation, re-deriving z_min from the boosted precision (a
    // z_min sized for the original target caps the asymptotic's
    // truncation accuracy no matter the working precision). Outside
    // the window |ψ| ≥ |ψ(5/4)| ≈ 2^-2.1, inside the guard's
    // reach.
    // Two cancellations route the positive branch through
    // cancellation_boosted: the near-root window (above), and the whole
    // Spouge regime (working_prec > STIRLING_REACH_THRESHOLD), whose S/S'
    // sum cancellation grows with the argument (~0.4·w at z≈1e6) and so
    // is uncoverable by a fixed margin (pf-0r1l verifier; the
    // spouge_digamma_scaled scale reports the depth and the iteration
    // recovers it). Below the threshold the shift-Stirling path has no
    // sum cancellation and runs directly (the differential_digamma lane,
    // p ≤ 256, stays on this fast path).
    if in_positive_root_window(x) || working_prec > STIRLING_REACH_THRESHOLD {
        return super::ziv::cancellation_boosted(working_prec, |w| {
            digamma_positive_at_w(x, z_min_for_target(w), w)
        });
    }
    digamma_positive_at_w(x, z_min, working_prec).0
}

/// `x ∈ [5/4, 7/4]`: the window around ψ's positive root. Exact
/// dyadic bounds compared on the original input, so the trigger is
/// precision-independent.
fn in_positive_root_window(x: &BigFloat) -> bool {
    let quarter = |n: i64| {
        BigFloat::try_from_i64_exact(n, 3)
            .expect("3 bits hold 5 and 7")
            .scale_by_pow2(-2)
            .0
    };
    matches!(
        x.partial_cmp(&quarter(5)).0,
        Some(Ordering::Greater | Ordering::Equal)
    ) && matches!(
        x.partial_cmp(&quarter(7)).0,
        Some(Ordering::Less | Ordering::Equal)
    )
}

/// The positive-branch evaluation at one working precision,
/// returning `(value, operand_scale)` for `cancellation_boosted`;
/// callers outside the root window discard the scale.
fn digamma_positive_at_w(x: &BigFloat, z_min: u32, working_prec: u32) -> (BigFloat, i64) {
    let x_w = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    // Past the 17-pair Stirling table's reach, dispatch to the Spouge
    // derivative (pf-0r1l, ADR-0110): reaching the table for small x
    // requires shifting up to z_min, which the target/log2(z_min)
    // sizing caps at 2^28 — a ~268M-term shift sum costing ~28 minutes
    // for the deep-root inputs whose Ziv working precision climbs past
    // ~900 bits. Spouge's cost is linear in a ∝ working_prec, correct
    // and cost-proportional past the table reach. The threshold matches
    // lgamma's; below it the shift-Stirling path stays the faster one
    // and keeps the differential_digamma lane (p ≤ 256) on Stirling.
    if working_prec > STIRLING_REACH_THRESHOLD {
        return super::gamma_stirling::spouge_digamma_scaled(&x_w, working_prec);
    }

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

    if approx_x >= z_min {
        let v = stirling_digamma(&x_w, working_prec);
        let scale = super::ziv::value_exponent(&v);
        (v, scale)
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
        let (diff, _) = psi_z.sub(&sum_recip, RoundingMode::NearestEven);
        let scale = super::ziv::value_exponent(&psi_z).max(super::ziv::value_exponent(&sum_recip));
        (diff, scale)
    }
}

/// Working-precision threshold above which `digamma_positive_at_w`
/// dispatches to the Spouge derivative ([`super::gamma_stirling::spouge_digamma_scaled`]).
/// Below it, the 17-pair Stirling asymptotic with upward shift remains
/// the faster path. Set to lgamma's value (the shared
/// 17-Bernoulli-pair table caps both at ~895 bits) and conservatively
/// below it, so the `differential_digamma` lane (p ≤ 256) and the
/// shallow Ziv rungs stay on Stirling (pf-0r1l, ADR-0110).
const STIRLING_REACH_THRESHOLD: u32 = 600;

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
