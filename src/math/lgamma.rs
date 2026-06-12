//! `lgamma(x) = ln |Γ(x)|`: the log-magnitude of the gamma
//! function. Defined for all real `x` except the non-positive
//! integers (where `Γ` has simple poles).
//!
//! Algorithm:
//!
//! - For `x ≤ 0` and `x` a non-positive integer: `+∞ +
//!   DIV_BY_ZERO` (pole).
//! - For `x < 0` non-integer: reflection
//!   `ln|Γ(x)| = ln(π) − ln|sin(πx)| − lgamma(1 − x)`.
//! - For `0 < x < z_min(target_precision)`: shift via
//!   `Γ(z+1) = z · Γ(z)` to land at `z ≥ z_min`, accumulate
//!   `ln(x · (x+1) · … · (x+n−1))`, then `lgamma(x) =
//!   lgamma(x+n) − ln(product)`. `z_min` is sized so the truncated
//!   Stirling series clears `target_precision + 32` bits.
//! - For `x ≥ z_min`: direct Stirling.
//!
//! Correctly rounded under every IEEE rounding mode via the shared
//! [`crate::math::ziv::ziv_round`] driver (slice p1.2, ADR-0022).
//! The Stirling-plus-product-shift composition for small positive
//! `x` and the reflection composition for negative non-integer `x`
//! both run at the working precision the Ziv driver supplies; the
//! Ziv interval test handles the composition error analytically.
//! The pre-Ziv `+512` cap on the working-precision guard is dropped:
//! the Ziv loop now owns guard growth (64, 128, 256, 512, 1024).

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

use super::gamma_stirling::{spouge_lgamma_scaled, stirling_lgamma};
use super::pi_at;
use super::ziv::ziv_round;
use super::ziv_calibration::LGAMMA_ERROR_GUARD;

impl BigFloat {
    /// `lgamma(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn lgamma(&self, mode: RoundingMode) -> (Self, Status) {
        self.lgamma_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `lgamma(self)` with explicit result precision.
    pub fn lgamma_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(lgamma_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `lgamma(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::lgamma`].
    #[must_use]
    pub fn lgamma(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().lgamma(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn lgamma_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            // lgamma(±0) = +∞ + DIV_BY_ZERO (pole at 0).
            let inf = BigFloat::try_new_infinity(Sign::Positive, target_precision)
                .expect("precision >= 1");
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
            // lgamma(−∞) = +∞ (the natural extension of the
            // reflection: |Γ(−∞)| grows without bound in absolute
            // value along non-integer paths).
            return (
                BigFloat::try_new_infinity(Sign::Positive, target_precision)
                    .expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Normal { .. } => {}
    }

    // Negative integer pole, before any working-precision work.
    if matches!(x.sign(), Sign::Negative) && is_integer(x) {
        let inf =
            BigFloat::try_new_infinity(Sign::Positive, target_precision).expect("precision >= 1");
        auto_raise(Status::DIV_BY_ZERO);
        return (inf, Status::DIV_BY_ZERO);
    }

    // Exactly-representable true value dispatch (pf-kk16, ADR-0039).
    // lgamma(1) = ln(Γ(1)) = ln(1) = 0 and lgamma(2) = ln(Γ(2)) =
    // ln(1!) = 0 are both exactly representable at any precision.
    // For n ≥ 3, lgamma(n) = ln((n−1)!) is irrational (the
    // factorial is not a power of two), so only n ∈ {1, 2} is in
    // the exact-value dispatch subset. The Stirling/Spouge
    // composition at x ∈ {1, 2} would return 0 + epsilon under
    // directed modes and tip rounding to the smallest representable
    // value away from zero (the gamma(7) defect-class shape
    // recorded at `feedback_exact_value_defeats_ziv`).
    if let Some(exact) = try_lgamma_small_pos_int_exact(x, target_precision) {
        return (exact, Status::OK);
    }

    // Negative non-integer and positive finite both feed the Ziv
    // driver. The closure dispatches on the sign at each retry.
    let z_min = z_min_for_target(target_precision);
    let (result, status) = ziv_round(
        |working_prec| lgamma_at_w(x, z_min, working_prec),
        target_precision,
        mode,
        LGAMMA_ERROR_GUARD,
    );
    // Defensive INEXACT guard (pf-umlm, ADR-0066): a finite-normal
    // fall-through (x ∉ {1, 2}, which are dispatched above) is ln|Γ(x)|,
    // irrational. The ADR-0065 sweep showed this path already flags
    // INEXACT everywhere, so the force is a no-op hardening against
    // regression; its worst-case soundness rests on the irrationality of
    // ln|Γ| at dyadic arguments, not proven for every dyadic.
    let status = super::force_transcendental_inexact(&result, status);
    auto_raise(status);
    (result, status)
}

/// If `x ∈ {1, 2}` then `lgamma(x) = 0` exactly at any target
/// precision; otherwise returns `None`. lgamma at other positive
/// integers is `ln((n−1)!)`, irrational because the factorial is
/// not a power of two — those route through the Ziv envelope.
fn try_lgamma_small_pos_int_exact(x: &BigFloat, target_precision: u32) -> Option<BigFloat> {
    if !matches!(x.sign(), Sign::Positive) || x.is_zero() {
        return None;
    }
    let one = BigFloat::try_from_i64_exact(1, target_precision).ok()?;
    let two = BigFloat::try_from_i64_exact(2, target_precision).ok()?;
    if matches!(x.partial_cmp(&one).0, Some(Ordering::Equal))
        || matches!(x.partial_cmp(&two).0, Some(Ordering::Equal))
    {
        return BigFloat::try_new_zero(Sign::Positive, target_precision).ok();
    }
    None
}

/// Evaluate `lgamma(x)` at the supplied working precision via the
/// reflection formula (for negative non-integer `x`) or the
/// Stirling-with-upward-shift composition (for positive `x`). The
/// caller's special-case handling has peeled off NaN, ±0, ±∞, and
/// the negative-integer pole. Returns the unrounded value; the Ziv
/// driver handles rounding to the caller's target precision and
/// mode.
fn lgamma_at_w(x: &BigFloat, z_min: u32, working_prec: u32) -> BigFloat {
    if matches!(x.sign(), Sign::Negative) {
        // ln|Γ| has roots on the negative axis where the reflection
        // ln(π) − ln|sin(πx)| − lgamma(1−x) is a near-total cancellation
        // of O(1) terms; boost the working precision by the realised
        // cancellation so the Ziv half-width stays sound (review
        // 2026-05-29, root cause 2).
        //
        // Separately, proximity to the negative-axis POLES (the
        // integers) is input-encoded and collapses inside π·x before
        // sin ever sees it: at x = −3 + 2^-k the product carries the
        // 2^-k offset only to the working width, sin's relative error
        // grows by 2^k, and ln amplifies it into the result while the
        // result itself stays O(k) — so the realised-cancellation
        // probe never fires (lgamma(−3+2^-80 @p84) → 53 certified a
        // value wrong from bit ~41; pf-pdda's deep-beta consumer hit
        // this through the exact-sum handoff, ADR-0098; found by the
        // slice's adversarial verification, pre-existing). The depth
        // is exactly computable up front (the gap to the nearest
        // integer is on x's own grid): pre-boost by it, the
        // asin/ln pattern of ADR-0097.
        let pole_boost = pole_proximity_depth(x).saturating_add(8);
        let working_prec = working_prec.saturating_add(pole_boost);
        return super::ziv::cancellation_boosted(working_prec, |w| {
            let one = BigFloat::try_from_i64_exact(1, w).expect("precision >= 1");
            // y = 1 − x, positive since x < 0.
            let (y, _) = one.sub(x, RoundingMode::NearestEven);
            let pi = pi_at(w);
            let (pi_x, _) = pi.mul(x, RoundingMode::NearestEven);
            let (sin_val, _) = pi_x.sin(RoundingMode::NearestEven);
            let abs_sin = sin_val.abs();
            let (ln_sin, _) = abs_sin.ln(RoundingMode::NearestEven);
            let (ln_pi, _) = pi.ln(RoundingMode::NearestEven);
            let (lgamma_y, _) = y
                .lgamma_round(w, RoundingMode::NearestEven)
                .expect("precision >= 1");
            let (mid, _) = ln_pi.sub(&ln_sin, RoundingMode::NearestEven);
            let (result, _) = mid.sub(&lgamma_y, RoundingMode::NearestEven);
            let op_scale =
                super::ziv::value_exponent(&mid).max(super::ziv::value_exponent(&lgamma_y));
            (result, op_scale)
        });
    }

    // Positive branch. Near the positive roots x = 1 and x = 2 the
    // composition (Stirling-with-shift's lgamma(z) − ln(product), or
    // Spouge's leading-minus-ln chain) is a near-total cancellation
    // of terms whose scale grows with z_min, while the result shrinks
    // with the input's proximity to the root — proximity the relative
    // half-width model cannot see. lgamma(2 + 2^-100) certified a
    // value with relative error 2.5e-3 (pf-wmv7, ADR-0097). Mirror
    // the negative branch: inside the root windows, boost by the
    // realised cancellation. The window [3/4, 5/4] ∪ [7/4, 9/4]
    // bounds |lgamma| below ~2^-3.4 outside it, so the un-boosted
    // path's cancellation stays inside the Ziv guard there; inside,
    // z_min must be re-derived from the boosted precision (a z_min
    // sized for the original target caps Stirling's truncation
    // accuracy no matter how far the working precision grows).
    if in_positive_root_window(x) {
        return super::ziv::cancellation_boosted(working_prec, |w| {
            lgamma_positive_at_w(x, z_min_for_target(w), w)
        });
    }
    lgamma_positive_at_w(x, z_min, working_prec).0
}

/// `x ∈ [3/4, 5/4] ∪ [7/4, 9/4]`: the windows around lgamma's
/// positive roots at 1 and 2. Exact dyadic bounds compared on the
/// original input, so the trigger is precision-independent.
fn in_positive_root_window(x: &BigFloat) -> bool {
    let quarter = |n: i64| {
        BigFloat::try_from_i64_exact(n, 4)
            .expect("4 bits hold 3..=9")
            .scale_by_pow2(-2)
            .0
    };
    let within = |lo: i64, hi: i64| {
        matches!(
            x.partial_cmp(&quarter(lo)).0,
            Some(Ordering::Greater | Ordering::Equal)
        ) && matches!(
            x.partial_cmp(&quarter(hi)).0,
            Some(Ordering::Less | Ordering::Equal)
        )
    };
    within(3, 5) || within(7, 9)
}

/// The positive-branch evaluation at one working precision,
/// returning `(value, operand_scale)` — the scale of the largest
/// term that cancelled to form the value, which is what
/// `cancellation_boosted` charges the realised cancellation
/// against. Callers outside the root windows discard the scale.
fn lgamma_positive_at_w(x: &BigFloat, z_min: u32, working_prec: u32) -> (BigFloat, i64) {
    let x_w = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    // Dispatch on target precision. The 17-Bernoulli-pair Stirling
    // table caps accuracy at ~895 bits regardless of shift (pf-l6s5);
    // for working precisions past that reach, Spouge's approximation
    // delivers the full target accuracy with cost linear in its
    // parameter `a`. The 600-bit threshold leaves Stirling room for
    // its `+32`-bit margin from z_min_for_target and the Ziv guard.
    if working_prec > STIRLING_REACH_THRESHOLD {
        return spouge_lgamma_scaled(&x_w, working_prec);
    }

    // Decide whether to shift. We want z = x + n ≥ z_min.
    let e_x = match &x_w.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    let approx_x: u32 = if e_x < 0 {
        0
    } else if e_x >= 30 {
        // x ≥ 2^30 already past any z_min; no shift.
        u32::MAX
    } else {
        // x ∈ [2^e_x, 2^(e_x+1)); take 2^(e_x+1) as an upper-bound
        // estimate sufficient for "should we shift?" only.
        1u32 << ((e_x + 1) as u32)
    };

    if approx_x >= z_min {
        let v = stirling_lgamma(&x_w, working_prec);
        let scale = super::ziv::value_exponent(&v);
        (v, scale)
    } else {
        // Number of upward shifts so x + n ≥ z_min.
        let shifts = z_min - approx_x;
        let shifted = shift_up(&x_w, shifts, working_prec);
        let lgamma_z = stirling_lgamma(&shifted, working_prec);
        let ln_prod = product_ln(&x_w, shifts, working_prec);
        let (diff, _) = lgamma_z.sub(&ln_prod, RoundingMode::NearestEven);
        let scale = super::ziv::value_exponent(&lgamma_z).max(super::ziv::value_exponent(&ln_prod));
        (diff, scale)
    }
}

/// Working-precision threshold above which the lgamma kernel
/// dispatches to Spouge's approximation. Below the threshold, the
/// 17-pair Stirling table delivers the required precision via the
/// upward-shift composition. The threshold is set conservatively
/// (well below the ~895-bit Stirling reach cap, pf-l6s5) to keep
/// the `differential_gamma` lane on the existing Stirling path at
/// every `TRANSCENDENTAL_PRECISIONS` entry.
const STIRLING_REACH_THRESHOLD: u32 = 600;

/// Picks the shift target `z_min` such that Stirling truncated at
/// the 17 hardcoded coefficients clears `target_precision + 32`
/// bits.
///
/// Derivation: the truncation error after 17 terms is bounded by
/// `|c_17| · z^(−33)`. With `|c_17| ≈ 2^28.5` (from the hardcoded
/// table), `log₂(z) ≥ (target + 60) / 33` makes the residual fall
/// below `2^(−target − 32)`. We round the required `log₂(z)` up
/// and exponentiate; the result caps at `2^28` so the shift count
/// stays bounded.
///
/// Targets past `~600` bits exceed the table's reach; the gamma
/// kernels still produce a value but with degraded accuracy.
fn z_min_for_target(target_precision: u32) -> u32 {
    let log_z_needed = (target_precision + 60).div_ceil(33);
    let shift = log_z_needed.min(28);
    let z_min = 1u32 << shift;
    z_min.max(25)
}

/// Returns `x + n` at `working_prec`.
fn shift_up(x: &BigFloat, n: u32, working_prec: u32) -> BigFloat {
    let n_big = BigFloat::try_from_i64_exact(i64::from(n), working_prec).expect("precision >= 1");
    x.add(&n_big, RoundingMode::NearestEven).0
}

/// Returns `ln(x · (x+1) · … · (x + n − 1))` at `working_prec`.
fn product_ln(x: &BigFloat, n: u32, working_prec: u32) -> BigFloat {
    let mut product = x.clone();
    for k in 1..n {
        let k_big =
            BigFloat::try_from_i64_exact(i64::from(k), working_prec).expect("precision >= 1");
        let (factor, _) = x.add(&k_big, RoundingMode::NearestEven);
        let (next, _) = product.mul(&factor, RoundingMode::NearestEven);
        product = next;
    }
    product.ln(RoundingMode::NearestEven).0
}

/// Bits of proximity from `x` to its nearest integer: `−exponent`
/// of the (exactly computed) gap, or 0 when `x` is at least 2^-1
/// away or sits exactly on an integer (the caller dispatched the
/// poles already). The gap lives on `x`'s own grid, so the
/// subtraction is exact; for `|x| < 1` both neighbouring integers
/// (0 and ±1) are checked. `pub(super)`: digamma's reflection has
/// the same pole structure (ADR-0098).
pub(super) fn pole_proximity_depth(x: &BigFloat) -> u32 {
    let e_x = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => return 0,
    };
    let gap_exponent = if e_x < 0 {
        // |x| < 1: nearest integers are 0 (gap |x|, exponent e_x)
        // and ±1 (gap 1 − |x|, Sterbenz-exact).
        let one = BigFloat::try_from_i64_exact(1, x.precision()).expect("precision >= 1");
        let (gap_to_one, _) = one.sub(&x.abs(), RoundingMode::NearestEven);
        match &gap_to_one.class {
            Class::Normal { exponent, .. } => e_x.min(*exponent),
            // x = ±1 exactly: a pole/root, dispatched upstream.
            _ => return 0,
        }
    } else {
        if e_x >= i64::from(x.precision()) {
            // ulp(x) ≥ 1: x is itself an integer (dispatched).
            return 0;
        }
        // Nearest integer: round x to e_x + 1 mantissa bits (ulp 1).
        let int_bits = u32::try_from(e_x + 1).expect("e_x < precision <= u32::MAX");
        let n = x
            .round_to_precision(int_bits, RoundingMode::NearestEven)
            .expect("precision >= 1")
            .0;
        let (gap, _) = x.sub(&n, RoundingMode::NearestEven);
        match &gap.class {
            Class::Normal { exponent, .. } => *exponent,
            _ => return 0,
        }
    };
    if gap_exponent >= -1 {
        return 0;
    }
    u32::try_from(gap_exponent.saturating_neg()).unwrap_or(u32::MAX)
}

/// Returns `true` if `x` is a finite integer.
fn is_integer(x: &BigFloat) -> bool {
    match &x.class {
        Class::Zero { .. } => true,
        Class::Normal { exponent, .. } => {
            let scale = *exponent - i64::from(x.precision) + 1;
            if scale >= 0 {
                return true;
            }
            // |x| < 1 cannot be a nonzero integer.
            if *exponent < 0 {
                return false;
            }
            // For scale < 0, check that the low |scale| bits of the
            // mantissa-as-integer are zero.
            let m_int = crate::ops::limbs::extract_as_integer(
                match &x.class {
                    Class::Normal { mantissa, .. } => mantissa,
                    _ => unreachable!(),
                },
                x.precision,
            );
            let abs_scale = (-scale) as u32;
            let full_limbs = (abs_scale / 64) as usize;
            let partial_bits = abs_scale % 64;
            for &limb in m_int.iter().take(full_limbs) {
                if limb != 0 {
                    return false;
                }
            }
            if partial_bits > 0 {
                let mask = (1u64 << partial_bits) - 1;
                if m_int.get(full_limbs).copied().unwrap_or(0) & mask != 0 {
                    return false;
                }
            }
            true
        }
        _ => false,
    }
}

// Re-export for testing in gamma.rs.
#[allow(dead_code)]
pub(super) fn is_integer_test(x: &BigFloat) -> bool {
    is_integer(x)
}

// Silence unused-import warning when `Ordering` isn't directly used.
#[allow(dead_code)]
const _USE_ORDERING: Option<Ordering> = None;

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
    fn lgamma_one_is_zero() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.lgamma(RoundingMode::NearestEven);
        assert!(
            r.is_zero()
                || close_at(
                    &r,
                    &BigFloat::try_new_zero(Sign::Positive, 113).unwrap(),
                    100
                )
        );
    }

    #[test]
    fn lgamma_one_two_is_zero_under_every_directed_mode() {
        // pf-kk16 pinning test: the exact-value pre-Ziv dispatch for
        // lgamma(1) and lgamma(2) returns exactly +0 under every
        // mode. Without the dispatch, the Stirling/Spouge composition
        // returns 0 + epsilon and TP would round up to the smallest
        // representable positive value (the gamma(7) defect shape
        // recorded at feedback_exact_value_defeats_ziv).
        for &mode in &[
            RoundingMode::NearestEven,
            RoundingMode::TowardPositive,
            RoundingMode::TowardNegative,
            RoundingMode::TowardZero,
            RoundingMode::NearestAway,
        ] {
            for &prec in &[24u32, 53, 113] {
                for &n in &[1i64, 2] {
                    let x = BigFloat::try_from_i64_exact(n, prec).unwrap();
                    let (r, status) = x.lgamma(mode);
                    assert!(status.is_ok(), "lgamma({n}) status under {mode:?}@p{prec}");
                    assert!(
                        r.is_zero() && !r.is_sign_negative(),
                        "lgamma({n}) should be +0 under {mode:?}@p{prec}, got {r:?}"
                    );
                    assert_eq!(r.precision(), prec);
                }
            }
        }
    }

    #[test]
    fn lgamma_two_is_zero() {
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (r, _) = two.lgamma(RoundingMode::NearestEven);
        let zero = BigFloat::try_new_zero(Sign::Positive, 113).unwrap();
        assert!(close_at(&r, &zero, 100));
    }

    #[test]
    fn lgamma_factorials() {
        // lgamma(6) = ln(5!) = ln(120) ≈ 4.78749174278204599.
        let six = BigFloat::try_from_i64_exact(6, 113).unwrap();
        let (r, _) = six.lgamma(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "4.787491742782045994247700934523243048399592315172",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 90));
    }

    #[test]
    fn lgamma_half_is_half_ln_pi() {
        // lgamma(1/2) = ln(√π) = (1/2) · ln(π) ≈ 0.5723649429247.
        let half = BigFloat::parse_str("0.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = half.lgamma(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "0.57236494292470008707171367567652935582364740645766",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn lgamma_zero_is_pos_inf_div() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, status) = z.lgamma(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
        assert!(status.div_by_zero());
    }

    #[test]
    fn lgamma_negative_integer_is_pos_inf_div() {
        let neg_three = BigFloat::try_from_i64_exact(-3, 53).unwrap();
        let (r, status) = neg_three.lgamma(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
        assert!(status.div_by_zero());
    }

    #[test]
    fn lgamma_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.lgamma(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn lgamma_neg_inf() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, _) = ni.lgamma(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn lgamma_negative_non_integer_via_reflection() {
        // lgamma(-0.5) = ln(|Γ(-0.5)|) = ln(2√π) ≈ 1.2655121234846454.
        let neg_half = BigFloat::parse_str("-0.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = neg_half.lgamma(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "1.2655121234846453964889457971347059238991475408179",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn lgamma_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.lgamma(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn lgamma_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.lgamma(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
