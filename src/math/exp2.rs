//! `exp2(x) = 2^x`: binary exponential.
//!
//! Composition: `2^x = exp(x · ln(2))`. The kernel computes the
//! product at working precision (with `ln(2)` from the shared
//! 1024-bit constant), then calls `exp`. All special cases flow
//! through composition.
//!
//! Correctly rounded under every IEEE 754-2019 rounding mode via
//! the shared [`crate::math::ziv::ziv_round`] driver (slice p1.24,
//! ADR-0038). The `exp · ln(2)` composition has no cancellation
//! regime; the Ziv envelope's working-precision growth certifies
//! the rounding-mode interval test on the final round.
//!
//! Special cases per IEEE 754-2019 §9.2 reduce to:
//!
//! - `exp2(±0) = 1`.
//! - `exp2(+∞) = +∞`, `exp2(−∞) = +0`.
//! - `exp2(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::ln_2_at;
use super::ziv::ziv_round;
use super::ziv_calibration::EXP2_ERROR_GUARD;

impl BigFloat {
    /// `2^self` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn exp2(&self, mode: RoundingMode) -> (Self, Status) {
        self.exp2_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `exp2(self)` with explicit result precision.
    pub fn exp2_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(exp2_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `exp2(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::exp2`].
    #[must_use]
    pub fn exp2(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().exp2(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn exp2_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
                BigFloat::try_new_infinity(Sign::Positive, target_precision)
                    .expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Infinity {
            sign: Sign::Negative,
        } => {
            return (
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Normal { .. } => {}
    }

    // Exact-input dispatch (pf-njs5, ADR-0060). 2^x is exactly
    // representable iff x is an integer: then 2^x is a single-bit
    // power of two, exact at every precision and exponent. For a
    // non-integer x, 2^x is irrational — Gelfond–Schneider gives
    // transcendence for irrational algebraic exponents, and a
    // rational non-integer p/q yields an irrational q-th root — so
    // the composition fall-through below forces INEXACT.
    if let Some(k) = super::pow::integer_exponent(x) {
        return (two_pow_at(k, target_precision), Status::OK);
    }

    // Exponent-rim triage (pf-qm0h, the ADR-0096 pattern; found by
    // the R1-merge CI failure). The composition below discards exp's
    // Status, so exp's mode-aware rim dispatch arrived here as a bare
    // +inf/+0 that half_width(non-Normal) = 0 certified with INEXACT
    // only. exp2 needs none of exp's certified-division machinery:
    // the result's binary exponent IS floor(x), exactly computable.
    // e_x ≤ 61 keeps the established path untouched (|floor(x)| <
    // 2^62, no rim interaction).
    let e_x = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => unreachable!("specials handled above"),
    };
    if e_x >= 62 {
        return exp2_extreme(x, e_x, target_precision, mode);
    }

    // Ziv-driven correct rounding under every IEEE mode. The
    // composition `exp(x · ln(2))` has no cancellation regime; the
    // Ziv driver's working-precision growth handles the rounding-
    // mode interval test at the final round to target.
    let (result, status) = ziv_round(
        |w| {
            let x_w = x
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let ln_2 = ln_2_at(w);
            let (product, _) = x_w.mul(&ln_2, RoundingMode::NearestEven);
            let (e_val, _) = product.exp(RoundingMode::NearestEven);
            e_val
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0
        },
        target_precision,
        mode,
        EXP2_ERROR_GUARD,
    );
    // Non-integer x ⟹ 2^x irrational ⟹ INEXACT, even where the
    // working-precision evaluation rounds onto a grid value.
    let status = status | Status::INEXACT;
    auto_raise(status);
    (result, status)
}

/// `2^k` at precision `precision`: `1.0` with its exponent shifted by
/// `k`. `2^k` is a single-bit power of two, exact at every precision;
/// [`super::pow::integer_exponent`] guarantees `|k| < 2^63`, so the
/// exponent shift never saturates.
fn two_pow_at(k: i64, precision: u32) -> BigFloat {
    let mut v = BigFloat::try_from_i64_exact(1, precision).expect("precision >= 1");
    if let Class::Normal { exponent, .. } = &mut v.class {
        *exponent = exponent.saturating_add(k);
    }
    v
}

/// `2^x` for `|x| ≥ 2^62` (pf-qm0h, ADR-0096 pattern): the result's
/// binary exponent is `n = floor(x)`, exact integer arithmetic on the
/// input's own grid — no certified division. Classification against
/// the i64 rim, reusing exp's mode-aware dispatch results:
///
/// - positive `x`, `e_x ≥ 63`: `n ≥ 2^63`, certain overflow.
/// - negative `x`, `e_x ≥ 64`: `n < i64::MIN − 1`, deep underflow.
/// - the bands in between classify on the exact `floor(x)`:
///   overflow / sliver / deep underflow at the rims, else the
///   representable window: Ziv-certify the unscaled `s = 2^frac`
///   (`frac = x − n` exact, `s ∈ [1, 2)`) and compose with the exact
///   `scale_by_pow2(k)`; `k` is clamped one step inside the rim so a
///   certified carry routes to the overflow dispatch, the ADR-0096
///   compose contract.
fn exp2_extreme(
    x: &BigFloat,
    e_x: i64,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    use super::exp::{
        exp_certain_overflow, exp_deep_underflow, exp_underflow_sliver, round_bigfloat_to_i64,
    };

    let positive = !x.is_sign_negative();
    if (positive && e_x >= 63) || (!positive && e_x >= 64) {
        return if positive {
            exp_certain_overflow(target_precision, mode)
        } else {
            exp_deep_underflow(target_precision, mode)
        };
    }

    // The bands: positive e_x = 62 (floor in [2^62, 2^63), always a
    // representable window) and negative e_x ∈ {62, 63} (floor in
    // (−2^64, −2^62], straddling i64::MIN). floor(x) is exact: round
    // x to e_x + 1 mantissa bits toward −∞ (ulp = 1 there).
    let int_bits = u32::try_from(e_x + 1).expect("e_x ∈ [62, 63]");
    let (n_big, _) = x
        .round_to_precision(int_bits.max(1), RoundingMode::TowardNegative)
        .expect("precision >= 1");

    // Classify against the rim using i128 (|floor| < 2^64 here).
    // round_bigfloat_to_i64 saturates, so compare via the big value:
    // build the rim constants exactly and partial_cmp.
    let min_big = {
        let (v, s) = BigFloat::try_from_i64_exact(1, 66)
            .expect("precision >= 1")
            .scale_by_pow2(63);
        debug_assert!(s.is_ok());
        v.negated() // −2^63 = i64::MIN
    };
    let min_minus_1 = {
        let one = BigFloat::try_from_i64_exact(1, 66).expect("precision >= 1");
        min_big.sub(&one, RoundingMode::NearestEven).0 // exact at 66 bits
    };
    let max_plus_1 = min_big.negated(); // +2^63 = i64::MAX + 1

    // Integers with |x| ≥ 2^63 are NOT caught by the exact dispatch
    // upstream (integer_exponent rejects magnitudes past i64), so
    // detect exactness here: floor(x) == x. Two rows need it (found
    // by the ADR-0101 adversarial verification): x = −2^63 is a
    // representable exact power (the window path's on-grid value
    // would defeat the Ziv interval test and force a spurious
    // INEXACT), and x = −2^63 − 1 puts the truth EXACTLY at the
    // MinPos/2 tie, where the sliver's strict-interior justification
    // fails and the tie must be resolved explicitly.
    let x_is_integer = {
        let (gap, _) = x.sub(&n_big, RoundingMode::NearestEven);
        gap.is_zero()
    };

    use core::cmp::Ordering;
    let ge = |a: &BigFloat, b: &BigFloat| {
        matches!(
            a.partial_cmp(b).0,
            Some(Ordering::Greater | Ordering::Equal)
        )
    };
    if ge(&n_big, &max_plus_1) {
        return exp_certain_overflow(target_precision, mode);
    }
    if !ge(&n_big, &min_minus_1) {
        return exp_deep_underflow(target_precision, mode);
    }
    if !ge(&n_big, &min_big) {
        if x_is_integer {
            // x = i64::MIN − 1 exactly: the truth 2^x is EXACTLY
            // MinPos/2, the to-nearest tie between +0 and MinPos.
            // Resolve NE to +0 (the IEEE analogue: zero carries the
            // even significand; pfloat's no-subnormal grid leaves
            // both candidates' mantissa lsbs degenerate, so the
            // convention is recorded here). Away/upward take MinPos,
            // inward +0. The truth is below every representable:
            // UNDERFLOW|INEXACT.
            let value = match mode {
                RoundingMode::NearestAway | RoundingMode::TowardPositive => {
                    BigFloat::try_new_zero(Sign::Positive, target_precision)
                        .expect("precision >= 1")
                        .next_up()
                        .0
                }
                _ => BigFloat::try_new_zero(Sign::Positive, target_precision)
                    .expect("precision >= 1"),
            };
            let status = Status::UNDERFLOW | Status::INEXACT;
            auto_raise(status);
            return (value, status);
        }
        // Non-integer x: the truth is strictly inside
        // (MinPos/2, MinPos) and the shared sliver dispatch holds.
        return exp_underflow_sliver(target_precision, mode);
    }

    if x_is_integer {
        // x = −2^63 (the only window integer reachable: positives
        // and shallower negatives with ≤62-bit spans were dispatched
        // exactly upstream): 2^x is exactly representable at every
        // precision. The window path's on-grid s = 0.5 would defeat
        // the Ziv interval test (feedback_exact_value_defeats_ziv)
        // and exhaust into a spurious INEXACT; return the exact
        // power directly, Status::OK under every mode.
        let n = round_bigfloat_to_i64(&n_big);
        let (value, scale_status) = BigFloat::try_from_i64_exact(1, target_precision)
            .expect("precision >= 1")
            .scale_by_pow2(n);
        debug_assert!(scale_status.is_ok(), "n ∈ [i64::MIN, i64::MAX] window");
        return (value, Status::OK);
    }

    // Representable window. k = floor(x) fits i64; clamp one inside
    // the rim (frac shifts by the clamp so x = k + frac still holds,
    // keeping s = 2^frac ∈ [1, 4) with exponent(s) ∈ {0, 1} and the
    // compose exact for in-range truths; a certified carry past the
    // top binade routes to the overflow dispatch).
    let n = round_bigfloat_to_i64(&n_big);
    let k = n.clamp(i64::MIN + 1, i64::MAX - 1);
    let k_big = BigFloat::try_from_i64_exact(k, 66).expect("precision >= 1");
    let (s, ziv_status) = ziv_round(
        |w| {
            let x_w = x
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            // frac = x − k: exact (both ends on x's grid, |frac| < 2).
            let (frac, _) = x_w.sub(&k_big, RoundingMode::NearestEven);
            let ln_2 = ln_2_at(w);
            let (product, _) = frac.mul(&ln_2, RoundingMode::NearestEven);
            let (e_val, _) = product.exp(RoundingMode::NearestEven);
            e_val
        },
        target_precision,
        mode,
        EXP2_ERROR_GUARD,
    );
    let (result, scale_status) = s.scale_by_pow2(k);
    if scale_status.overflow() {
        return exp_certain_overflow(target_precision, mode);
    }
    let status = ziv_status | scale_status | Status::INEXACT;
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
    fn exp2_zero_is_one() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, _) = z.exp2(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    #[test]
    fn exp2_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.exp2(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn exp2_neg_inf_is_zero() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, _) = ni.exp2(RoundingMode::NearestEven);
        assert!(r.is_zero() && r.is_sign_positive());
    }

    #[test]
    fn exp2_one_is_two() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.exp2(RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        assert!(close_at(&r, &two, 113 - 8));
    }

    #[test]
    fn exp2_ten_is_1024() {
        let ten = BigFloat::try_from_i64_exact(10, 113).unwrap();
        let (r, _) = ten.exp2(RoundingMode::NearestEven);
        let k = BigFloat::try_from_i64_exact(1024, 113).unwrap();
        assert!(close_at(&r, &k, 113 - 12));
    }

    #[test]
    fn exp2_negative_ten() {
        // 2^-10 = 1/1024
        let neg_ten = BigFloat::try_from_i64_exact(-10, 113).unwrap();
        let (r, _) = neg_ten.exp2(RoundingMode::NearestEven);
        let recip = BigFloat::parse_str("0.0009765625", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        assert!(close_at(&r, &recip, 113 - 12));
    }

    #[test]
    fn exp2_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.exp2(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn exp2_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.exp2(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixed_exp2() {
        let one = FixedFloat::<53>::try_from_i64_exact(1).unwrap();
        let (r, _) = one.exp2(RoundingMode::NearestEven);
        let two = FixedFloat::<53>::try_from_i64_exact(2).unwrap();
        assert_eq!(r.partial_cmp(&two).0, Some(Ordering::Equal));
    }
}
