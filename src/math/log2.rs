//! `log2(x)`: base-2 logarithm.
//!
//! Composition: `log2(x) = ln(x) / ln(2)` at working precision
//! `target + 64`. Special cases flow through composition: `ln(x)`
//! handles the entire IEEE 754-2019 §9.2 dispatch (NaN, `±0` →
//! `−∞ + DIV_BY_ZERO`, negative finite → `qNaN + INVALID`, `±∞`),
//! and the subsequent division by the positive `ln(2)` preserves
//! sign, infinity, and NaN.
//!
//! pfloat stores normal `BigFloat` as `m × 2^e` with the top bit of
//! the mantissa set. An exact-power-of-two fast path would return
//! `e` directly without rounding; slice 3d skips that optimization
//! and lets the composition produce the same value via `ln`. The
//! result on a power-of-two input rounds to the exact integer at
//! the 64-bit guard precision.

use crate::big::{BigFloat, BuildError, Parts};
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::ln_2_at;

impl BigFloat {
    /// `log2(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn log2(&self, mode: RoundingMode) -> (Self, Status) {
        self.log2_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `log2(self)` with explicit result precision.
    pub fn log2_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(log2_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `log2(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::log2`].
    #[must_use]
    pub fn log2(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().log2(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn log2_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    // Power-of-two exact dispatch (pf-kk16, ADR-0039). For
    // x = 2^k with integer k, log2(2^k) = k exactly. The
    // composition `ln(x) / ln(2)` at working_prec returns
    // k + epsilon (epsilon ~ 2^-w from the division of two
    // working-precision approximations to ln(2^k) and ln(2));
    // under directed modes the rounded result would land 1 ULP
    // away from the exact integer. The pre-composition dispatch
    // detects x via its mantissa-bit-pattern and returns k
    // directly; if k does not fit at target_precision, fall
    // through to the composition (which converges correctly via
    // its 64-bit guard for k representable at target_precision).
    if let Some(k) = power_of_two_exponent(x) {
        if let Ok(result) = BigFloat::try_from_i64_exact(k, target_precision) {
            return (result, Status::OK);
        }
    }

    let working_prec = target_precision.saturating_add(64);
    let (ln_x, ln_status) = x
        .ln_round(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1");
    let ln_2 = ln_2_at(working_prec);
    let (ratio, div_status) = ln_x.div(&ln_2, RoundingMode::NearestEven);
    let (rounded, round_status) = ratio
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    let mut status = ln_status | div_status | round_status;
    // A finite normal result for x outside the exact set is an
    // irrational log2 rounded onto the grid ⟹ INEXACT, even where the
    // composition lands exactly on an integer (for a dyadic x = m·2^e
    // with odd m > 1, log2 x = e + log2 m and log2 of an odd m > 1 is
    // irrational; pf-njs5 over-report). The exact x = 2^k case is
    // dispatched above. Non-finite or domain results (qNaN + INVALID
    // from x < 0 or x = NaN, ±∞ / DIV_BY_ZERO from x = ±0 or x = +∞)
    // flow through the composition and keep their status untouched.
    if matches!(rounded.parts(), Parts::Normal { .. }) {
        status |= Status::INEXACT;
    }
    auto_raise(status);
    (rounded, status)
}

/// If `x` is `+2^k` for some `k ∈ i64`, returns `Some(k)`;
/// otherwise returns `None`. A `BigFloat` stores its mantissa in
/// little-endian limbs with the top bit of the most-significant
/// limb always set (ADR-0001). A power of two stores as the
/// most-significant limb equal to `1u64 << 63` with every other
/// limb zero; the value `2^k` then has `BigFloat::exponent = k`
/// (the integer interpretation `2^(precision-1)` scaled by
/// `2^(exponent - precision + 1)` collapses to `2^exponent`).
/// The check is `O(limbs)`.
fn power_of_two_exponent(x: &BigFloat) -> Option<i64> {
    match x.parts() {
        Parts::Normal {
            sign: Sign::Positive,
            exponent,
            mantissa,
            ..
        } => {
            let top_limb_idx = mantissa.len().checked_sub(1)?;
            const TOP_BIT_ONLY: u64 = 1u64 << 63;
            if mantissa[top_limb_idx] != TOP_BIT_ONLY {
                return None;
            }
            for &limb in &mantissa[..top_limb_idx] {
                if limb != 0 {
                    return None;
                }
            }
            Some(exponent)
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sign::Sign;
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
    fn log2_one_is_zero() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.log2(RoundingMode::NearestEven);
        assert!(r.is_zero());
    }

    #[test]
    fn log2_two_is_one() {
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (r, _) = two.log2(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        assert!(close_at(&r, &one, 113 - 12));
    }

    #[test]
    fn log2_eight_is_three() {
        let eight = BigFloat::try_from_i64_exact(8, 113).unwrap();
        let (r, _) = eight.log2(RoundingMode::NearestEven);
        let three = BigFloat::try_from_i64_exact(3, 113).unwrap();
        assert!(close_at(&r, &three, 113 - 12));
    }

    #[test]
    fn log2_powers_of_two_are_exact_under_every_directed_mode() {
        // pf-kk16 pinning test: the power-of-two pre-composition
        // dispatch returns the exact integer k under every mode for
        // x = 2^k. Without the dispatch, the ln(x)/ln(2) composition
        // returns k + epsilon and TP rounds up to k + ULP for
        // sufficiently small target precision
        // (feedback_exact_value_defeats_ziv).
        for &mode in &[
            RoundingMode::NearestEven,
            RoundingMode::TowardPositive,
            RoundingMode::TowardNegative,
            RoundingMode::TowardZero,
            RoundingMode::NearestAway,
        ] {
            for &prec in &[24u32, 53, 113] {
                for &k in &[0i64, 1, 2, 3, 8, 10, 20] {
                    let two_to_k = if k == 0 {
                        BigFloat::try_from_i64_exact(1, prec).unwrap()
                    } else {
                        BigFloat::try_from_i64_exact(1i64 << k, prec).unwrap()
                    };
                    let (r, status) = two_to_k.log2(mode);
                    assert!(
                        status.is_ok(),
                        "log2(2^{k}) status under {mode:?}@p{prec}: {status:?}"
                    );
                    let expected = BigFloat::try_from_i64_exact(k, prec).unwrap();
                    assert_eq!(
                        r.partial_cmp(&expected).0,
                        Some(Ordering::Equal),
                        "log2(2^{k}) = {k} expected under {mode:?}@p{prec}, got {r:?}"
                    );
                    assert_eq!(r.precision(), prec);
                }
            }
        }
    }

    #[test]
    fn log2_one_is_zero_under_every_directed_mode() {
        // pf-kk16: log2(1) = log2(2^0) = 0 exactly under every mode.
        // Covered by the power-of-two dispatch (k=0). Also covered
        // by ln(1)=0 through the composition path; this test pins
        // the pre-composition dispatch.
        for &mode in &[
            RoundingMode::NearestEven,
            RoundingMode::TowardPositive,
            RoundingMode::TowardNegative,
            RoundingMode::TowardZero,
            RoundingMode::NearestAway,
        ] {
            for &prec in &[24u32, 53, 113] {
                let one = BigFloat::try_from_i64_exact(1, prec).unwrap();
                let (r, _) = one.log2(mode);
                assert!(
                    r.is_zero() && !r.is_sign_negative(),
                    "log2(1) should be +0 under {mode:?}@p{prec}, got {r:?}"
                );
            }
        }
    }

    #[test]
    fn log2_non_power_of_two_falls_through() {
        // pf-kk16 regression guard: non-power-of-two inputs still
        // route through the ln(x)/ln(2) composition (i.e., they are
        // NOT silently snapped to an integer by the dispatch).
        // log2(3) ≈ 1.585; check the result is close to 1.585, not
        // exactly 1 or 2.
        let three = BigFloat::try_from_i64_exact(3, 113).unwrap();
        let (r, _) = three.log2(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        assert_ne!(
            r.partial_cmp(&one).0,
            Some(Ordering::Equal),
            "log2(3) must not snap to 1"
        );
        assert_ne!(
            r.partial_cmp(&two).0,
            Some(Ordering::Equal),
            "log2(3) must not snap to 2"
        );
    }

    #[test]
    fn log2_zero_is_neg_inf() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, status) = z.log2(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative());
        assert!(status.div_by_zero());
    }

    #[test]
    fn log2_negative_is_invalid() {
        let neg = BigFloat::try_from_i64_exact(-5, 53).unwrap();
        let (r, status) = neg.log2(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn log2_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.log2(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn log2_neg_inf_is_invalid() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, status) = ni.log2(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn log2_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.log2(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn log2_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.log2(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn log2_round_trip_with_exp2() {
        for n in &[2i64, 3, 5, 7, 10, 100] {
            let x = BigFloat::try_from_i64_exact(*n, 113).unwrap();
            let (l, _) = x.log2(RoundingMode::NearestEven);
            let (back, _) = l.exp2(RoundingMode::NearestEven);
            assert!(close_at(&back, &x, 113 - 12), "log2/exp2 round-trip at {n}");
        }
    }
}
