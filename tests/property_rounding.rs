//! Property-based tests for the slice 1b rounding pipeline.
//!
//! Verifies the round-trip and direction invariants of
//! [`pfloat::BigFloat::round_to_precision`] and
//! [`pfloat::BigFloat::try_from_i64_round`] over arbitrary inputs.

#![cfg(feature = "big")]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Sign};
use proptest::prelude::*;

fn arb_mode() -> impl Strategy<Value = RoundingMode> {
    prop_oneof![
        Just(RoundingMode::NearestEven),
        Just(RoundingMode::NearestAway),
        Just(RoundingMode::TowardZero),
        Just(RoundingMode::TowardPositive),
        Just(RoundingMode::TowardNegative),
    ]
}

fn arb_precision() -> impl Strategy<Value = u32> {
    1u32..=256
}

proptest! {
    /// `try_from_i64_round` at exact-fit precision is exact under
    /// every mode (status is OK).
    #[test]
    fn from_i64_round_at_64_precision_exact(n in any::<i64>(), mode in arb_mode()) {
        let (_v, status) = BigFloat::try_from_i64_round(n, 64, mode).unwrap();
        prop_assert!(status.is_ok());
    }

    /// `try_from_i64_round` matches `try_from_i64_exact` whenever
    /// the latter succeeds.
    #[test]
    fn round_matches_exact_when_exact_fits(
        n in any::<i64>(),
        precision in arb_precision(),
        mode in arb_mode(),
    ) {
        if let Ok(exact) = BigFloat::try_from_i64_exact(n, precision) {
            let (rounded, status) = BigFloat::try_from_i64_round(n, precision, mode).unwrap();
            prop_assert!(status.is_ok());
            prop_assert_eq!(rounded, exact);
        }
    }

    /// `round_to_precision` extension to a wider precision is
    /// numerically equal and exact (status OK).
    #[test]
    fn extension_is_value_equal(n in any::<i64>(), precision in 1u32..=64, extra in 1u32..=128) {
        if let Ok(v) = BigFloat::try_from_i64_exact(n, precision) {
            let target = precision.saturating_add(extra);
            let (extended, status) = v.round_to_precision(target, RoundingMode::NearestEven).unwrap();
            prop_assert!(status.is_ok());
            prop_assert_eq!(extended.partial_cmp(&v).0, Some(Ordering::Equal));
            prop_assert_eq!(extended.precision(), target);
        }
    }

    /// `round_to_precision` from precision P to a wider precision
    /// then back to P (under any mode) is identity for finite
    /// values.
    #[test]
    fn extend_then_narrow_round_trip(
        n in any::<i64>(),
        precision in 2u32..=128,
        mode in arb_mode(),
    ) {
        if let Ok(v) = BigFloat::try_from_i64_exact(n, precision) {
            let wider = precision + 64;
            let (extended, _) = v.round_to_precision(wider, mode).unwrap();
            let (back, status) = extended.round_to_precision(precision, mode).unwrap();
            prop_assert!(status.is_ok(), "round-trip should be exact");
            prop_assert_eq!(back, v);
        }
    }

    /// `round_to_precision(self, self.precision, mode)` is identity
    /// for finite values (status OK).
    #[test]
    fn rounding_to_same_precision_is_identity(
        n in any::<i64>(),
        precision in arb_precision(),
        mode in arb_mode(),
    ) {
        if let Ok(v) = BigFloat::try_from_i64_exact(n, precision) {
            let (r, status) = v.round_to_precision(precision, mode).unwrap();
            prop_assert!(status.is_ok());
            prop_assert_eq!(r, v);
        }
    }

    /// Toward-zero rounding never produces a value larger in
    /// magnitude than the unrounded value.
    #[test]
    fn toward_zero_does_not_inflate_magnitude(
        n in any::<i64>(),
        precision in 1u32..=64,
    ) {
        // Build at a wider precision, round to narrower under TowardZero.
        let wide = precision + 8;
        if let Ok(wide_v) = BigFloat::try_from_i64_exact(n, wide) {
            if let Ok((narrow_v, _)) = wide_v.round_to_precision(precision, RoundingMode::TowardZero) {
                // |narrow_v| <= |wide_v|.
                let abs_narrow = narrow_v.abs();
                let abs_wide = wide_v.abs();
                let cmp = abs_narrow.partial_cmp(&abs_wide).0;
                prop_assert!(matches!(cmp, Some(Ordering::Less | Ordering::Equal)),
                    "TowardZero should not increase magnitude (got {:?})", cmp);
            }
        }
    }

    /// Toward-positive rounding never produces a value smaller
    /// than the unrounded value.
    #[test]
    fn toward_positive_does_not_decrease(
        n in any::<i64>(),
        precision in 1u32..=64,
    ) {
        let wide = precision + 8;
        if let Ok(wide_v) = BigFloat::try_from_i64_exact(n, wide) {
            if let Ok((narrow_v, _)) = wide_v.round_to_precision(precision, RoundingMode::TowardPositive) {
                let cmp = narrow_v.partial_cmp(&wide_v).0;
                prop_assert!(matches!(cmp, Some(Ordering::Greater | Ordering::Equal)),
                    "TowardPositive should round up (got {:?})", cmp);
            }
        }
    }

    /// Toward-negative rounding never produces a value larger
    /// than the unrounded value.
    #[test]
    fn toward_negative_does_not_increase(
        n in any::<i64>(),
        precision in 1u32..=64,
    ) {
        let wide = precision + 8;
        if let Ok(wide_v) = BigFloat::try_from_i64_exact(n, wide) {
            if let Ok((narrow_v, _)) = wide_v.round_to_precision(precision, RoundingMode::TowardNegative) {
                let cmp = narrow_v.partial_cmp(&wide_v).0;
                prop_assert!(matches!(cmp, Some(Ordering::Less | Ordering::Equal)),
                    "TowardNegative should round down (got {:?})", cmp);
            }
        }
    }

    /// `round_to_precision` preserves the kind for special values
    /// (Zero stays Zero, Infinity stays Infinity, NaN stays NaN).
    #[test]
    fn special_values_keep_kind(
        new_prec in arb_precision(),
        sign in prop_oneof![Just(Sign::Positive), Just(Sign::Negative)],
        mode in arb_mode(),
    ) {
        let pz = BigFloat::try_new_zero(sign, 53).unwrap();
        let pi = BigFloat::try_new_infinity(sign, 53).unwrap();
        let qn = BigFloat::try_new_quiet_nan(sign, 53, &[]).unwrap();
        let sn = BigFloat::try_new_signaling_nan(sign, 53, &[]).unwrap();

        let (r_z, _) = pz.round_to_precision(new_prec, mode).unwrap();
        prop_assert!(r_z.is_zero());
        prop_assert_eq!(r_z.sign(), sign);

        let (r_i, _) = pi.round_to_precision(new_prec, mode).unwrap();
        prop_assert!(r_i.is_infinite());
        prop_assert_eq!(r_i.sign(), sign);

        let (r_q, _) = qn.round_to_precision(new_prec, mode).unwrap();
        prop_assert!(r_q.is_quiet_nan());
        prop_assert_eq!(r_q.sign(), sign);

        let (r_s, _) = sn.round_to_precision(new_prec, mode).unwrap();
        prop_assert!(r_s.is_signaling_nan());
        prop_assert_eq!(r_s.sign(), sign);
    }
}
