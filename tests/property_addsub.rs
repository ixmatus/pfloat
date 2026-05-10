//! Property-based tests for [`pfloat::BigFloat::add`] and
//! [`pfloat::BigFloat::sub`].
//!
//! Verifies IEEE 754-2019 §6.5 add/sub identities and edge cases
//! over arbitrary inputs.

#![cfg(feature = "big")]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Sign};
use proptest::prelude::*;

fn arb_sign() -> impl Strategy<Value = Sign> {
    prop_oneof![Just(Sign::Positive), Just(Sign::Negative)]
}

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
    prop_oneof![Just(53u32), Just(64u32), Just(113u32), Just(256u32)]
}

/// Generate any finite `BigFloat` constructable from `i64` at any
/// of a few useful precisions.
fn arb_finite_bigfloat() -> impl Strategy<Value = BigFloat> {
    (any::<i64>(), arb_precision()).prop_filter_map("exact-fit", |(n, p)| {
        BigFloat::try_from_i64_exact(n, p).ok()
    })
}

fn arb_any_bigfloat() -> impl Strategy<Value = BigFloat> {
    (arb_sign(), arb_precision(), 0..5u8, any::<i64>()).prop_filter_map(
        "construct",
        |(sign, prec, kind, n)| match kind {
            0 => BigFloat::try_new_zero(sign, prec).ok(),
            1 => BigFloat::try_new_infinity(sign, prec).ok(),
            2 => BigFloat::try_new_quiet_nan(sign, prec, &[]).ok(),
            3 => BigFloat::try_new_signaling_nan(sign, prec, &[]).ok(),
            _ => BigFloat::try_from_i64_exact(n, prec).ok(),
        },
    )
}

proptest! {
    /// `a + 0 == a` for any finite `a`.
    #[test]
    fn add_zero_identity(a in arb_finite_bigfloat(), mode in arb_mode()) {
        let zero = BigFloat::try_new_zero(Sign::Positive, a.precision()).unwrap();
        let (sum, _) = a.add(&zero, mode);
        prop_assert_eq!(sum.partial_cmp(&a).0, Some(Ordering::Equal));
    }

    /// `0 + a == a` for any finite `a` (commutativity check at zero).
    #[test]
    fn zero_plus_a(a in arb_finite_bigfloat(), mode in arb_mode()) {
        let zero = BigFloat::try_new_zero(Sign::Positive, a.precision()).unwrap();
        let (sum, _) = zero.add(&a, mode);
        prop_assert_eq!(sum.partial_cmp(&a).0, Some(Ordering::Equal));
    }

    /// `a - a == 0` for any finite `a`.
    #[test]
    fn self_minus_self_is_zero(a in arb_finite_bigfloat(), mode in arb_mode()) {
        let (diff, _) = a.sub(&a, mode);
        prop_assert!(diff.is_zero());
    }

    /// `a + (-a) == 0` for any finite `a` under round-to-nearest.
    #[test]
    fn add_negation_is_zero(a in arb_finite_bigfloat()) {
        let neg_a = a.negated();
        let (sum, _) = a.add(&neg_a, RoundingMode::NearestEven);
        prop_assert!(sum.is_zero());
    }

    /// Commutativity: `a + b == b + a` under any rounding mode.
    #[test]
    fn add_commutative(
        a in arb_finite_bigfloat(),
        b in arb_finite_bigfloat(),
        mode in arb_mode(),
    ) {
        let (ab, _) = a.add(&b, mode);
        let (ba, _) = b.add(&a, mode);
        prop_assert_eq!(ab.partial_cmp(&ba).0, Some(Ordering::Equal));
    }

    /// `a - b == -(b - a)` for any finite operands.
    #[test]
    fn sub_anti_symmetric(
        a in arb_finite_bigfloat(),
        b in arb_finite_bigfloat(),
        mode in arb_mode(),
    ) {
        let (ab, _) = a.sub(&b, mode);
        let (ba, _) = b.sub(&a, mode);
        let neg_ba = ba.negated();
        prop_assert_eq!(ab.partial_cmp(&neg_ba).0, Some(Ordering::Equal));
    }

    /// `a + b - b == a` when the operation is exact (precision wide
    /// enough). At sufficiently high precision (much wider than
    /// either operand's significant bits), the round-trip is exact.
    #[test]
    fn add_then_sub_round_trip_at_high_precision(
        a in arb_finite_bigfloat(),
        b in arb_finite_bigfloat(),
    ) {
        // Round both to a single high precision so the add/sub
        // result has room.
        let p = 512u32;
        let (a_hi, _) = a.round_to_precision(p, RoundingMode::NearestEven).unwrap();
        let (b_hi, _) = b.round_to_precision(p, RoundingMode::NearestEven).unwrap();
        let (sum, _) = a_hi.add(&b_hi, RoundingMode::NearestEven);
        let (back, _) = sum.sub(&b_hi, RoundingMode::NearestEven);
        prop_assert_eq!(back.partial_cmp(&a_hi).0, Some(Ordering::Equal));
    }

    /// NaN inputs propagate: any add/sub with a NaN operand returns
    /// a NaN.
    #[test]
    fn nan_propagates(a in arb_any_bigfloat(), b in arb_any_bigfloat(), mode in arb_mode()) {
        let nan_in = a.is_nan() || b.is_nan();
        if nan_in {
            let (sum, _) = a.add(&b, mode);
            prop_assert!(sum.is_nan(), "NaN in → NaN out");
        }
    }

    /// Signaling-NaN input raises INVALID.
    #[test]
    fn snan_raises_invalid(a in arb_any_bigfloat(), b in arb_any_bigfloat(), mode in arb_mode()) {
        if a.is_signaling_nan() || b.is_signaling_nan() {
            let (_sum, status) = a.add(&b, mode);
            prop_assert!(status.invalid());
        }
    }

    /// Inf + Inf same sign = Inf (same sign).
    #[test]
    fn inf_plus_inf_same_sign(s in arb_sign(), p in arb_precision()) {
        let inf = BigFloat::try_new_infinity(s, p).unwrap();
        let (sum, status) = inf.add(&inf, RoundingMode::NearestEven);
        prop_assert!(sum.is_infinite());
        prop_assert_eq!(sum.sign(), s);
        prop_assert!(status.is_ok());
    }

    /// Inf - Inf (same-sign sub) = NaN + INVALID.
    #[test]
    fn inf_minus_inf_is_invalid(s in arb_sign(), p in arb_precision()) {
        let inf = BigFloat::try_new_infinity(s, p).unwrap();
        let (diff, status) = inf.sub(&inf, RoundingMode::NearestEven);
        prop_assert!(diff.is_quiet_nan());
        prop_assert!(status.invalid());
    }

    /// Result precision equals max of input precisions.
    #[test]
    fn result_precision_is_max(
        a in arb_finite_bigfloat(),
        b in arb_finite_bigfloat(),
    ) {
        let (sum, _) = a.add(&b, RoundingMode::NearestEven);
        prop_assert_eq!(sum.precision(), a.precision().max(b.precision()));
    }
}
