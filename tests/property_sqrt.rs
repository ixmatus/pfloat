//! Property-based tests for [`pfloat::BigFloat::sqrt`].

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

fn arb_positive_int() -> impl Strategy<Value = i64> {
    1i64..(1i64 << 30)
}

fn arb_finite_positive() -> impl Strategy<Value = BigFloat> {
    (arb_positive_int(), arb_precision()).prop_filter_map("exact-fit", |(n, p)| {
        BigFloat::try_from_i64_exact(n, p).ok()
    })
}

proptest! {
    /// sqrt of any finite positive is positive (or +0).
    #[test]
    fn sqrt_is_nonneg(a in arb_finite_positive(), mode in arb_mode()) {
        let (r, _) = a.sqrt(mode);
        prop_assert!(!r.is_sign_negative() || r.is_zero(),
            "sqrt result should be non-negative, got {r:?}");
    }

    /// sqrt(x²) == x for finite positive x (when no rounding loss).
    #[test]
    fn sqrt_of_square_returns_original(x in arb_finite_positive()) {
        // Use a precision wide enough for the square to be exact.
        let target = 512u32;
        let (x_hi, _) = x.round_to_precision(target, RoundingMode::NearestEven).unwrap();
        let (square, _) = x_hi.mul(&x_hi, RoundingMode::NearestEven);
        let (sqrt_back, _) = square.sqrt(RoundingMode::NearestEven);
        prop_assert_eq!(sqrt_back.partial_cmp(&x_hi).0, Some(Ordering::Equal));
    }

    /// sqrt is monotonic: a <= b → sqrt(a) <= sqrt(b) for positive a, b.
    #[test]
    fn sqrt_is_monotonic(a in arb_finite_positive(), b in arb_finite_positive()) {
        // Lift to a common high precision so the sqrt is reliable.
        let p = 256u32;
        let (a_hi, _) = a.round_to_precision(p, RoundingMode::NearestEven).unwrap();
        let (b_hi, _) = b.round_to_precision(p, RoundingMode::NearestEven).unwrap();
        let (sa, _) = a_hi.sqrt(RoundingMode::NearestEven);
        let (sb, _) = b_hi.sqrt(RoundingMode::NearestEven);
        let original_ord = a_hi.partial_cmp(&b_hi).0;
        let sqrt_ord = sa.partial_cmp(&sb).0;
        match (original_ord, sqrt_ord) {
            (Some(Ordering::Less), Some(ord)) =>
                prop_assert!(ord != Ordering::Greater, "monotonicity violated"),
            (Some(Ordering::Equal), Some(ord)) =>
                prop_assert_eq!(ord, Ordering::Equal),
            (Some(Ordering::Greater), Some(ord)) =>
                prop_assert!(ord != Ordering::Less, "monotonicity violated"),
            _ => {}
        }
    }

    /// sqrt of a negative finite is qNaN + INVALID.
    #[test]
    fn sqrt_negative_invalid(x in arb_finite_positive()) {
        let neg = x.negated();
        let (r, status) = neg.sqrt(RoundingMode::NearestEven);
        prop_assert!(r.is_quiet_nan());
        prop_assert!(status.invalid());
    }

    /// sqrt(±0) = ±0.
    #[test]
    fn sqrt_zero(sign in arb_sign(), p in arb_precision()) {
        let z = BigFloat::try_new_zero(sign, p).unwrap();
        let (r, status) = z.sqrt(RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
        prop_assert_eq!(r.sign(), sign);
        prop_assert!(status.is_ok());
    }
}
