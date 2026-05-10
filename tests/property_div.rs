//! Property-based tests for [`pfloat::BigFloat::div`].

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

fn arb_nonzero_small_int() -> impl Strategy<Value = i64> {
    prop_oneof![-(1i64 << 20)..=-1i64, 1i64..=(1i64 << 20),]
}

fn arb_finite_nonzero() -> impl Strategy<Value = BigFloat> {
    (arb_nonzero_small_int(), arb_precision()).prop_filter_map("exact-fit", |(n, p)| {
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
    /// `a / a == 1` for finite non-zero `a`.
    #[test]
    fn self_div_self_is_one(a in arb_finite_nonzero()) {
        let (q, status) = a.div(&a, RoundingMode::NearestEven);
        prop_assert!(status.is_ok() || status.inexact());
        let one = BigFloat::try_from_i64_exact(1, a.precision()).unwrap();
        prop_assert_eq!(q.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    /// `a / 1 == a` for finite `a`.
    #[test]
    fn div_by_one_identity(a in arb_finite_nonzero(), mode in arb_mode()) {
        let one = BigFloat::try_from_i64_exact(1, a.precision()).unwrap();
        let (q, _) = a.div(&one, mode);
        prop_assert_eq!(q.partial_cmp(&a).0, Some(Ordering::Equal));
    }

    /// Sign rule: `sign(a / b) == sign(a) XOR sign(b)` for non-zero
    /// non-NaN finite operands.
    #[test]
    fn div_sign_rule(a in arb_finite_nonzero(), b in arb_finite_nonzero()) {
        let (q, _) = a.div(&b, RoundingMode::NearestEven);
        if !q.is_zero() && !q.is_nan() {
            prop_assert_eq!(q.sign(), a.sign().xor(b.sign()));
        }
    }

    /// Mul/div round-trip: `(a × b) / b == a` (exact when there's no
    /// rounding loss; at sufficient precision, exact for integer
    /// inputs).
    #[test]
    fn mul_div_round_trip(
        a in arb_finite_nonzero(),
        b in arb_finite_nonzero(),
    ) {
        // Round both to high precision so the round-trip is exact.
        let p = 512u32;
        let (a_hi, _) = a.round_to_precision(p, RoundingMode::NearestEven).unwrap();
        let (b_hi, _) = b.round_to_precision(p, RoundingMode::NearestEven).unwrap();
        let (prod, _) = a_hi.mul(&b_hi, RoundingMode::NearestEven);
        let (back, _) = prod.div(&b_hi, RoundingMode::NearestEven);
        prop_assert_eq!(back.partial_cmp(&a_hi).0, Some(Ordering::Equal));
    }

    /// `finite_nonzero / 0` raises `DIV_BY_ZERO` and returns ±Inf.
    #[test]
    fn div_by_zero_flag(a in arb_finite_nonzero(), zero_sign in arb_sign()) {
        let z = BigFloat::try_new_zero(zero_sign, a.precision()).unwrap();
        let (q, status) = a.div(&z, RoundingMode::NearestEven);
        prop_assert!(q.is_infinite());
        prop_assert!(status.div_by_zero());
        prop_assert_eq!(q.sign(), a.sign().xor(zero_sign));
    }

    /// `0 / 0` is qNaN + INVALID.
    #[test]
    fn zero_div_zero_invalid(s1 in arb_sign(), s2 in arb_sign(), p in arb_precision()) {
        let z1 = BigFloat::try_new_zero(s1, p).unwrap();
        let z2 = BigFloat::try_new_zero(s2, p).unwrap();
        let (q, status) = z1.div(&z2, RoundingMode::NearestEven);
        prop_assert!(q.is_quiet_nan());
        prop_assert!(status.invalid());
    }

    /// `Inf / Inf` is qNaN + INVALID.
    #[test]
    fn inf_div_inf_invalid(s1 in arb_sign(), s2 in arb_sign(), p in arb_precision()) {
        let i1 = BigFloat::try_new_infinity(s1, p).unwrap();
        let i2 = BigFloat::try_new_infinity(s2, p).unwrap();
        let (q, status) = i1.div(&i2, RoundingMode::NearestEven);
        prop_assert!(q.is_quiet_nan());
        prop_assert!(status.invalid());
    }

    /// `Inf / finite` is signed Inf with combined sign.
    #[test]
    fn inf_div_finite(s_inf in arb_sign(), b in arb_finite_nonzero()) {
        let inf = BigFloat::try_new_infinity(s_inf, b.precision()).unwrap();
        let (q, _) = inf.div(&b, RoundingMode::NearestEven);
        prop_assert!(q.is_infinite());
        prop_assert_eq!(q.sign(), s_inf.xor(b.sign()));
    }

    /// `finite / Inf` is signed zero with combined sign.
    #[test]
    fn finite_div_inf(a in arb_finite_nonzero(), s_inf in arb_sign()) {
        let inf = BigFloat::try_new_infinity(s_inf, a.precision()).unwrap();
        let (q, _) = a.div(&inf, RoundingMode::NearestEven);
        prop_assert!(q.is_zero());
        prop_assert_eq!(q.sign(), a.sign().xor(s_inf));
    }

    /// NaN propagation.
    #[test]
    fn nan_propagates(a in arb_any_bigfloat(), b in arb_any_bigfloat(), mode in arb_mode()) {
        if a.is_nan() || b.is_nan() {
            let (q, _) = a.div(&b, mode);
            prop_assert!(q.is_nan());
        }
    }

    /// sNaN raises INVALID.
    #[test]
    fn snan_raises_invalid(a in arb_any_bigfloat(), b in arb_any_bigfloat(), mode in arb_mode()) {
        if a.is_signaling_nan() || b.is_signaling_nan() {
            let (_q, status) = a.div(&b, mode);
            prop_assert!(status.invalid());
        }
    }

    /// Result precision equals max of input precisions.
    #[test]
    fn div_result_precision_is_max(
        a in arb_finite_nonzero(),
        b in arb_finite_nonzero(),
    ) {
        let (q, _) = a.div(&b, RoundingMode::NearestEven);
        prop_assert_eq!(q.precision(), a.precision().max(b.precision()));
    }
}
