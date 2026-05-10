//! Property-based tests for [`pfloat::BigFloat::mul`].

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

fn arb_small_finite() -> impl Strategy<Value = BigFloat> {
    // Values whose product comfortably fits in 53-bit precision.
    let small_int = -(1i64 << 20)..(1i64 << 20);
    (small_int, arb_precision()).prop_filter_map("exact-fit", |(n, p)| {
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
    /// `1 × a == a` for finite `a`.
    #[test]
    fn mul_identity(a in arb_small_finite(), mode in arb_mode()) {
        let one = BigFloat::try_from_i64_exact(1, a.precision()).unwrap();
        let (p, _) = a.mul(&one, mode);
        prop_assert_eq!(p.partial_cmp(&a).0, Some(Ordering::Equal));
    }

    /// `0 × a == 0` (with combined sign) for any finite `a`.
    #[test]
    fn mul_zero(a in arb_small_finite()) {
        let zero = BigFloat::try_new_zero(Sign::Positive, a.precision()).unwrap();
        let (p, status) = a.mul(&zero, RoundingMode::NearestEven);
        prop_assert!(status.is_ok());
        prop_assert!(p.is_zero());
    }

    /// Commutativity: `a × b == b × a` (exact for any precision wide
    /// enough; non-exact paths still produce equal results because
    /// the kernel is symmetric in its arguments).
    #[test]
    fn mul_commutative(
        a in arb_small_finite(),
        b in arb_small_finite(),
        mode in arb_mode(),
    ) {
        let (ab, _) = a.mul(&b, mode);
        let (ba, _) = b.mul(&a, mode);
        prop_assert_eq!(ab.partial_cmp(&ba).0, Some(Ordering::Equal));
    }

    /// Sign rule: `sign(a × b) == sign(a) XOR sign(b)` for non-zero
    /// non-NaN operands.
    #[test]
    fn mul_sign_rule(a in arb_small_finite(), b in arb_small_finite()) {
        if a.is_zero() || b.is_zero() {
            return Ok(()); // sign of zero handled separately
        }
        let (p, _) = a.mul(&b, RoundingMode::NearestEven);
        if !p.is_zero() {
            let expected = a.sign().xor(b.sign());
            prop_assert_eq!(p.sign(), expected);
        }
    }

    /// `a × (-1) == -a` for any finite `a`.
    #[test]
    fn mul_by_neg_one_negates(a in arb_small_finite()) {
        let neg_one = BigFloat::try_from_i64_exact(-1, a.precision()).unwrap();
        let (p, _) = a.mul(&neg_one, RoundingMode::NearestEven);
        let neg_a = a.negated();
        prop_assert_eq!(p.partial_cmp(&neg_a).0, Some(Ordering::Equal));
    }

    /// Result precision equals max of input precisions.
    #[test]
    fn mul_result_precision_is_max(a in arb_small_finite(), b in arb_small_finite()) {
        let (p, _) = a.mul(&b, RoundingMode::NearestEven);
        prop_assert_eq!(p.precision(), a.precision().max(b.precision()));
    }

    /// NaN propagation: any mul with NaN operand returns NaN.
    #[test]
    fn nan_propagates(a in arb_any_bigfloat(), b in arb_any_bigfloat(), mode in arb_mode()) {
        if a.is_nan() || b.is_nan() {
            let (p, _) = a.mul(&b, mode);
            prop_assert!(p.is_nan());
        }
    }

    /// sNaN raises INVALID.
    #[test]
    fn snan_raises_invalid(a in arb_any_bigfloat(), b in arb_any_bigfloat(), mode in arb_mode()) {
        if a.is_signaling_nan() || b.is_signaling_nan() {
            let (_p, status) = a.mul(&b, mode);
            prop_assert!(status.invalid());
        }
    }

    /// `Inf × 0` (either order) raises INVALID and returns qNaN.
    #[test]
    fn inf_times_zero_invalid(sign_inf in arb_sign(), sign_zero in arb_sign(), p in arb_precision()) {
        let inf = BigFloat::try_new_infinity(sign_inf, p).unwrap();
        let zero = BigFloat::try_new_zero(sign_zero, p).unwrap();
        let (r, status) = inf.mul(&zero, RoundingMode::NearestEven);
        prop_assert!(r.is_quiet_nan());
        prop_assert!(status.invalid());
        let (r2, status2) = zero.mul(&inf, RoundingMode::NearestEven);
        prop_assert!(r2.is_quiet_nan());
        prop_assert!(status2.invalid());
    }

    /// `Inf × finite_nonzero` is `±Inf` with combined sign.
    #[test]
    fn inf_times_finite(
        sign_inf in arb_sign(),
        a in arb_small_finite(),
    ) {
        // Skip zero operands.
        prop_assume!(!a.is_zero());
        let inf = BigFloat::try_new_infinity(sign_inf, a.precision()).unwrap();
        let (r, _) = inf.mul(&a, RoundingMode::NearestEven);
        prop_assert!(r.is_infinite());
        prop_assert_eq!(r.sign(), sign_inf.xor(a.sign()));
    }

    /// Squaring a value is positive (sign-rule property).
    #[test]
    fn square_is_positive(a in arb_small_finite()) {
        if a.is_zero() {
            return Ok(());
        }
        let (p, _) = a.mul(&a, RoundingMode::NearestEven);
        prop_assert!(p.is_sign_positive());
    }
}
