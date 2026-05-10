//! Property-based tests for [`pfloat::BigFloat::fma`].

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

fn arb_small_int() -> impl Strategy<Value = i64> {
    -(1i64 << 20)..(1i64 << 20)
}

fn arb_finite() -> impl Strategy<Value = BigFloat> {
    (arb_small_int(), arb_precision()).prop_filter_map("exact-fit", |(n, p)| {
        BigFloat::try_from_i64_exact(n, p).ok()
    })
}

fn arb_any() -> impl Strategy<Value = BigFloat> {
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
    /// `fma(a, b, 0) == a × b` (when c is +0).
    #[test]
    fn fma_with_zero_c_equals_mul(a in arb_finite(), b in arb_finite(), mode in arb_mode()) {
        let zero = BigFloat::try_new_zero(Sign::Positive, a.precision().max(b.precision())).unwrap();
        let (fma_r, _) = a.fma(&b, &zero, mode);
        let (mul_r, _) = a.mul(&b, mode);
        // Result precisions may differ if fma's target_precision is
        // larger; round both to the same comparison precision.
        let p = fma_r.precision().max(mul_r.precision());
        let (fma_p, _) = fma_r.round_to_precision(p, mode).unwrap();
        let (mul_p, _) = mul_r.round_to_precision(p, mode).unwrap();
        prop_assert_eq!(fma_p.partial_cmp(&mul_p).0, Some(Ordering::Equal));
    }

    /// `fma(0, b, c) == c` (and `fma(a, 0, c) == c` for finite c).
    #[test]
    fn fma_with_zero_a_equals_c(b in arb_finite(), c in arb_finite(), mode in arb_mode()) {
        let zero = BigFloat::try_new_zero(Sign::Positive, b.precision()).unwrap();
        let (fma_r, _) = zero.fma(&b, &c, mode);
        let target_p = zero.precision().max(b.precision()).max(c.precision());
        let (c_promoted, _) = c.round_to_precision(target_p, mode).unwrap();
        prop_assert_eq!(fma_r.partial_cmp(&c_promoted).0, Some(Ordering::Equal));
    }

    /// `fma(1, b, c) == b + c`.
    #[test]
    fn fma_with_one_a_equals_b_plus_c(b in arb_finite(), c in arb_finite(), mode in arb_mode()) {
        let one = BigFloat::try_from_i64_exact(1, b.precision().max(c.precision())).unwrap();
        let (fma_r, _) = one.fma(&b, &c, mode);
        let (add_r, _) = b.add(&c, mode);
        let p = fma_r.precision().max(add_r.precision());
        let (fma_p, _) = fma_r.round_to_precision(p, mode).unwrap();
        let (add_p, _) = add_r.round_to_precision(p, mode).unwrap();
        prop_assert_eq!(fma_p.partial_cmp(&add_p).0, Some(Ordering::Equal));
    }

    /// `fma(a, b, c) == fma(b, a, c)` (commutativity in the product).
    #[test]
    fn fma_commutative_in_product(
        a in arb_finite(),
        b in arb_finite(),
        c in arb_finite(),
        mode in arb_mode(),
    ) {
        let (ab_c, _) = a.fma(&b, &c, mode);
        let (ba_c, _) = b.fma(&a, &c, mode);
        prop_assert_eq!(ab_c.partial_cmp(&ba_c).0, Some(Ordering::Equal));
    }

    /// NaN in any operand propagates.
    #[test]
    fn nan_propagates(a in arb_any(), b in arb_any(), c in arb_any(), mode in arb_mode()) {
        if a.is_nan() || b.is_nan() || c.is_nan() {
            let (r, _) = a.fma(&b, &c, mode);
            prop_assert!(r.is_nan());
        }
    }

    /// sNaN raises INVALID.
    #[test]
    fn snan_raises_invalid(a in arb_any(), b in arb_any(), c in arb_any(), mode in arb_mode()) {
        if a.is_signaling_nan() || b.is_signaling_nan() || c.is_signaling_nan() {
            let (_r, status) = a.fma(&b, &c, mode);
            prop_assert!(status.invalid());
        }
    }

    /// `0 × ∞` with c-not-NaN raises INVALID.
    #[test]
    fn zero_inf_with_finite_c_invalid(
        sz in arb_sign(),
        si in arb_sign(),
        c in arb_finite(),
    ) {
        let z = BigFloat::try_new_zero(sz, 53).unwrap();
        let inf = BigFloat::try_new_infinity(si, 53).unwrap();
        let (r, status) = z.fma(&inf, &c, RoundingMode::NearestEven);
        prop_assert!(r.is_quiet_nan());
        prop_assert!(status.invalid());
    }

    /// `0 × ∞` with c-is-qNaN propagates the qNaN without raising
    /// INVALID (per IEEE 754-2019 §7.2 note).
    #[test]
    fn zero_inf_with_qnan_c_suppresses_invalid(
        sz in arb_sign(),
        si in arb_sign(),
    ) {
        let z = BigFloat::try_new_zero(sz, 53).unwrap();
        let inf = BigFloat::try_new_infinity(si, 53).unwrap();
        let qnan = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = z.fma(&inf, &qnan, RoundingMode::NearestEven);
        prop_assert!(r.is_quiet_nan());
        prop_assert!(!status.invalid());
    }

    /// Result precision equals max of three input precisions.
    #[test]
    fn fma_result_precision_is_max(a in arb_finite(), b in arb_finite(), c in arb_finite()) {
        let (r, _) = a.fma(&b, &c, RoundingMode::NearestEven);
        prop_assert_eq!(
            r.precision(),
            a.precision().max(b.precision()).max(c.precision())
        );
    }
}
