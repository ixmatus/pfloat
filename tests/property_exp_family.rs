//! Property-based tests for the exp-family wrappers shipped in
//! slice 3d: `expm1`, `exp2`, `exp10`.

#![cfg(all(feature = "big", feature = "exp-log"))]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Sign};
use proptest::prelude::*;

fn arb_precision() -> impl Strategy<Value = u32> {
    prop_oneof![Just(53u32), Just(64u32), Just(113u32)]
}

fn close_within(a: &BigFloat, b: &BigFloat, bits: u32) -> bool {
    let (diff, _) = a.sub(b, RoundingMode::NearestEven);
    let abs_diff = diff.abs();
    if abs_diff.is_zero() {
        return true;
    }
    let p = a.precision().max(b.precision());
    let two = BigFloat::try_from_i64_exact(2, p).unwrap();
    let one = BigFloat::try_from_i64_exact(1, p).unwrap();
    let abs_b = b.abs();
    let mut bound = if abs_b.is_zero() { one } else { abs_b };
    for _ in 0..bits {
        bound = bound.div(&two, RoundingMode::NearestEven).0;
    }
    matches!(
        abs_diff.partial_cmp(&bound).0,
        Some(Ordering::Less | Ordering::Equal)
    )
}

proptest! {
    // ---- expm1 ----

    /// expm1(0) = 0 at any precision.
    #[test]
    fn expm1_zero(p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.expm1(RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
    }

    /// expm1(x) + 1 ≈ exp(x) for x ≥ 0 where the addition does not
    /// cancel. For x < 0 with large |x|, expm1(x) → −1 and the
    /// `+ 1` operation cancels ~|x · log2 e| bits; that regime is
    /// covered by the inline `expm1_neg_inf_is_neg_one` test, not
    /// this property.
    #[test]
    fn expm1_plus_one_matches_exp(n in 0i64..=20) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (r, _) = x.expm1(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (back, _) = r.add(&one, RoundingMode::NearestEven);
        let (e_x, _) = x.exp(RoundingMode::NearestEven);
        prop_assert!(
            close_within(&back, &e_x, p.saturating_sub(12)),
            "expm1({n}) + 1 = {back}, exp({n}) = {e_x}",
        );
    }

    /// expm1 is sign-preserving on small inputs (expm1(x) and x agree in sign).
    #[test]
    fn expm1_preserves_sign_small(n in -10i64..=10) {
        if n == 0 { return Ok(()); }
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        // Divide by 2^20 to get a small but nonzero value.
        let mut small = x;
        for _ in 0..20 {
            small = small.div(&two, RoundingMode::NearestEven).0;
        }
        let (r, _) = small.expm1(RoundingMode::NearestEven);
        prop_assert_eq!(r.is_sign_negative(), small.is_sign_negative());
    }

    // ---- exp2 ----

    /// exp2(0) = 1.
    #[test]
    fn exp2_zero_is_one(p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.exp2(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        prop_assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    /// exp2(k) = 2^k for small integer k.
    #[test]
    fn exp2_integer_matches_two_to_the_k(k in 0i64..=20) {
        let p = 113u32;
        let kbig = BigFloat::try_from_i64_exact(k, p).unwrap();
        let (r, _) = kbig.exp2(RoundingMode::NearestEven);
        let expected = BigFloat::try_from_i64_exact(1i64 << k, p).unwrap();
        prop_assert!(
            close_within(&r, &expected, p.saturating_sub(12)),
            "exp2({k}) = {r}, expected 2^{k}",
        );
    }

    /// exp2(a + b) = exp2(a) · exp2(b).
    #[test]
    fn exp2_additive(a in -10i64..=10, b in -10i64..=10) {
        let p = 113u32;
        let av = BigFloat::try_from_i64_exact(a, p).unwrap();
        let bv = BigFloat::try_from_i64_exact(b, p).unwrap();
        let (sum, _) = av.add(&bv, RoundingMode::NearestEven);
        let (lhs, _) = sum.exp2(RoundingMode::NearestEven);
        let (ea, _) = av.exp2(RoundingMode::NearestEven);
        let (eb, _) = bv.exp2(RoundingMode::NearestEven);
        let (rhs, _) = ea.mul(&eb, RoundingMode::NearestEven);
        prop_assert!(
            close_within(&lhs, &rhs, p.saturating_sub(12)),
            "exp2({a}+{b}) = {lhs}, exp2({a})·exp2({b}) = {rhs}",
        );
    }

    // ---- exp10 ----

    /// exp10(0) = 1.
    #[test]
    fn exp10_zero_is_one(p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.exp10(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        prop_assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    /// exp10(k) = 10^k for small integer k.
    #[test]
    fn exp10_integer_matches_ten_to_the_k(k in 0i64..=8) {
        let p = 113u32;
        let kbig = BigFloat::try_from_i64_exact(k, p).unwrap();
        let (r, _) = kbig.exp10(RoundingMode::NearestEven);
        let mut expected = BigFloat::try_from_i64_exact(1, p).unwrap();
        for _ in 0..k {
            let ten = BigFloat::try_from_i64_exact(10, p).unwrap();
            expected = expected.mul(&ten, RoundingMode::NearestEven).0;
        }
        prop_assert!(
            close_within(&r, &expected, p.saturating_sub(12)),
            "exp10({k}) = {r}, expected 10^{k}",
        );
    }

    /// exp10(a + b) = exp10(a) · exp10(b).
    #[test]
    fn exp10_additive(a in -5i64..=5, b in -5i64..=5) {
        let p = 113u32;
        let av = BigFloat::try_from_i64_exact(a, p).unwrap();
        let bv = BigFloat::try_from_i64_exact(b, p).unwrap();
        let (sum, _) = av.add(&bv, RoundingMode::NearestEven);
        let (lhs, _) = sum.exp10(RoundingMode::NearestEven);
        let (ea, _) = av.exp10(RoundingMode::NearestEven);
        let (eb, _) = bv.exp10(RoundingMode::NearestEven);
        let (rhs, _) = ea.mul(&eb, RoundingMode::NearestEven);
        prop_assert!(
            close_within(&lhs, &rhs, p.saturating_sub(16)),
            "exp10({a}+{b}) = {lhs}, exp10({a})·exp10({b}) = {rhs}",
        );
    }
}
