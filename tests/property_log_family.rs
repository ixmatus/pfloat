//! Property-based tests for the log-family wrappers shipped in
//! slice 3d: `log1p`, `log2`, `log10`.

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
    // ---- log1p ----

    /// log1p(0) = 0.
    #[test]
    fn log1p_zero(p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.log1p(RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
    }

    /// log1p(x) = ln(1 + x) for moderate x > -1.
    #[test]
    fn log1p_matches_ln_one_plus_x(n in 0i64..=100) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (lhs, _) = x.log1p(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (one_plus_x, _) = one.add(&x, RoundingMode::NearestEven);
        let (rhs, _) = one_plus_x.ln(RoundingMode::NearestEven);
        prop_assert!(
            close_within(&lhs, &rhs, p.saturating_sub(12)),
            "log1p({n}) = {lhs}, ln(1+{n}) = {rhs}",
        );
    }

    /// log1p / expm1 round-trip: expm1(log1p(x)) ≈ x for x > -1 of
    /// moderate magnitude.
    #[test]
    fn log1p_expm1_round_trip(n in 0i64..=20) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (l, _) = x.log1p(RoundingMode::NearestEven);
        let (back, _) = l.expm1(RoundingMode::NearestEven);
        prop_assert!(
            close_within(&back, &x, p.saturating_sub(16)),
            "expm1(log1p({n})) = {back}",
        );
    }

    // ---- log2 ----

    /// log2(1) = 0.
    #[test]
    fn log2_one_is_zero(p in arb_precision()) {
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (r, _) = one.log2(RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
    }

    /// log2(2^k) ≈ k for small k ≥ 0.
    #[test]
    fn log2_power_of_two(k in 0u32..=20) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(1i64 << k, p).unwrap();
        let (r, _) = x.log2(RoundingMode::NearestEven);
        let expected = BigFloat::try_from_i64_exact(i64::from(k), p).unwrap();
        prop_assert!(
            close_within(&r, &expected, p.saturating_sub(12)),
            "log2(2^{k}) = {r}",
        );
    }

    /// log2 monotonic on positive integers.
    #[test]
    fn log2_monotonic(n in 1i64..=1000) {
        let p = 113u32;
        let a = BigFloat::try_from_i64_exact(n, p).unwrap();
        let b = BigFloat::try_from_i64_exact(n + 1, p).unwrap();
        let (la, _) = a.log2(RoundingMode::NearestEven);
        let (lb, _) = b.log2(RoundingMode::NearestEven);
        prop_assert_eq!(lb.partial_cmp(&la).0, Some(Ordering::Greater));
    }

    /// log2(a · b) = log2(a) + log2(b).
    #[test]
    fn log2_product_is_sum(a in 1i64..=100, b in 1i64..=100) {
        let p = 113u32;
        let av = BigFloat::try_from_i64_exact(a, p).unwrap();
        let bv = BigFloat::try_from_i64_exact(b, p).unwrap();
        let (prod, _) = av.mul(&bv, RoundingMode::NearestEven);
        let (lhs, _) = prod.log2(RoundingMode::NearestEven);
        let (la, _) = av.log2(RoundingMode::NearestEven);
        let (lb, _) = bv.log2(RoundingMode::NearestEven);
        let (rhs, _) = la.add(&lb, RoundingMode::NearestEven);
        prop_assert!(
            close_within(&lhs, &rhs, p.saturating_sub(12)),
            "log2({a}·{b}) = {lhs}, sum = {rhs}",
        );
    }

    // ---- log10 ----

    /// log10(1) = 0.
    #[test]
    fn log10_one_is_zero(p in arb_precision()) {
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (r, _) = one.log10(RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
    }

    /// log10(10^k) ≈ k.
    #[test]
    fn log10_power_of_ten(k in 0i64..=8) {
        let p = 113u32;
        let mut x = BigFloat::try_from_i64_exact(1, p).unwrap();
        for _ in 0..k {
            let ten = BigFloat::try_from_i64_exact(10, p).unwrap();
            x = x.mul(&ten, RoundingMode::NearestEven).0;
        }
        let (r, _) = x.log10(RoundingMode::NearestEven);
        let expected = BigFloat::try_from_i64_exact(k, p).unwrap();
        prop_assert!(
            close_within(&r, &expected, p.saturating_sub(16)),
            "log10(10^{k}) = {r}",
        );
    }

    /// log10(a · b) = log10(a) + log10(b).
    #[test]
    fn log10_product_is_sum(a in 1i64..=100, b in 1i64..=100) {
        let p = 113u32;
        let av = BigFloat::try_from_i64_exact(a, p).unwrap();
        let bv = BigFloat::try_from_i64_exact(b, p).unwrap();
        let (prod, _) = av.mul(&bv, RoundingMode::NearestEven);
        let (lhs, _) = prod.log10(RoundingMode::NearestEven);
        let (la, _) = av.log10(RoundingMode::NearestEven);
        let (lb, _) = bv.log10(RoundingMode::NearestEven);
        let (rhs, _) = la.add(&lb, RoundingMode::NearestEven);
        prop_assert!(
            close_within(&lhs, &rhs, p.saturating_sub(16)),
            "log10({a}·{b}) = {lhs}, sum = {rhs}",
        );
    }
}
