//! Property-based tests for the gamma family shipped in slice 4b:
//! `gamma`, `lgamma`.

#![cfg(all(feature = "big", feature = "specials"))]

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
    /// gamma(1) = 1.
    #[test]
    fn gamma_one(p in arb_precision()) {
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (r, _) = one.gamma(RoundingMode::NearestEven);
        prop_assert!(close_within(&r, &one, p.saturating_sub(16)));
    }

    /// gamma(n+1) = n · gamma(n). Recurrence identity.
    #[test]
    fn gamma_recurrence(n in 1i64..=10) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let x_plus_1 = BigFloat::try_from_i64_exact(n + 1, p).unwrap();
        let (gx, _) = x.gamma(RoundingMode::NearestEven);
        let (gx_plus_1, _) = x_plus_1.gamma(RoundingMode::NearestEven);
        let (lhs, _) = x.mul(&gx, RoundingMode::NearestEven);
        prop_assert!(
            close_within(&lhs, &gx_plus_1, p.saturating_sub(16)),
            "gamma({})·{} = {}, gamma({}+1) = {}",
            n, n, lhs, n, gx_plus_1,
        );
    }

    /// gamma(n) = (n-1)! for positive integer n. Spot-check the
    /// factorial identity by reconstructing the factorial via
    /// repeated multiplication.
    #[test]
    fn gamma_matches_factorial(n in 1i64..=12) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (r, _) = x.gamma(RoundingMode::NearestEven);
        let mut expected = BigFloat::try_from_i64_exact(1, p).unwrap();
        for k in 1..n {
            let k_big = BigFloat::try_from_i64_exact(k, p).unwrap();
            expected = expected.mul(&k_big, RoundingMode::NearestEven).0;
        }
        prop_assert!(
            close_within(&r, &expected, p.saturating_sub(16)),
            "gamma({n}) = {r}, expected (n-1)! = {expected}",
        );
    }

    /// lgamma(n+1) = lgamma(n) + ln(n).
    #[test]
    fn lgamma_recurrence(n in 1i64..=20) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let x_plus_1 = BigFloat::try_from_i64_exact(n + 1, p).unwrap();
        let (lx, _) = x.lgamma(RoundingMode::NearestEven);
        let (lx_plus_1, _) = x_plus_1.lgamma(RoundingMode::NearestEven);
        let (ln_n, _) = x.ln(RoundingMode::NearestEven);
        let (lhs, _) = lx.add(&ln_n, RoundingMode::NearestEven);
        prop_assert!(
            close_within(&lhs, &lx_plus_1, p.saturating_sub(20)),
            "lgamma({}) + ln({}) = {}, lgamma({}+1) = {}",
            n, n, lhs, n, lx_plus_1,
        );
    }

    /// gamma(x) > 0 for x > 0.
    #[test]
    fn gamma_positive_for_positive(num in 1i64..=20) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(num, p).unwrap();
        let (r, _) = x.gamma(RoundingMode::NearestEven);
        prop_assert!(r.is_sign_positive());
    }

    /// gamma sign alternates on (−n−1, −n) intervals: positive for
    /// x ∈ (−2, −1), negative for x ∈ (−1, 0), etc.
    #[test]
    fn gamma_sign_on_negative_intervals(k in 1u32..=4) {
        let p = 113u32;
        // Pick x = -(k - 1) - 0.5 = -k + 0.5. For k=1: x=-0.5 in (-1,0), Γ<0.
        // For k=2: x=-1.5 in (-2,-1), Γ>0. Etc.
        let half = BigFloat::parse_str("0.5", p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let k_big = BigFloat::try_from_i64_exact(i64::from(k), p).unwrap();
        let (x_pos, _) = k_big.sub(&half, RoundingMode::NearestEven);
        let x = x_pos.negated();
        let (r, _) = x.gamma(RoundingMode::NearestEven);
        let expect_negative = k % 2 == 1;
        prop_assert_eq!(r.is_sign_negative(), expect_negative);
    }
}
