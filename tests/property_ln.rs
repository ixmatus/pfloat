//! Property-based tests for [`pfloat::BigFloat::ln`].

#![cfg(all(feature = "big", feature = "exp-log"))]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Sign};
use proptest::prelude::*;

fn arb_precision() -> impl Strategy<Value = u32> {
    prop_oneof![Just(53u32), Just(64u32), Just(113u32)]
}

fn arb_positive_input() -> impl Strategy<Value = BigFloat> {
    (1i64..=10_000, arb_precision()).prop_filter_map("exact-fit", |(n, p)| {
        BigFloat::try_from_i64_exact(n, p).ok()
    })
}

/// `|a - b| <= 2^(expected_exp - bits)` where `expected_exp ≈
/// log2(|expected|)`. Loose tolerance because parsed-from-decimal
/// expected values lose accuracy past the round-trip digit count.
fn close_within(a: &BigFloat, b: &BigFloat, bits: u32) -> bool {
    let (diff, _) = a.sub(b, RoundingMode::NearestEven);
    let abs_diff = diff.abs();
    if abs_diff.is_zero() {
        return true;
    }
    // Use a coarse tolerance bound via repeated halving of one.
    let p = a.precision().max(b.precision());
    let two = BigFloat::try_from_i64_exact(2, p).unwrap();
    let one = BigFloat::try_from_i64_exact(1, p).unwrap();
    let abs_b = b.abs();
    let mut bound = if abs_b.is_zero() { one.clone() } else { abs_b };
    for _ in 0..bits {
        bound = bound.div(&two, RoundingMode::NearestEven).0;
    }
    matches!(
        abs_diff.partial_cmp(&bound).0,
        Some(Ordering::Less | Ordering::Equal)
    )
}

proptest! {
    /// ln(1) = 0.
    #[test]
    fn ln_one_zero(p in arb_precision()) {
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (r, status) = one.ln(RoundingMode::NearestEven);
        prop_assert!(status.is_ok());
        prop_assert!(r.is_zero());
    }

    /// ln(positive_int) is positive for integers > 1, negative for
    /// integers < 1 (which we don't generate here).
    #[test]
    fn ln_sign_matches_input(n in 2i64..=10_000) {
        let x = BigFloat::try_from_i64_exact(n, 113).unwrap();
        let (r, _) = x.ln(RoundingMode::NearestEven);
        prop_assert!(r.is_sign_positive());
    }

    /// ln(negative) = qNaN + INVALID.
    #[test]
    fn ln_negative_invalid(n in 1i64..=10_000) {
        let neg = BigFloat::try_from_i64_exact(-n, 113).unwrap();
        let (r, status) = neg.ln(RoundingMode::NearestEven);
        prop_assert!(r.is_quiet_nan());
        prop_assert!(status.invalid());
    }

    /// ln(0) = -∞ + DIV_BY_ZERO.
    #[test]
    fn ln_zero_div_by_zero(_dummy in 0u32..1) {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, status) = z.ln(RoundingMode::NearestEven);
        prop_assert!(r.is_infinite());
        prop_assert!(r.is_sign_negative());
        prop_assert!(status.div_by_zero());
    }

    /// ln is monotonic on positive inputs.
    #[test]
    fn ln_is_monotonic(n in 1i64..=1000) {
        let p = 113u32;
        let a = BigFloat::try_from_i64_exact(n, p).unwrap();
        let b = BigFloat::try_from_i64_exact(n + 1, p).unwrap();
        let (la, _) = a.ln(RoundingMode::NearestEven);
        let (lb, _) = b.ln(RoundingMode::NearestEven);
        prop_assert_eq!(lb.partial_cmp(&la).0, Some(Ordering::Greater));
    }

    /// exp(ln(x)) ≈ x for positive x. The round-trip relative
    /// error is bounded by the input's precision (we can't recover
    /// bits lost when storing x at its precision); allow a small
    /// slack for accumulated rounding.
    #[test]
    fn exp_ln_round_trip(x in arb_positive_input()) {
        let (ln_x, _) = x.ln(RoundingMode::NearestEven);
        let (back, _) = ln_x.exp(RoundingMode::NearestEven);
        let tol = x.precision().saturating_sub(8);
        prop_assert!(
            close_within(&back, &x, tol),
            "exp(ln({x})) = {back}",
        );
    }

    /// ln(exp(x)) ≈ x for moderate x.
    #[test]
    fn ln_exp_round_trip(n in -20i64..=20) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (e_x, _) = x.exp(RoundingMode::NearestEven);
        let (back, _) = e_x.ln(RoundingMode::NearestEven);
        prop_assert!(
            close_within(&back, &x, p.saturating_sub(8)),
            "ln(exp({n})) = {back}",
        );
    }

    /// ln(a * b) ≈ ln(a) + ln(b) (additive identity).
    #[test]
    fn ln_product_is_sum(a in 1i64..=100, b in 1i64..=100) {
        let p = 113u32;
        let av = BigFloat::try_from_i64_exact(a, p).unwrap();
        let bv = BigFloat::try_from_i64_exact(b, p).unwrap();
        let (prod, _) = av.mul(&bv, RoundingMode::NearestEven);
        let (ln_prod, _) = prod.ln(RoundingMode::NearestEven);
        let (ln_a, _) = av.ln(RoundingMode::NearestEven);
        let (ln_b, _) = bv.ln(RoundingMode::NearestEven);
        let (sum, _) = ln_a.add(&ln_b, RoundingMode::NearestEven);
        prop_assert!(
            close_within(&ln_prod, &sum, p.saturating_sub(8)),
            "ln({a}*{b}) = {ln_prod}, ln({a})+ln({b}) = {sum}",
        );
    }
}
