//! Property-based tests for [`pfloat::BigFloat::pow`].

#![cfg(all(feature = "big", feature = "exp-log"))]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Sign};
use proptest::prelude::*;

fn arb_precision() -> impl Strategy<Value = u32> {
    prop_oneof![Just(53u32), Just(64u32), Just(113u32)]
}

/// `|a - b| <= 2^(-bits) · max(|b|, 1)` via repeated halving — same
/// pattern as `property_ln.rs::close_within`.
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
    /// pow(x, 0) = 1 for any finite positive x.
    #[test]
    fn pow_x_zero_is_one(n in 1i64..=1000, p in arb_precision()) {
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let zero = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = x.pow(&zero, RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        prop_assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    /// pow(1, y) = 1 for any finite y.
    #[test]
    fn pow_one_y_is_one(n in -1000i64..=1000, p in arb_precision()) {
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let y = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (r, _) = one.pow(&y, RoundingMode::NearestEven);
        let expected = BigFloat::try_from_i64_exact(1, p).unwrap();
        prop_assert_eq!(r.partial_cmp(&expected).0, Some(Ordering::Equal));
    }

    /// pow(x, 1) = x.
    #[test]
    fn pow_x_one_is_x(n in 1i64..=1000) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (r, _) = x.pow(&one, RoundingMode::NearestEven);
        prop_assert!(
            close_within(&r, &x, p.saturating_sub(8)),
            "pow({n}, 1) = {r}",
        );
    }

    /// pow(x, 2) = x * x.
    #[test]
    fn pow_x_two_is_square(n in 1i64..=100) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let (r, _) = x.pow(&two, RoundingMode::NearestEven);
        let (sq, _) = x.mul(&x, RoundingMode::NearestEven);
        prop_assert!(
            close_within(&r, &sq, p.saturating_sub(8)),
            "pow({n}, 2) = {r}, {n}*{n} = {sq}",
        );
    }

    /// pow(x, -1) ≈ 1 / x for x > 0.
    #[test]
    fn pow_x_neg_one_is_reciprocal(n in 2i64..=1000) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let neg_one = BigFloat::try_from_i64_exact(-1, p).unwrap();
        let (r, _) = x.pow(&neg_one, RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (recip, _) = one.div(&x, RoundingMode::NearestEven);
        prop_assert!(
            close_within(&r, &recip, p.saturating_sub(8)),
            "pow({n}, -1) = {r}, 1/{n} = {recip}",
        );
    }

    /// pow(x, a) · pow(x, b) ≈ pow(x, a + b).
    #[test]
    fn pow_exponent_addition(a in 1i64..=10, b in 1i64..=10) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(3, p).unwrap();
        let av = BigFloat::try_from_i64_exact(a, p).unwrap();
        let bv = BigFloat::try_from_i64_exact(b, p).unwrap();
        let (pa, _) = x.pow(&av, RoundingMode::NearestEven);
        let (pb, _) = x.pow(&bv, RoundingMode::NearestEven);
        let (prod, _) = pa.mul(&pb, RoundingMode::NearestEven);

        let (sum, _) = av.add(&bv, RoundingMode::NearestEven);
        let (p_sum, _) = x.pow(&sum, RoundingMode::NearestEven);

        prop_assert!(
            close_within(&prod, &p_sum, p.saturating_sub(12)),
            "pow(3, {a})·pow(3, {b}) = {prod}, pow(3, {a}+{b}) = {p_sum}",
        );
    }

    /// pow(negative_base, even_integer) is positive.
    #[test]
    fn pow_neg_base_even_int_positive(n in 1i64..=100, k in 1i64..=10) {
        let p = 113u32;
        let neg_x = BigFloat::try_from_i64_exact(-n, p).unwrap();
        let two_k = BigFloat::try_from_i64_exact(2 * k, p).unwrap();
        let (r, _) = neg_x.pow(&two_k, RoundingMode::NearestEven);
        prop_assert!(r.is_sign_positive());
    }

    /// pow(negative_base, odd_integer) is negative (for n > 0).
    #[test]
    fn pow_neg_base_odd_int_negative(n in 1i64..=100, k in 0i64..=5) {
        let p = 113u32;
        let neg_x = BigFloat::try_from_i64_exact(-n, p).unwrap();
        let odd = BigFloat::try_from_i64_exact(2 * k + 1, p).unwrap();
        let (r, _) = neg_x.pow(&odd, RoundingMode::NearestEven);
        prop_assert!(r.is_sign_negative());
    }

    /// pow(negative, non-integer) raises INVALID and returns qNaN.
    #[test]
    fn pow_neg_base_non_int_is_invalid(n in 1i64..=100, num in 1i64..=10, den in 2i64..=10) {
        let p = 113u32;
        // num/den; if it accidentally is an integer (den | num), skip.
        if num % den == 0 {
            return Ok(());
        }
        let neg_x = BigFloat::try_from_i64_exact(-n, p).unwrap();
        let numerator = BigFloat::try_from_i64_exact(num, p).unwrap();
        let denominator = BigFloat::try_from_i64_exact(den, p).unwrap();
        let (y, _) = numerator.div(&denominator, RoundingMode::NearestEven);
        let (r, status) = neg_x.pow(&y, RoundingMode::NearestEven);
        prop_assert!(r.is_quiet_nan());
        prop_assert!(status.invalid());
    }

    /// pow(x, y) for positive integer y matches repeated multiplication.
    #[test]
    fn pow_matches_repeated_mul(n in 2i64..=10, k in 1i64..=8) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let y = BigFloat::try_from_i64_exact(k, p).unwrap();
        let (via_pow, _) = x.pow(&y, RoundingMode::NearestEven);

        let mut via_mul = BigFloat::try_from_i64_exact(1, p).unwrap();
        for _ in 0..k {
            via_mul = via_mul.mul(&x, RoundingMode::NearestEven).0;
        }

        prop_assert!(
            close_within(&via_pow, &via_mul, p.saturating_sub(12)),
            "pow({n}, {k}) = {via_pow}, repeated = {via_mul}",
        );
    }
}
