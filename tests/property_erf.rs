//! Property-based tests for the erf family shipped in slice 4a:
//! `erf`, `erfc`.

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
    /// erf(0) = 0 at any precision.
    #[test]
    fn erf_zero(p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.erf(RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
    }

    /// erfc(0) = 1 at any precision.
    #[test]
    fn erfc_zero(p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.erfc(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        prop_assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    /// erf is odd.
    #[test]
    fn erf_is_odd(num in 1i64..=20, den in 1i64..=10) {
        let p = 113u32;
        let n = BigFloat::try_from_i64_exact(num, p).unwrap();
        let d = BigFloat::try_from_i64_exact(den, p).unwrap();
        let (x, _) = n.div(&d, RoundingMode::NearestEven);
        let neg_x = x.negated();
        let (a, _) = x.erf(RoundingMode::NearestEven);
        let (b, _) = neg_x.erf(RoundingMode::NearestEven);
        let neg_a = a.negated();
        prop_assert!(close_within(&b, &neg_a, p.saturating_sub(12)));
    }

    /// erf is monotonic on positive inputs.
    #[test]
    fn erf_monotonic(n in 1i64..=30) {
        let p = 113u32;
        let a = BigFloat::try_from_i64_exact(n, p).unwrap();
        let b = BigFloat::try_from_i64_exact(n + 1, p).unwrap();
        let (la, _) = a.erf(RoundingMode::NearestEven);
        let (lb, _) = b.erf(RoundingMode::NearestEven);
        // For large n erf(n) = 1 at target precision; equality is fine.
        let ord = lb.partial_cmp(&la).0;
        prop_assert!(matches!(ord, Some(Ordering::Greater | Ordering::Equal)));
    }

    /// erfc is monotonic decreasing on positive inputs.
    #[test]
    fn erfc_monotonic_decreasing(n in 0i64..=10) {
        let p = 113u32;
        let a = BigFloat::try_from_i64_exact(n, p).unwrap();
        let b = BigFloat::try_from_i64_exact(n + 1, p).unwrap();
        let (la, _) = a.erfc(RoundingMode::NearestEven);
        let (lb, _) = b.erfc(RoundingMode::NearestEven);
        let ord = lb.partial_cmp(&la).0;
        prop_assert!(matches!(ord, Some(Ordering::Less | Ordering::Equal)));
    }

    /// erf(x) + erfc(x) = 1.
    #[test]
    fn erf_plus_erfc_is_one(num in -10i64..=10) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(num, p).unwrap();
        let (a, _) = x.erf(RoundingMode::NearestEven);
        let (b, _) = x.erfc(RoundingMode::NearestEven);
        let (sum, _) = a.add(&b, RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        prop_assert!(
            close_within(&sum, &one, p.saturating_sub(16)),
            "erf({num}) + erfc({num}) = {sum}",
        );
    }

    /// erfc(-x) = 2 - erfc(x).
    #[test]
    fn erfc_reflection(num in 1i64..=20) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(num, p).unwrap();
        let neg_x = x.negated();
        let (a, _) = neg_x.erfc(RoundingMode::NearestEven);
        let (b, _) = x.erfc(RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let (expected, _) = two.sub(&b, RoundingMode::NearestEven);
        prop_assert!(
            close_within(&a, &expected, p.saturating_sub(16)),
            "erfc(-{num}) = {a}, 2 − erfc({num}) = {expected}",
        );
    }

    /// erf bounded in (-1, 1) for any finite x.
    #[test]
    fn erf_bounded(num in -50i64..=50) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(num, p).unwrap();
        let (r, _) = x.erf(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let neg_one = BigFloat::try_from_i64_exact(-1, p).unwrap();
        // Allow equality at the extreme — for |x| past the
        // precision threshold, erf rounds to ±1 exactly.
        let upper = r.partial_cmp(&one).0;
        let lower = r.partial_cmp(&neg_one).0;
        prop_assert!(matches!(upper, Some(Ordering::Less | Ordering::Equal)));
        prop_assert!(matches!(lower, Some(Ordering::Greater | Ordering::Equal)));
    }
}
