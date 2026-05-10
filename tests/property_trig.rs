//! Property-based tests for the forward trig family shipped in
//! slice 3f: `sin`, `cos`, `tan`.

#![cfg(all(feature = "big", feature = "trig"))]

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
    /// sin(0) = 0 at any precision.
    #[test]
    fn sin_zero(p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.sin(RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
    }

    /// cos(0) = 1 at any precision.
    #[test]
    fn cos_zero(p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.cos(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        prop_assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    /// sin is odd: sin(-x) = -sin(x).
    #[test]
    fn sin_is_odd(n in -10i64..=10) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let neg_x = BigFloat::try_from_i64_exact(-n, p).unwrap();
        let (a, _) = x.sin(RoundingMode::NearestEven);
        let (b, _) = neg_x.sin(RoundingMode::NearestEven);
        let neg_a = a.negated();
        prop_assert!(close_within(&b, &neg_a, p.saturating_sub(12)));
    }

    /// cos is even: cos(-x) = cos(x).
    #[test]
    fn cos_is_even(n in -10i64..=10) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let neg_x = BigFloat::try_from_i64_exact(-n, p).unwrap();
        let (a, _) = x.cos(RoundingMode::NearestEven);
        let (b, _) = neg_x.cos(RoundingMode::NearestEven);
        prop_assert!(close_within(&a, &b, p.saturating_sub(12)));
    }

    /// sin² + cos² = 1. The identity holds with high relative
    /// precision because both summands are non-negative.
    #[test]
    fn pythagorean_identity(n in -10i64..=10) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (s, _) = x.sin(RoundingMode::NearestEven);
        let (c, _) = x.cos(RoundingMode::NearestEven);
        let (s2, _) = s.mul(&s, RoundingMode::NearestEven);
        let (c2, _) = c.mul(&c, RoundingMode::NearestEven);
        let (sum, _) = s2.add(&c2, RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        prop_assert!(
            close_within(&sum, &one, p.saturating_sub(16)),
            "sin²({n}) + cos²({n}) = {sum}",
        );
    }

    /// tan = sin / cos.
    #[test]
    fn tan_matches_ratio(n in -5i64..=5) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (t, _) = x.tan(RoundingMode::NearestEven);
        let (s, _) = x.sin(RoundingMode::NearestEven);
        let (c, _) = x.cos(RoundingMode::NearestEven);
        let (expected, _) = s.div(&c, RoundingMode::NearestEven);
        prop_assert!(close_within(&t, &expected, p.saturating_sub(16)));
    }

    /// tan is odd.
    #[test]
    fn tan_is_odd(n in -5i64..=5) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let neg_x = BigFloat::try_from_i64_exact(-n, p).unwrap();
        let (a, _) = x.tan(RoundingMode::NearestEven);
        let (b, _) = neg_x.tan(RoundingMode::NearestEven);
        let neg_a = a.negated();
        prop_assert!(close_within(&b, &neg_a, p.saturating_sub(16)));
    }

    /// sin(x + 2π) ≈ sin(x) in absolute terms. Periodicity is
    /// gated by both the accuracy of the parsed 2π and the
    /// reduction's behavior on the shifted argument; the absolute
    /// error floor is governed by ULP at unity at target precision.
    #[test]
    fn sin_periodic(n in -3i64..=3) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let pi_str = "3.14159265358979323846264338327950288419716939937510582097";
        let pi = BigFloat::parse_str(pi_str, p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let two_pi = pi.add(&pi, RoundingMode::NearestEven).0;
        let shifted = x.add(&two_pi, RoundingMode::NearestEven).0;
        let (a, _) = x.sin(RoundingMode::NearestEven);
        let (b, _) = shifted.sin(RoundingMode::NearestEven);
        // Absolute tolerance: difference bounded by ULP(1) · slack.
        let (diff, _) = a.sub(&b, RoundingMode::NearestEven);
        let abs_diff = diff.abs();
        // Tolerance bound = 2^-(p-16).
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let mut bound = one;
        for _ in 0..(p - 16) {
            bound = bound.div(&two, RoundingMode::NearestEven).0;
        }
        prop_assert!(
            matches!(
                abs_diff.partial_cmp(&bound).0,
                Some(Ordering::Less | Ordering::Equal)
            ),
            "sin({n}) = {a}, sin({n}+2π) = {b}, |diff| = {abs_diff}",
        );
    }
}
