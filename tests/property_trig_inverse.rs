//! Property-based tests for the inverse trig family shipped in
//! slice 3g: `asin`, `acos`, `atan`, `atan2`.

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
    // ---- atan ----

    /// atan(0) = 0.
    #[test]
    fn atan_zero(p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.atan(RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
    }

    /// atan is odd.
    #[test]
    fn atan_is_odd(n in -100i64..=100) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let neg_x = BigFloat::try_from_i64_exact(-n, p).unwrap();
        let (a, _) = x.atan(RoundingMode::NearestEven);
        let (b, _) = neg_x.atan(RoundingMode::NearestEven);
        let neg_a = a.negated();
        prop_assert!(close_within(&b, &neg_a, p.saturating_sub(12)));
    }

    /// atan(tan(x)) = x for x ∈ (−π/2, π/2). Restrict to small x
    /// where tan(x) is finite and within reasonable magnitude.
    #[test]
    fn atan_tan_round_trip(n in -1i64..=1) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (t, _) = x.tan(RoundingMode::NearestEven);
        let (back, _) = t.atan(RoundingMode::NearestEven);
        prop_assert!(
            close_within(&back, &x, p.saturating_sub(16)),
            "atan(tan({n})) = {back}",
        );
    }

    // ---- asin ----

    /// asin(0) = 0.
    #[test]
    fn asin_zero(p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.asin(RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
    }

    /// asin is odd.
    #[test]
    fn asin_is_odd(num in 1i64..=9) {
        let p = 113u32;
        let den = BigFloat::try_from_i64_exact(10, p).unwrap();
        let n = BigFloat::try_from_i64_exact(num, p).unwrap();
        let (x, _) = n.div(&den, RoundingMode::NearestEven);
        let neg_x = x.negated();
        let (a, _) = x.asin(RoundingMode::NearestEven);
        let (b, _) = neg_x.asin(RoundingMode::NearestEven);
        let neg_a = a.negated();
        prop_assert!(close_within(&b, &neg_a, p.saturating_sub(12)));
    }

    /// asin(sin(x)) = x for x ∈ [−π/2, π/2].
    #[test]
    fn asin_sin_round_trip(n in -1i64..=1) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (s, _) = x.sin(RoundingMode::NearestEven);
        let (back, _) = s.asin(RoundingMode::NearestEven);
        prop_assert!(close_within(&back, &x, p.saturating_sub(16)));
    }

    // ---- acos ----

    /// acos(cos(x)) = x for x ∈ [0, π]. Restrict to small x.
    #[test]
    fn acos_cos_round_trip(n in 0i64..=3) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (c, _) = x.cos(RoundingMode::NearestEven);
        let (back, _) = c.acos(RoundingMode::NearestEven);
        prop_assert!(
            close_within(&back, &x, p.saturating_sub(16)),
            "acos(cos({n})) = {back}",
        );
    }

    /// asin(x) + acos(x) = π/2. Holds for any x ∈ [−1, 1].
    #[test]
    fn asin_plus_acos_is_pi_over_2(num in -9i64..=9) {
        let p = 113u32;
        let den = BigFloat::try_from_i64_exact(10, p).unwrap();
        let n = BigFloat::try_from_i64_exact(num, p).unwrap();
        let (x, _) = n.div(&den, RoundingMode::NearestEven);
        let (a, _) = x.asin(RoundingMode::NearestEven);
        let (c, _) = x.acos(RoundingMode::NearestEven);
        let (sum, _) = a.add(&c, RoundingMode::NearestEven);
        // π/2 reference: via the public sin → asin path on a known value.
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (pi_over_2, _) = one.asin(RoundingMode::NearestEven);
        prop_assert!(
            close_within(&sum, &pi_over_2, p.saturating_sub(16)),
            "asin({num}/10) + acos({num}/10) = {sum}",
        );
    }

    // ---- atan2 ----

    /// atan2(y, x) = atan(y/x) for x > 0. The dispatch path
    /// matches the single-argument atan in the first/fourth
    /// quadrants.
    #[test]
    fn atan2_matches_atan_for_positive_x(y in -10i64..=10, x in 1i64..=10) {
        let p = 113u32;
        let y_b = BigFloat::try_from_i64_exact(y, p).unwrap();
        let x_b = BigFloat::try_from_i64_exact(x, p).unwrap();
        let (r, _) = y_b.atan2(&x_b, RoundingMode::NearestEven);
        let (ratio, _) = y_b.div(&x_b, RoundingMode::NearestEven);
        let (expected, _) = ratio.atan(RoundingMode::NearestEven);
        prop_assert!(close_within(&r, &expected, p.saturating_sub(16)));
    }

    /// atan2(y, x) for the second quadrant (x < 0, y > 0): result
    /// equals atan(y/x) + π and lies in (π/2, π).
    #[test]
    fn atan2_quadrant_two(y in 1i64..=10, x in 1i64..=10) {
        let p = 113u32;
        let y_b = BigFloat::try_from_i64_exact(y, p).unwrap();
        let neg_x = BigFloat::try_from_i64_exact(-x, p).unwrap();
        let (r, _) = y_b.atan2(&neg_x, RoundingMode::NearestEven);
        // Expected: π − atan(y/x). Use the public surface for π:
        // 2·acos(0).
        let zero = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (pi_over_2, _) = zero.acos(RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let (pi, _) = pi_over_2.mul(&two, RoundingMode::NearestEven);
        let (ratio, _) = y_b.div(&BigFloat::try_from_i64_exact(x, p).unwrap(), RoundingMode::NearestEven);
        let (atan_pos, _) = ratio.atan(RoundingMode::NearestEven);
        let (expected, _) = pi.sub(&atan_pos, RoundingMode::NearestEven);
        prop_assert!(close_within(&r, &expected, p.saturating_sub(16)));
    }
}
