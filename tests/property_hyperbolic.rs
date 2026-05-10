//! Property-based tests for the hyperbolic family shipped in slice
//! 3e: `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`.

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
    // ---- sinh ----

    /// sinh(0) = 0 at any precision.
    #[test]
    fn sinh_zero(p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.sinh(RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
    }

    /// sinh is odd: sinh(-x) = -sinh(x).
    #[test]
    fn sinh_is_odd(n in -10i64..=10) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let neg_x = BigFloat::try_from_i64_exact(-n, p).unwrap();
        let (sa, _) = x.sinh(RoundingMode::NearestEven);
        let (sb, _) = neg_x.sinh(RoundingMode::NearestEven);
        let neg_sa = sa.negated();
        prop_assert!(close_within(&sb, &neg_sa, p.saturating_sub(12)));
    }

    // ---- cosh ----

    /// cosh(0) = 1.
    #[test]
    fn cosh_zero_is_one(p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.cosh(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        prop_assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    /// cosh is even: cosh(x) = cosh(-x).
    #[test]
    fn cosh_is_even(n in -10i64..=10) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let neg_x = BigFloat::try_from_i64_exact(-n, p).unwrap();
        let (a, _) = x.cosh(RoundingMode::NearestEven);
        let (b, _) = neg_x.cosh(RoundingMode::NearestEven);
        prop_assert!(close_within(&a, &b, p.saturating_sub(12)));
    }

    /// cosh²(x) - sinh²(x) = 1. The subtraction cancels ~2·|x|/ln(2)
    /// bits because cosh²(x) and sinh²(x) both scale like e^(2|x|)/4,
    /// matching in their leading bits. Restrict to small |x| so the
    /// cancellation budget stays inside the target precision.
    #[test]
    fn pythagorean_identity(n in -3i64..=3) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (c, _) = x.cosh(RoundingMode::NearestEven);
        let (s, _) = x.sinh(RoundingMode::NearestEven);
        let (c2, _) = c.mul(&c, RoundingMode::NearestEven);
        let (s2, _) = s.mul(&s, RoundingMode::NearestEven);
        let (diff, _) = c2.sub(&s2, RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        prop_assert!(
            close_within(&diff, &one, p.saturating_sub(24)),
            "cosh²({n}) - sinh²({n}) = {diff}",
        );
    }

    // ---- tanh ----

    /// tanh(0) = 0.
    #[test]
    fn tanh_zero(p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.tanh(RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
    }

    /// tanh is odd: tanh(-x) = -tanh(x).
    #[test]
    fn tanh_is_odd(n in -10i64..=10) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let neg_x = BigFloat::try_from_i64_exact(-n, p).unwrap();
        let (a, _) = x.tanh(RoundingMode::NearestEven);
        let (b, _) = neg_x.tanh(RoundingMode::NearestEven);
        let neg_a = a.negated();
        prop_assert!(close_within(&b, &neg_a, p.saturating_sub(12)));
    }

    /// tanh = sinh / cosh.
    #[test]
    fn tanh_matches_ratio(n in -8i64..=8) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (t, _) = x.tanh(RoundingMode::NearestEven);
        let (s, _) = x.sinh(RoundingMode::NearestEven);
        let (c, _) = x.cosh(RoundingMode::NearestEven);
        let (expected, _) = s.div(&c, RoundingMode::NearestEven);
        prop_assert!(close_within(&t, &expected, p.saturating_sub(16)));
    }

    // ---- asinh ----

    /// asinh(sinh(x)) = x for moderate x.
    #[test]
    fn asinh_sinh_round_trip(n in -5i64..=5) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (s, _) = x.sinh(RoundingMode::NearestEven);
        let (back, _) = s.asinh(RoundingMode::NearestEven);
        prop_assert!(
            close_within(&back, &x, p.saturating_sub(16)),
            "asinh(sinh({n})) = {back}",
        );
    }

    /// asinh is odd.
    #[test]
    fn asinh_is_odd(n in -50i64..=50) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let neg_x = BigFloat::try_from_i64_exact(-n, p).unwrap();
        let (a, _) = x.asinh(RoundingMode::NearestEven);
        let (b, _) = neg_x.asinh(RoundingMode::NearestEven);
        let neg_a = a.negated();
        prop_assert!(close_within(&b, &neg_a, p.saturating_sub(12)));
    }

    // ---- acosh ----

    /// acosh(cosh(x)) = |x| for x ≥ 0.
    #[test]
    fn acosh_cosh_round_trip(n in 0i64..=5) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (c, _) = x.cosh(RoundingMode::NearestEven);
        let (back, _) = c.acosh(RoundingMode::NearestEven);
        prop_assert!(
            close_within(&back, &x, p.saturating_sub(16)),
            "acosh(cosh({n})) = {back}",
        );
    }

    /// acosh(1) = 0 at any precision.
    #[test]
    fn acosh_one_zero(p in arb_precision()) {
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (r, _) = one.acosh(RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
    }

    // ---- atanh ----

    /// atanh(tanh(x)) = x for moderate x (avoid x near edge where
    /// tanh(x) is too close to ±1 for inverse to recover bits).
    #[test]
    fn atanh_tanh_round_trip(n in -3i64..=3) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let (t, _) = x.tanh(RoundingMode::NearestEven);
        let (back, _) = t.atanh(RoundingMode::NearestEven);
        prop_assert!(
            close_within(&back, &x, p.saturating_sub(16)),
            "atanh(tanh({n})) = {back}",
        );
    }

    /// atanh is odd.
    #[test]
    fn atanh_is_odd(num in 1i64..=9, den in 10i64..=10) {
        // x = num/den with |x| < 1.
        let p = 113u32;
        let n = BigFloat::try_from_i64_exact(num, p).unwrap();
        let d = BigFloat::try_from_i64_exact(den, p).unwrap();
        let (x, _) = n.div(&d, RoundingMode::NearestEven);
        let neg_x = x.negated();
        let (a, _) = x.atanh(RoundingMode::NearestEven);
        let (b, _) = neg_x.atanh(RoundingMode::NearestEven);
        let neg_a = a.negated();
        prop_assert!(close_within(&b, &neg_a, p.saturating_sub(12)));
    }
}
