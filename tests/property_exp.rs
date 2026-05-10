//! Property-based tests for [`pfloat::BigFloat::exp`].

#![cfg(all(feature = "big", feature = "exp-log"))]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Sign};
use proptest::prelude::*;

fn arb_mode() -> impl Strategy<Value = RoundingMode> {
    prop_oneof![
        Just(RoundingMode::NearestEven),
        Just(RoundingMode::NearestAway),
        Just(RoundingMode::TowardZero),
        Just(RoundingMode::TowardPositive),
        Just(RoundingMode::TowardNegative),
    ]
}

fn arb_small_input() -> impl Strategy<Value = BigFloat> {
    // Integers in [-30, 30] — well within the Taylor series'
    // accuracy at 53-bit precision with 64-bit guard.
    (
        -30i64..=30,
        prop_oneof![Just(53u32), Just(64u32), Just(113u32)],
    )
        .prop_filter_map("exact-fit", |(n, p)| {
            BigFloat::try_from_i64_exact(n, p).ok()
        })
}

/// Returns true when |a - b| <= 2^-bits at any magnitude (very loose
/// tolerance suitable for relative-error checks via subtraction).
fn within_relative(a: &BigFloat, b: &BigFloat, bits: u32) -> bool {
    let (diff, _) = a.sub(b, RoundingMode::NearestEven);
    let abs_diff = diff.abs();
    if abs_diff.is_zero() {
        return true;
    }
    // Reach for a string comparison instead of accessing internal fields.
    // Format both at low precision and compare.
    let p = a.precision().max(b.precision());
    let scaled = abs_diff
        .round_to_precision(p, RoundingMode::NearestEven)
        .unwrap()
        .0;
    // |diff| < 2^-bits iff log2(|diff|) < -bits.
    // Use the cheap "compare with 2^-bits" approach: build the bound
    // and compare via partial_cmp.
    if bits == 0 {
        return true;
    }
    // Build 2^-bits BigFloat: parse_str of "2e-N" is awkward; build via
    // exact integer + round_to_precision tricks instead.
    let one = BigFloat::try_from_i64_exact(1, p).unwrap();
    // Halve `one` `bits` times: this is equivalent to setting its
    // exponent to -bits. But we don't have public exponent access.
    // Use repeated division by 2.
    let two = BigFloat::try_from_i64_exact(2, p).unwrap();
    let mut bound = one;
    let chunks = bits / 32;
    let two_pow_32 = {
        let mut v = BigFloat::try_from_i64_exact(1, p).unwrap();
        for _ in 0..32 {
            v = v.div(&two, RoundingMode::NearestEven).0;
        }
        v
    };
    for _ in 0..chunks {
        bound = bound.mul(&two_pow_32, RoundingMode::NearestEven).0;
    }
    for _ in 0..(bits % 32) {
        bound = bound.div(&two, RoundingMode::NearestEven).0;
    }
    matches!(
        scaled.partial_cmp(&bound).0,
        Some(Ordering::Less | Ordering::Equal)
    )
}

proptest! {
    /// exp(0) = 1 for any precision.
    #[test]
    fn exp_zero_one(p in prop_oneof![Just(53u32), Just(64u32), Just(113u32), Just(256u32)]) {
        let zero = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = zero.exp(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        prop_assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    /// exp(x) is always positive (or +0 for very negative x).
    #[test]
    fn exp_is_nonneg(x in arb_small_input(), mode in arb_mode()) {
        let (r, _) = x.exp(mode);
        if !r.is_nan() {
            prop_assert!(!r.is_sign_negative() || r.is_zero(),
                "exp({x}) = {r} should be non-negative");
        }
    }

    /// exp(qNaN) propagates the qNaN; exp(sNaN) raises INVALID.
    #[test]
    fn nan_handling(_dummy in 0u32..1) {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = q.exp(RoundingMode::NearestEven);
        prop_assert!(r.is_quiet_nan());
        prop_assert!(!status.invalid());

        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r2, status2) = sn.exp(RoundingMode::NearestEven);
        prop_assert!(r2.is_quiet_nan());
        prop_assert!(status2.invalid());
    }

    /// exp(+∞) = +∞, exp(-∞) = +0.
    #[test]
    fn infinity_inputs(_dummy in 0u32..1) {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.exp(RoundingMode::NearestEven);
        prop_assert!(r.is_infinite());
        prop_assert!(r.is_sign_positive());

        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r2, _) = ni.exp(RoundingMode::NearestEven);
        prop_assert!(r2.is_zero());
        prop_assert!(r2.is_sign_positive());
    }

    /// exp(x) is monotonically increasing.
    #[test]
    fn exp_is_monotonic(n in -20i64..20) {
        let p = 113u32;
        let a = BigFloat::try_from_i64_exact(n, p).unwrap();
        let b = BigFloat::try_from_i64_exact(n + 1, p).unwrap();
        let (ea, _) = a.exp(RoundingMode::NearestEven);
        let (eb, _) = b.exp(RoundingMode::NearestEven);
        prop_assert_eq!(eb.partial_cmp(&ea).0, Some(Ordering::Greater),
            "exp({}) <= exp({})", n + 1, n);
    }

    /// exp(-x) * exp(x) ≈ 1 (reciprocal identity).
    /// At 113-bit precision with ≤30 magnitude, the accumulated
    /// error stays comfortably below 2^-90.
    #[test]
    fn exp_negation_reciprocal(n in 1i64..=20) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let neg_x = x.negated();
        let (e_pos, _) = x.exp(RoundingMode::NearestEven);
        let (e_neg, _) = neg_x.exp(RoundingMode::NearestEven);
        let (product, _) = e_pos.mul(&e_neg, RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        prop_assert!(within_relative(&product, &one, 90),
            "exp({n}) * exp(-{n}) = {product}, expected ≈ 1");
    }
}
