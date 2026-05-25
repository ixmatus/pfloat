//! Property tests for Airy `Bi` (slice 6n). `Bi` is positive and
//! strictly increasing on `[0, ∞)`, with `Bi(0) ≈ 0.615`; plus
//! precision self-consistency. The Wronskian binding `Bi` to the
//! family lives in `property_ai.rs`.

#![cfg(all(feature = "big", feature = "airy"))]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Sign};
use proptest::prelude::*;

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
    #![proptest_config(ProptestConfig::with_cases(16))]

    /// Bi(0) is a finite positive normal value (≈ 0.615).
    #[test]
    fn bi_zero_is_finite_positive(p in prop_oneof![Just(53u32), Just(113u32)]) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.bi(RoundingMode::NearestEven);
        prop_assert!(!r.is_nan() && !r.is_infinite() && !r.is_zero());
        prop_assert!(!r.is_sign_negative());
    }

    /// On [0, ∞) Bi is strictly positive and strictly increasing.
    #[test]
    fn bi_positive_increasing(a in 0i64..=6, b in 7i64..=15) {
        let p = 80u32;
        let xa = BigFloat::try_from_i64_exact(a, p).unwrap();
        let xb = BigFloat::try_from_i64_exact(b, p).unwrap();
        let (fa, _) = xa.bi(RoundingMode::NearestEven);
        let (fb, _) = xb.bi(RoundingMode::NearestEven);
        let zero = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        prop_assert_eq!(fa.partial_cmp(&zero).0, Some(Ordering::Greater));
        prop_assert_eq!(fb.partial_cmp(&fa).0, Some(Ordering::Greater));
    }

    /// Precision self-consistency: Bi at p agrees with Bi at p+48
    /// rounded to p. The argument is an exact dyadic rational
    /// (denominator a power of two) so both precisions evaluate the
    /// *same* real point — a non-dyadic argument would differ in its
    /// low bits between the two precisions, and near a zero of Bi
    /// (the first at x ≈ −1.174) that argument mismatch is amplified
    /// by |Bi′/Bi| into a spurious failure even at the loose
    /// `p − 12` tolerance (the property of the test construction the
    /// pf-ok9 lesson identifies, also encoded in
    /// `property_yn::self_consistent`, `property_ik::self_consistent`,
    /// `property_zeta::self_consistent`,
    /// `property_jn::self_consistent`, and
    /// `property_ai::ai_self_consistent` per ADR-0036).
    #[test]
    fn bi_self_consistent(
        num in -10i64..=10,
        den in prop_oneof![Just(1i64), Just(2), Just(4)],
    ) {
        let p = 80u32;
        let lo_x = {
            let n = BigFloat::try_from_i64_exact(num, p).unwrap();
            let d = BigFloat::try_from_i64_exact(den, p).unwrap();
            n.div(&d, RoundingMode::NearestEven).0
        };
        let hi_x = {
            let n = BigFloat::try_from_i64_exact(num, p + 48).unwrap();
            let d = BigFloat::try_from_i64_exact(den, p + 48).unwrap();
            n.div(&d, RoundingMode::NearestEven).0
        };
        let (lo, _) = lo_x.bi(RoundingMode::NearestEven);
        let (hi_raw, _) = hi_x.bi(RoundingMode::NearestEven);
        let hi = hi_raw
            .round_to_precision(p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        prop_assert!(close_within(&lo, &hi, p - 12));
    }
}
