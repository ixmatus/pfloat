//! Property-based tests for `agm` shipped in slice 6l.

#![cfg(all(feature = "big", feature = "agm"))]

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
    /// agm(x, x) = x at any precision.
    #[test]
    fn agm_fixed_point(x in 1i64..=1000, p in arb_precision()) {
        let v = BigFloat::try_from_i64_exact(x, p).unwrap();
        let (r, _) = v.agm(&v, RoundingMode::NearestEven);
        prop_assert_eq!(r.partial_cmp(&v).0, Some(Ordering::Equal));
    }

    /// agm(0, x) = 0 at any precision.
    #[test]
    fn agm_zero_kills(x in 1i64..=1000, p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let v = BigFloat::try_from_i64_exact(x, p).unwrap();
        let (r, _) = z.agm(&v, RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
        let (r2, _) = v.agm(&z, RoundingMode::NearestEven);
        prop_assert!(r2.is_zero());
    }

    /// agm is symmetric: agm(a, b) = agm(b, a).
    #[test]
    fn agm_is_symmetric(a in 1i64..=1000, b in 1i64..=1000) {
        let p = 113u32;
        let av = BigFloat::try_from_i64_exact(a, p).unwrap();
        let bv = BigFloat::try_from_i64_exact(b, p).unwrap();
        let (ab, _) = av.agm(&bv, RoundingMode::NearestEven);
        let (ba, _) = bv.agm(&av, RoundingMode::NearestEven);
        prop_assert!(close_within(&ab, &ba, p.saturating_sub(8)));
    }

    /// Sandwich: min(a, b) <= agm(a, b) <= max(a, b).
    #[test]
    fn agm_sandwiched(a in 1i64..=1000, b in 1i64..=1000) {
        let p = 113u32;
        let av = BigFloat::try_from_i64_exact(a, p).unwrap();
        let bv = BigFloat::try_from_i64_exact(b, p).unwrap();
        let (r, _) = av.agm(&bv, RoundingMode::NearestEven);
        let lo = if a <= b { &av } else { &bv };
        let hi = if a <= b { &bv } else { &av };
        let lo_ord = r.partial_cmp(lo).0;
        let hi_ord = r.partial_cmp(hi).0;
        prop_assert!(matches!(lo_ord, Some(Ordering::Greater | Ordering::Equal)));
        prop_assert!(matches!(hi_ord, Some(Ordering::Less | Ordering::Equal)));
    }

    /// Step invariance: agm(a, b) = agm((a + b) / 2, sqrt(a · b)).
    #[test]
    fn agm_step_invariance(a in 1i64..=200, b in 1i64..=200) {
        let p = 113u32;
        let av = BigFloat::try_from_i64_exact(a, p).unwrap();
        let bv = BigFloat::try_from_i64_exact(b, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let (sum, _) = av.add(&bv, RoundingMode::NearestEven);
        let (am, _) = sum.div(&two, RoundingMode::NearestEven);
        let (prod, _) = av.mul(&bv, RoundingMode::NearestEven);
        let (gm, _) = prod.sqrt(RoundingMode::NearestEven);
        let (direct, _) = av.agm(&bv, RoundingMode::NearestEven);
        let (one_step, _) = am.agm(&gm, RoundingMode::NearestEven);
        prop_assert!(close_within(&direct, &one_step, p.saturating_sub(16)));
    }
}
