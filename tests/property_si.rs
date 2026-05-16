//! Property tests for `Si` (slice 6m): oddness and `Si(0) = 0`.

#![cfg(all(feature = "big", feature = "integrals"))]

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
    /// Si(0) = 0 at any precision.
    #[test]
    fn si_zero(p in arb_precision()) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.si(RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
    }

    /// Si is odd: Si(−x) = −Si(x).
    #[test]
    fn si_is_odd(num in 1i64..=40, den in 1i64..=8) {
        let p = 128u32;
        let n = BigFloat::try_from_i64_exact(num, p).unwrap();
        let d = BigFloat::try_from_i64_exact(den, p).unwrap();
        let (x, _) = n.div(&d, RoundingMode::NearestEven);
        let (a, _) = x.si(RoundingMode::NearestEven);
        let (b, _) = x.negated().si(RoundingMode::NearestEven);
        prop_assert!(close_within(&b, &a.negated(), p - 12));
    }
}
