//! Property tests for `Ei` (slice 6m): the `Ei(x) = li(eˣ)`
//! identity and the `Ei(0)` pole.

#![cfg(all(feature = "big", feature = "integrals"))]

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
    /// Ei(x) = li(eˣ) for real x ≠ 0 (DLMF 6.2.8). Composition,
    /// so a few-ULP tolerance.
    #[test]
    fn ei_equals_li_of_exp(num in -30i64..=30, den in 1i64..=8) {
        prop_assume!(num != 0);
        let p = 160u32;
        let n = BigFloat::try_from_i64_exact(num, p).unwrap();
        let d = BigFloat::try_from_i64_exact(den, p).unwrap();
        let (x, _) = n.div(&d, RoundingMode::NearestEven);
        let (ei, _) = x.ei(RoundingMode::NearestEven);
        let (ex, _) = x.exp(RoundingMode::NearestEven);
        let (li_ex, _) = ex.li(RoundingMode::NearestEven);
        prop_assert!(
            close_within(&ei, &li_ex, p - 16),
            "Ei({num}/{den})={ei} vs li(e^x)={li_ex}"
        );
    }

    /// Ei(±0) = −∞ raising DIV_BY_ZERO (a pole).
    #[test]
    fn ei_zero_is_pole(p in prop_oneof![Just(53u32), Just(113u32)]) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, status) = z.ei(RoundingMode::NearestEven);
        prop_assert!(r.is_infinite() && r.is_sign_negative());
        prop_assert!(status.div_by_zero());
    }
}
