//! Property tests for `li` (slice 6m): the `li(eˣ) = Ei(x)`
//! identity and the domain rules (`li(0)=0`, `li(1)` pole,
//! `li(x<0)` invalid).

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
    /// li(eˣ) = Ei(x) (the dual of the Ei identity).
    #[test]
    fn li_of_exp_equals_ei(num in -30i64..=30, den in 1i64..=8) {
        prop_assume!(num != 0);
        let p = 160u32;
        let n = BigFloat::try_from_i64_exact(num, p).unwrap();
        let d = BigFloat::try_from_i64_exact(den, p).unwrap();
        let (x, _) = n.div(&d, RoundingMode::NearestEven);
        let (ex, _) = x.exp(RoundingMode::NearestEven);
        let (li_ex, _) = ex.li(RoundingMode::NearestEven);
        let (ei, _) = x.ei(RoundingMode::NearestEven);
        prop_assert!(
            close_within(&li_ex, &ei, p - 16),
            "li(e^({num}/{den}))={li_ex} vs Ei={ei}"
        );
    }

    /// li(0) = 0 (not a pole).
    #[test]
    fn li_zero_is_zero(p in prop_oneof![Just(53u32), Just(113u32)]) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, status) = z.li(RoundingMode::NearestEven);
        prop_assert!(r.is_zero());
        prop_assert!(status.is_ok());
    }

    /// li(1) = −∞ + DIV_BY_ZERO (pole).
    #[test]
    fn li_one_is_pole(p in prop_oneof![Just(53u32), Just(113u32)]) {
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (r, status) = one.li(RoundingMode::NearestEven);
        prop_assert!(r.is_infinite() && r.is_sign_negative());
        prop_assert!(status.div_by_zero());
    }

    /// li(x < 0) = NaN + INVALID.
    #[test]
    fn li_negative_is_nan_invalid(num in 1i64..=1000) {
        let p = 64u32;
        let x = BigFloat::try_from_i64_exact(-num, p).unwrap();
        let (r, status) = x.li(RoundingMode::NearestEven);
        prop_assert!(r.is_nan());
        prop_assert!(status.invalid());
    }
}
