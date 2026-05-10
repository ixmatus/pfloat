//! Property-based tests for [`pfloat::BigFloat::partial_cmp`].
//!
//! `partial_cmp` returns `(Option<Ordering>, Status)`. The
//! `Option<Ordering>` is `None` exactly when at least one operand is
//! NaN; the [`pfloat::Status`] flag set raises `INVALID` exactly
//! when at least one operand is a signaling NaN.

#![cfg(feature = "big")]

use core::cmp::Ordering;

use pfloat::{BigFloat, Sign};
use proptest::prelude::*;

fn arb_sign() -> impl Strategy<Value = Sign> {
    prop_oneof![Just(Sign::Positive), Just(Sign::Negative)]
}

fn arb_precision() -> impl Strategy<Value = u32> {
    prop_oneof![Just(1u32), Just(53u32), Just(113u32)]
}

fn arb_bigfloat() -> impl Strategy<Value = BigFloat> {
    (arb_sign(), arb_precision(), 0..5u8, any::<i64>()).prop_filter_map(
        "construct BigFloat",
        |(sign, prec, kind, n)| match kind {
            0 => BigFloat::try_new_zero(sign, prec).ok(),
            1 => BigFloat::try_new_infinity(sign, prec).ok(),
            2 => BigFloat::try_new_quiet_nan(sign, prec, &[]).ok(),
            3 => BigFloat::try_new_signaling_nan(sign, prec, &[]).ok(),
            _ => BigFloat::try_from_i64_exact(n, prec).ok(),
        },
    )
}

proptest! {
    /// `partial_cmp` returns `None` iff at least one operand is NaN.
    #[test]
    fn none_iff_nan(a in arb_bigfloat(), b in arb_bigfloat()) {
        let (ord, _status) = a.partial_cmp(&b);
        prop_assert_eq!(ord.is_none(), a.is_nan() || b.is_nan());
    }

    /// `partial_cmp` raises `INVALID` iff at least one operand is a
    /// signaling NaN.
    #[test]
    fn invalid_iff_signaling_nan(a in arb_bigfloat(), b in arb_bigfloat()) {
        let (_ord, status) = a.partial_cmp(&b);
        let expected_invalid = a.is_signaling_nan() || b.is_signaling_nan();
        prop_assert_eq!(status.invalid(), expected_invalid);
    }

    /// On non-NaN inputs, `partial_cmp` agrees with `total_cmp`
    /// except for the ±0 == ±0 numeric-equality rule.
    #[test]
    fn agrees_with_total_cmp_off_nan_off_zero(a in arb_bigfloat(), b in arb_bigfloat()) {
        if !(a.is_nan() || b.is_nan() || (a.is_zero() && b.is_zero())) {
            let (ord, _) = a.partial_cmp(&b);
            prop_assert_eq!(ord, Some(a.total_cmp(&b)));
        }
    }

    /// `partial_cmp(±0, ±0)` is `Some(Equal)` for both ±0/±0
    /// pairings. (totalOrder distinguishes them; numeric does not.)
    #[test]
    fn zeros_compare_equal_numerically(s1 in arb_sign(), s2 in arb_sign()) {
        let p = 53;
        let z1 = BigFloat::try_new_zero(s1, p).unwrap();
        let z2 = BigFloat::try_new_zero(s2, p).unwrap();
        let (ord, status) = z1.partial_cmp(&z2);
        prop_assert_eq!(ord, Some(Ordering::Equal));
        prop_assert!(status.is_ok());
    }

    /// Antisymmetry on the ordered values: if `a < b`, then `b > a`.
    #[test]
    fn antisymmetric_on_ordered(a in arb_bigfloat(), b in arb_bigfloat()) {
        let (ab, _) = a.partial_cmp(&b);
        let (ba, _) = b.partial_cmp(&a);
        match (ab, ba) {
            (Some(x), Some(y)) => prop_assert_eq!(x, y.reverse()),
            (None, None) => {} // both NaN-flavored, fine
            (Some(_), None) | (None, Some(_)) => {
                prop_assert!(false, "asymmetric NaN handling: ab={:?} ba={:?}", ab, ba);
            }
        }
    }
}
