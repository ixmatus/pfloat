//! Property-based tests for [`pfloat::BigFloat`] classification
//! predicates.
//!
//! These verify the IEEE 754-2019 §6.2 invariants over arbitrary
//! generated values, complementing the unit tests in
//! `src/classify.rs` that cover specific constants.

#![cfg(feature = "big")]

use pfloat::{BigFloat, IeeeClass, Sign};
use proptest::prelude::*;

fn arb_sign() -> impl Strategy<Value = Sign> {
    prop_oneof![Just(Sign::Positive), Just(Sign::Negative)]
}

fn arb_precision() -> impl Strategy<Value = u32> {
    prop_oneof![
        Just(1u32),
        Just(2u32),
        Just(53u32),
        Just(113u32),
        Just(256u32)
    ]
}

/// Generate a `BigFloat` in any of the four kinds the user can build
/// in slice 1a (zero, infinity, NaN, finite-from-i64).
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
    /// Every value falls into exactly one of `is_nan`, `is_infinite`,
    /// `is_finite` (the latter union of zero and normal). The three
    /// predicates partition the value space.
    #[test]
    fn classification_partitions(v in arb_bigfloat()) {
        let n = u32::from(v.is_nan());
        let i = u32::from(v.is_infinite());
        let f = u32::from(v.is_finite());
        prop_assert_eq!(n + i + f, 1);
    }

    /// `is_quiet_nan` and `is_signaling_nan` partition NaNs.
    #[test]
    fn nan_quiet_signaling_partition(v in arb_bigfloat()) {
        if v.is_nan() {
            let q = u32::from(v.is_quiet_nan());
            let s = u32::from(v.is_signaling_nan());
            prop_assert_eq!(q + s, 1);
        } else {
            prop_assert!(!v.is_quiet_nan());
            prop_assert!(!v.is_signaling_nan());
        }
    }

    /// `is_zero` and `is_normal` partition the finite values.
    #[test]
    fn finite_zero_normal_partition(v in arb_bigfloat()) {
        if v.is_finite() {
            let z = u32::from(v.is_zero());
            let n = u32::from(v.is_normal());
            prop_assert_eq!(z + n, 1);
        }
    }

    /// pfloat never produces a subnormal value.
    #[test]
    fn never_subnormal(v in arb_bigfloat()) {
        prop_assert!(!v.is_subnormal());
    }

    /// `is_sign_positive` and `is_sign_negative` partition every value.
    #[test]
    fn sign_partition(v in arb_bigfloat()) {
        let p = u32::from(v.is_sign_positive());
        let n = u32::from(v.is_sign_negative());
        prop_assert_eq!(p + n, 1);
    }

    /// `abs(v)` is always sign-positive.
    #[test]
    fn abs_is_positive(v in arb_bigfloat()) {
        prop_assert!(v.abs().is_sign_positive());
    }

    /// `negated(negated(v)) == v` (involution).
    #[test]
    fn negated_is_involution(v in arb_bigfloat()) {
        prop_assert_eq!(v.negated().negated(), v);
    }

    /// `negated(v)` flips the sign and preserves the kind.
    #[test]
    fn negated_preserves_kind(v in arb_bigfloat()) {
        let n = v.negated();
        prop_assert_eq!(n.is_nan(), v.is_nan());
        prop_assert_eq!(n.is_infinite(), v.is_infinite());
        prop_assert_eq!(n.is_zero(), v.is_zero());
        prop_assert_eq!(n.is_normal(), v.is_normal());
        prop_assert_eq!(n.is_quiet_nan(), v.is_quiet_nan());
        prop_assert_eq!(n.is_signaling_nan(), v.is_signaling_nan());
        prop_assert_ne!(n.is_sign_positive(), v.is_sign_positive());
    }

    /// `copysign(self, src).sign() == src.sign()`.
    #[test]
    fn copysign_takes_sign_from_arg(a in arb_bigfloat(), b in arb_bigfloat()) {
        let r = a.copysign(&b);
        prop_assert_eq!(r.sign(), b.sign());
    }

    /// `copysign(self, src)` preserves the kind of `self`.
    #[test]
    fn copysign_preserves_kind(a in arb_bigfloat(), b in arb_bigfloat()) {
        let r = a.copysign(&b);
        prop_assert_eq!(r.is_nan(), a.is_nan());
        prop_assert_eq!(r.is_infinite(), a.is_infinite());
        prop_assert_eq!(r.is_zero(), a.is_zero());
        prop_assert_eq!(r.is_normal(), a.is_normal());
    }

    /// `ieee_class()` returns one of the eight pfloat-reachable
    /// variants (the two subnormal variants are unreachable).
    #[test]
    fn ieee_class_reachable_variants(v in arb_bigfloat()) {
        let c = v.ieee_class();
        prop_assert!(!matches!(c, IeeeClass::NegativeSubnormal | IeeeClass::PositiveSubnormal));
    }
}
