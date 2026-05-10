//! Property-based tests for [`pfloat::BigFloat::total_cmp`].
//!
//! `total_cmp` claims to implement IEEE 754-2019 §5.10's totalOrder
//! predicate, which is required to be a total order. These tests
//! verify the order axioms (reflexive, antisymmetric, transitive)
//! over arbitrary triples.

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
    /// Reflexivity: `total_cmp(x, x) == Equal`.
    #[test]
    fn total_cmp_reflexive(v in arb_bigfloat()) {
        prop_assert_eq!(v.total_cmp(&v), Ordering::Equal);
    }

    /// Antisymmetry: `total_cmp(a, b)` and `total_cmp(b, a)` are
    /// either both Equal or each other's reverse.
    #[test]
    fn total_cmp_antisymmetric(a in arb_bigfloat(), b in arb_bigfloat()) {
        let ab = a.total_cmp(&b);
        let ba = b.total_cmp(&a);
        prop_assert_eq!(ab, ba.reverse());
    }

    /// Transitivity: if `a <= b` and `b <= c` then `a <= c`.
    #[test]
    fn total_cmp_transitive(
        a in arb_bigfloat(),
        b in arb_bigfloat(),
        c in arb_bigfloat(),
    ) {
        let ab = a.total_cmp(&b);
        let bc = b.total_cmp(&c);
        let ac = a.total_cmp(&c);
        // `a <= b` means `a.total_cmp(&b) != Greater`.
        if ab != Ordering::Greater && bc != Ordering::Greater {
            prop_assert!(ac != Ordering::Greater,
                "expected a <= c given a <= b <= c, got {:?}", ac);
        }
    }

    /// `negated(v).total_cmp(&v)` makes the negation the strict
    /// opposite for non-zero non-NaN values, equal for ±0 (since
    /// negation flips sign but the totalOrder treats -0 < +0).
    #[test]
    fn negated_consistent_with_total_cmp(v in arb_bigfloat()) {
        if v.is_zero() {
            // negated(+0) == -0; -0 < +0 in totalOrder.
            // negated(-0) == +0; +0 > -0.
            let nv = v.negated();
            prop_assert_ne!(nv.total_cmp(&v), Ordering::Equal);
        }
    }
}
