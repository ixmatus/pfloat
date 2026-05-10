//! Property-based tests for [`pfloat::BigFloat`]'s `Display` and
//! `to_decimal_string`, focused on round-trip identity through
//! `parse_str`.

#![cfg(feature = "big")]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode};
use proptest::prelude::*;

fn arb_precision() -> impl Strategy<Value = u32> {
    prop_oneof![Just(53u32), Just(64u32), Just(113u32), Just(256u32)]
}

proptest! {
    /// Round-trip: `parse_str(value.to_string(), p)` recovers the
    /// exact value when `p` is the value's own precision and at
    /// least the round-trip digit count of decimal digits are
    /// emitted (which `Display` does by default).
    #[test]
    fn display_parse_round_trip_integers(n in any::<i64>(), p in arb_precision()) {
        // Need >= 64-bit precision to hold an i64 exactly.
        let target = p.max(64);
        let v = BigFloat::try_from_i64_exact(n, target).unwrap();
        let s = v.to_string();
        let parsed = BigFloat::parse_str(&s, target, RoundingMode::NearestEven).unwrap().0;
        prop_assert_eq!(parsed.partial_cmp(&v).0, Some(Ordering::Equal));
    }

    /// Round-trip with the exact-fit values from parsing a decimal
    /// literal: parse "1.5e3" (or similar small-integer literal),
    /// format, re-parse, and confirm equality.
    #[test]
    fn parse_display_round_trip(
        n in -10_000i64..=10_000i64,
        k in 0u32..=5,
        p in arb_precision(),
    ) {
        let s = format!("{n}e{k}");
        let v = BigFloat::parse_str(&s, p, RoundingMode::NearestEven).unwrap().0;
        let formatted = v.to_string();
        let back = BigFloat::parse_str(&formatted, p, RoundingMode::NearestEven).unwrap().0;
        prop_assert_eq!(back.partial_cmp(&v).0, Some(Ordering::Equal));
    }

    /// Sign in the string matches the value's sign.
    #[test]
    fn sign_in_display(n in 1i64..=1_000_000) {
        let neg = BigFloat::try_from_i64_exact(-n, 53).unwrap();
        let s = neg.to_string();
        prop_assert!(s.starts_with('-'), "negative value should start with '-': {s}");
    }

    /// Zero displays as "0" or "-0".
    #[test]
    fn zero_display(_dummy in 0u32..1) {
        use pfloat::Sign;
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        prop_assert_eq!(pz.to_string(), "0");
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        prop_assert_eq!(nz.to_string(), "-0");
    }
}
