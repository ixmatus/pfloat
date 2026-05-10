//! Property-based tests for [`pfloat::BigFloat::parse_str`].

#![cfg(feature = "big")]

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

fn arb_precision() -> impl Strategy<Value = u32> {
    prop_oneof![Just(53u32), Just(64u32), Just(113u32), Just(256u32)]
}

proptest! {
    /// Parsing the decimal representation of an `i64` at sufficient
    /// precision matches the exact integer constructor.
    #[test]
    fn parse_integer_matches_exact(n in any::<i64>(), p in arb_precision()) {
        let s = n.to_string();
        // Need enough precision to hold an arbitrary i64 exactly: 64 bits.
        let target = p.max(64);
        let (parsed, status) = BigFloat::parse_str(&s, target, RoundingMode::NearestEven).unwrap();
        prop_assert!(status.is_ok());
        let exact = BigFloat::try_from_i64_exact(n, target).unwrap();
        prop_assert_eq!(parsed.partial_cmp(&exact).0, Some(Ordering::Equal));
    }

    /// Parsing `"<n>e<k>"` equals parsing `"<n>" * 10^k` (when the
    /// shifted value fits in i64 exactly).
    #[test]
    fn parse_exponent_shift(n in -1000i64..=1000i64, k in 0u32..=10) {
        let s = format!("{n}e{k}");
        let p = 256u32;
        let (parsed, _) = BigFloat::parse_str(&s, p, RoundingMode::NearestEven).unwrap();
        // n × 10^k as an integer.
        let shifted: i64 = n.saturating_mul(10i64.saturating_pow(k));
        if shifted == i64::MAX {
            return Ok(()); // saturated; skip
        }
        let expected = BigFloat::try_from_i64_exact(shifted, p).unwrap();
        prop_assert_eq!(parsed.partial_cmp(&expected).0, Some(Ordering::Equal));
    }

    /// Parsing under any rounding mode gives a result whose
    /// precision matches the requested precision.
    #[test]
    fn parsed_precision_matches_requested(
        n in any::<i64>(),
        p in arb_precision(),
        mode in arb_mode(),
    ) {
        let s = n.to_string();
        let (parsed, _) = BigFloat::parse_str(&s, p, mode).unwrap();
        prop_assert_eq!(parsed.precision(), p);
    }

    /// Parsing a signed number preserves the sign (for non-zero).
    #[test]
    fn parse_sign_preserved(n in 1i64..=1_000_000_i64) {
        let s = format!("-{n}");
        let (parsed, _) = BigFloat::parse_str(&s, 53, RoundingMode::NearestEven).unwrap();
        prop_assert!(parsed.is_sign_negative());
    }

    /// `nan` / `Inf` / `-Inf` parse to the right kind.
    #[test]
    fn special_tokens(_dummy in 0u32..1) {
        let (n, _) = BigFloat::parse_str("nan", 53, RoundingMode::NearestEven).unwrap();
        prop_assert!(n.is_quiet_nan());
        let (pi, _) = BigFloat::parse_str("inf", 53, RoundingMode::NearestEven).unwrap();
        prop_assert!(pi.is_infinite());
        prop_assert_eq!(pi.sign(), Sign::Positive);
        let (ni, _) = BigFloat::parse_str("-inf", 53, RoundingMode::NearestEven).unwrap();
        prop_assert!(ni.is_infinite());
        prop_assert_eq!(ni.sign(), Sign::Negative);
    }
}
