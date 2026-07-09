//! Regression guards for the pf-nt21 fmt/parse review-tail findings
//! adjudicated CONFIRMED and fixed:
//!
//! - `boundary-io/type/3` (fmt.rs): `to_decimal_string(0, ..)` on a nonzero
//!   value returned an all-zeros string in RELEASE (guarded only by a debug
//!   assert). `digits` is now clamped to at least 1.
//! - `boundary-io/type/2` (parse.rs): leading-zero stripping via repeated
//!   `Vec::remove(0)` was O(n^2); a multi-megabyte all-zeros input blew the
//!   ADR-0031 parse budget. Now a single count-then-drain pass.
//!
//! Run: `cargo test --release --features std,big,fmt
//! --test regression_nt21_fmt_parse`.

#![cfg(all(feature = "big", feature = "fmt"))]

use pfloat::{BigFloat, RoundingMode};

const NE: RoundingMode = RoundingMode::NearestEven;

#[test]
fn to_decimal_string_zero_digits_is_clamped_not_junk() {
    // Every build (release included) must render at least one significant
    // digit, never the old all-zeros junk that re-parsed to 0.
    for s in ["1.5", "7", "0.001", "-2.5", "9.99e40"] {
        let v = BigFloat::parse_str(s, 53, NE).unwrap().0;
        let rendered = v.to_decimal_string(0, NE);
        let (reparsed, _) = BigFloat::parse_str(&rendered, 53, NE).unwrap();
        assert!(
            !reparsed.is_zero(),
            "{s}: to_decimal_string(0) = {rendered:?} must not render a nonzero value as 0"
        );
        // Clamped to 1 significant digit == the 1-digit rounding of the value.
        let one_digit = v.to_decimal_string(1, NE);
        assert_eq!(rendered, one_digit, "{s}: digits=0 must behave as digits=1");
    }
}

#[test]
fn leading_zeros_parse_linearly_to_zero() {
    // A large all-zeros mantissa parses to 0 quickly (the O(n^2) strip took
    // minutes at 2 MB; the linear strip is milliseconds). Assert
    // correctness plus a generous wall-clock ceiling as a coarse guard
    // against a quadratic regression.
    use std::time::Instant;
    let s = "0".repeat(2_000_000);
    let t = Instant::now();
    let (v, _) = BigFloat::parse_str(&s, 53, NE).unwrap();
    let elapsed = t.elapsed();
    assert!(v.is_zero() && v.is_sign_positive(), "2M zeros -> +0");
    assert!(
        elapsed.as_secs() < 2,
        "parsing 2M leading zeros took {elapsed:?}; the O(n^2) strip would take minutes"
    );
    // Leading zeros before a significant digit are stripped without changing
    // the value.
    let (a, _) = BigFloat::parse_str(&format!("{}5", "0".repeat(10_000)), 53, NE).unwrap();
    let (b, _) = BigFloat::parse_str("5", 53, NE).unwrap();
    assert_eq!(a.total_cmp(&b), core::cmp::Ordering::Equal, "0…05 == 5");
}
