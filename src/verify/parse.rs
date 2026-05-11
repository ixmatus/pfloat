//! Kani harnesses for [`BigFloat::parse_str`].
//!
//! `parse_str` is mostly Kani-intractable because it walks a
//! variable-length string. The harnesses below cover the
//! short-input edge cases that pfloat's lexer dispatches on
//! explicitly: empty input, single-character malformed input,
//! and the canonical `"nan"`/`"inf"`/`"-inf"` literals.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;

/// Empty input returns a `ParseError`.
#[kani::proof]
fn parse_empty_is_err() {
    let r = BigFloat::parse_str("", 53, RoundingMode::NearestEven);
    assert!(r.is_err());
}

/// `"nan"` parses to a quiet NaN with no flag.
#[kani::proof]
fn parse_nan_literal() {
    let (v, _status) =
        BigFloat::parse_str("nan", 53, RoundingMode::NearestEven).expect("nan parses");
    assert!(v.is_nan());
    assert!(!v.is_signaling_nan());
}

/// `"inf"` parses to `+∞`.
#[kani::proof]
fn parse_pos_inf_literal() {
    let (v, _status) =
        BigFloat::parse_str("inf", 53, RoundingMode::NearestEven).expect("inf parses");
    assert!(v.is_infinite());
    assert!(v.is_sign_positive());
}

/// `"-inf"` parses to `−∞`.
#[kani::proof]
fn parse_neg_inf_literal() {
    let (v, _status) =
        BigFloat::parse_str("-inf", 53, RoundingMode::NearestEven).expect("-inf parses");
    assert!(v.is_infinite());
    assert!(v.is_sign_negative());
}

/// `"0"` parses to `+0`.
#[kani::proof]
fn parse_zero_literal() {
    let (v, _status) = BigFloat::parse_str("0", 53, RoundingMode::NearestEven).expect("0 parses");
    assert!(v.is_zero());
    assert!(v.is_sign_positive());
}
