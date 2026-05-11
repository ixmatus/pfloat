//! Kani harnesses for [`BigFloat`] formatting via
//! [`core::fmt::Display`].
//!
//! Display is variable-length and string-shaped; Kani cannot
//! enumerate the output meaningfully. The harnesses below assert
//! that Display does not panic on the canonical constants and
//! that the output starts with the expected sign / class prefix.

use crate::big::BigFloat;
use crate::sign::Sign;

extern crate alloc;
use alloc::string::ToString;

#[kani::proof]
fn fmt_quiet_nan_starts_with_nan() {
    let v = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let s = v.to_string();
    assert_eq!(s, "nan");
}

#[kani::proof]
fn fmt_neg_quiet_nan_starts_with_minus() {
    let v = BigFloat::try_new_quiet_nan(Sign::Negative, 53, &[]).expect("precision >= 1");
    let s = v.to_string();
    assert_eq!(s, "-nan");
}

#[kani::proof]
fn fmt_pos_inf_is_inf() {
    let v = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let s = v.to_string();
    assert_eq!(s, "inf");
}

#[kani::proof]
fn fmt_neg_inf_is_minus_inf() {
    let v = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let s = v.to_string();
    assert_eq!(s, "-inf");
}

#[kani::proof]
fn fmt_pos_zero_is_zero() {
    let v = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let s = v.to_string();
    assert_eq!(s, "0");
}

#[kani::proof]
fn fmt_neg_zero_is_minus_zero() {
    let v = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let s = v.to_string();
    assert_eq!(s, "-0");
}
