//! Kani harnesses for [`BigFloat::asinh`].
//!
//! `asinh(NaN) = NaN`, `asinh(±∞) = ±∞`, `asinh(±0) = ±0`.
//! Domain is the full real line; no INVALID emissions outside
//! the sNaN trigger.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn asinh_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.asinh(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn asinh_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.asinh(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn asinh_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.asinh(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

#[kani::proof]
fn asinh_neg_inf_is_neg_inf() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.asinh(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.is_ok());
}

#[kani::proof]
fn asinh_neg_zero_is_neg_zero() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.asinh(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_negative());
    assert!(status.is_ok());
}
