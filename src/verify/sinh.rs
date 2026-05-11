//! Kani harnesses for [`BigFloat::sinh`].
//!
//! `sinh(NaN) = NaN`, `sinh(±∞) = ±∞`, `sinh(±0) = ±0`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn sinh_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.sinh(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn sinh_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.sinh(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.contains(Status::INVALID));
}

#[kani::proof]
fn sinh_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.sinh(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_empty());
}

#[kani::proof]
fn sinh_neg_inf_is_neg_inf() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.sinh(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.is_empty());
}

#[kani::proof]
fn sinh_neg_zero_is_neg_zero() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.sinh(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_negative());
    assert!(status.is_empty());
}
