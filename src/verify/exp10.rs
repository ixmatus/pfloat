//! Kani harnesses for [`BigFloat::exp10`].
//!
//! `exp10(x) = 10^x`. Special cases match `exp`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn exp10_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.exp10(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn exp10_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.exp10(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.contains(Status::INVALID));
}

#[kani::proof]
fn exp10_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.exp10(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_empty());
}

#[kani::proof]
fn exp10_neg_inf_is_pos_zero() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.exp10(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_positive());
    assert!(status.is_empty());
}
