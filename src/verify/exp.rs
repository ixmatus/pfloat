//! Kani harnesses for [`BigFloat::exp`].
//!
//! Coverage: NaN propagation, signaling-NaN INVALID,
//! `exp(+∞) = +∞`, `exp(−∞) = +0`, and `exp(±0) = +1`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn exp_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.exp(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn exp_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.exp(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(!r.is_signaling_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn exp_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.exp(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

#[kani::proof]
fn exp_neg_inf_is_pos_zero() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.exp(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}
