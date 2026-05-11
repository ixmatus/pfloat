//! Kani harnesses for [`BigFloat::lgamma`].
//!
//! `lgamma(NaN) = NaN`, `lgamma(±0) = +∞ + DIV_BY_ZERO`,
//! `lgamma(+∞) = +∞`, `lgamma(non-positive integer) → +∞ +
//! DIV_BY_ZERO`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn lgamma_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.lgamma(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn lgamma_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.lgamma(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn lgamma_pos_zero_is_pos_inf_div_by_zero() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.lgamma(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn lgamma_neg_zero_is_pos_inf_div_by_zero() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.lgamma(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn lgamma_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.lgamma(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}
