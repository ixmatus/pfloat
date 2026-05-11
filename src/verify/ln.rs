//! Kani harnesses for [`BigFloat::ln`].
//!
//! Coverage: NaN propagation, signaling-NaN INVALID,
//! `ln(±0) = −∞ + DIV_BY_ZERO`, `ln(negative_finite) → qNaN +
//! INVALID`, `ln(−∞) → qNaN + INVALID`, `ln(+∞) = +∞`,
//! `ln(+1) = +0`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn ln_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.ln(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn ln_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.ln(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.contains(Status::INVALID));
}

#[kani::proof]
fn ln_pos_zero_is_neg_inf_div_by_zero() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.ln(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.contains(Status::DIV_BY_ZERO));
}

#[kani::proof]
fn ln_neg_zero_is_neg_inf_div_by_zero() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.ln(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.contains(Status::DIV_BY_ZERO));
}

#[kani::proof]
fn ln_neg_finite_is_nan_invalid() {
    let a = BigFloat::try_from_i64_exact(-1, 53).expect("-1 fits");
    let (r, status) = a.ln(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.contains(Status::INVALID));
}

#[kani::proof]
fn ln_neg_inf_is_nan_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.ln(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.contains(Status::INVALID));
}

#[kani::proof]
fn ln_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.ln(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_empty());
}
