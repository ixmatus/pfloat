//! Kani harnesses for [`BigFloat::atanh`].
//!
//! `atanh(NaN) = NaN`, `atanh(x) → qNaN + INVALID` for `|x| > 1`,
//! `atanh(±1) = ±∞ + DIV_BY_ZERO`, `atanh(±0) = ±0`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn atanh_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.atanh(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn atanh_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.atanh(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn atanh_above_one_is_nan_invalid() {
    let a = BigFloat::try_from_i64_exact(2, 53).expect("2 fits");
    let (r, status) = a.atanh(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn atanh_pos_inf_is_nan_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.atanh(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn atanh_pos_one_is_pos_inf_div_by_zero() {
    let a = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, status) = a.atanh(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn atanh_neg_one_is_neg_inf_div_by_zero() {
    let a = BigFloat::try_from_i64_exact(-1, 53).expect("-1 fits");
    let (r, status) = a.atanh(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn atanh_neg_zero_is_neg_zero() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.atanh(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_negative());
    assert!(status.is_ok());
}
