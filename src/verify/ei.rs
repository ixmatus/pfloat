//! Kani harnesses for [`BigFloat::ei`].
//!
//! `Ei(NaN) = NaN`, `Ei(±0) = −∞ + DIV_BY_ZERO`, `Ei(+∞) = +∞`,
//! `Ei(−∞) = −0`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn ei_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.ei(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn ei_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.ei(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn ei_pos_zero_is_neg_inf_div_by_zero() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.ei(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn ei_neg_zero_is_neg_inf_div_by_zero() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.ei(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn ei_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.ei(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(!r.is_sign_negative());
    assert!(status.is_ok());
}

#[kani::proof]
fn ei_neg_inf_is_zero() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.ei(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(status.is_ok());
}
