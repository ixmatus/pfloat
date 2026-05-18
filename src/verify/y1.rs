//! Kani harnesses for [`BigFloat::y1`] (Bessel `Y1`, real only for
//! `x > 0`).
//!
//! Same domain table as [`super::y0`]: `Y1(NaN) = NaN`;
//! `Y1(+0) = −∞` raising `DIV_BY_ZERO` (a pole, DLMF 10.8.1);
//! `Y1(−0)`, `Y1(x<0)`, `Y1(−∞)` are `qNaN` + `INVALID` (complex off
//! the positive axis); `Y1(+∞) = +0` (decaying-envelope).

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn y1_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.y1(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn y1_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.y1(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn y1_pos_zero_is_pole() {
    // Y1(+0) = −∞ + DIV_BY_ZERO (a pole, DLMF 10.8.1).
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.y1(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn y1_neg_zero_is_invalid() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.y1(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn y1_negative_arg_is_invalid() {
    let a = BigFloat::try_from_i64_exact(-3, 53).expect("precision >= 1");
    let (r, status) = a.y1(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn y1_pos_inf_is_pos_zero() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, _status) = a.y1(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
}

#[kani::proof]
fn y1_neg_inf_is_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.y1(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}
