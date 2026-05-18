//! Kani harnesses for [`BigFloat::k1`] (modified Bessel `K1`, real
//! only for `x > 0`).
//!
//! `K1(NaN) = NaN`; `K1(+0) = +∞` raising `DIV_BY_ZERO` (a pole,
//! DLMF 10.30.2); `K1(−0)`, `K1(x<0)`, `K1(−∞)` are `qNaN` +
//! `INVALID` (complex off the positive axis); `K1(+∞) = +0` (a
//! genuine exponential-decay limit, `Status::OK`). The Normal
//! series is out of Kani scope.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn k1_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.k1(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn k1_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.k1(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn k1_pos_zero_is_pole() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.k1(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(!r.is_sign_negative());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn k1_neg_zero_is_invalid() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.k1(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn k1_negative_arg_is_invalid() {
    let a = BigFloat::try_from_i64_exact(-3, 53).expect("precision >= 1");
    let (r, status) = a.k1(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn k1_pos_inf_is_pos_zero() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.k1(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
    assert!(!status.invalid());
}

#[kani::proof]
fn k1_neg_inf_is_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.k1(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}
