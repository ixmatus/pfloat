//! Kani harnesses for [`BigFloat::yn`] (Bessel `Yn`, integer order,
//! real only for `x > 0`), including negative order.
//!
//! `Yn(NaN) = NaN`; `Yn(+0) = −∞` raising `DIV_BY_ZERO` (a pole,
//! DLMF 10.8.1) for every order; `Yn(−0)`, `Yn(x<0)`, `Yn(−∞)` are
//! `qNaN` + `INVALID` (complex off the positive axis); `Yn(+∞) = +0`
//! for every order. The negative-order reduction `Y₋ₙ = (−1)ⁿ Yₙ`
//! (DLMF 10.4.1) is exercised on the special-value arms (the Normal
//! series is out of Kani scope).

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn yn_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.yn(2, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn yn_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.yn(2, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn yn_pos_zero_is_pole() {
    // Yn(+0) = −∞ + DIV_BY_ZERO for every order (DLMF 10.8.1).
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.yn(3, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn yn_negative_order_pos_zero_is_pole() {
    // Y₋₃(+0): the m = |n| = 3 reduction still yields the pole.
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.yn(-3, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn yn_neg_zero_is_invalid() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.yn(3, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn yn_negative_arg_is_invalid() {
    let a = BigFloat::try_from_i64_exact(-3, 53).expect("precision >= 1");
    let (r, status) = a.yn(2, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn yn_pos_inf_is_pos_zero() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, _status) = a.yn(5, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
}

#[kani::proof]
fn yn_negative_order_inf_is_pos_zero() {
    // Yn(+∞) = +0 for every order, negative order included.
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, _status) = a.yn(-5, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
}

#[kani::proof]
fn yn_neg_inf_is_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.yn(2, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}
