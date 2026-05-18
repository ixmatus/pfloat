//! Kani harnesses for [`BigFloat::y0`] (Bessel `Y0`, real only for
//! `x > 0`).
//!
//! `Y0(NaN) = NaN`; `Y0(+0) = −∞` raising `DIV_BY_ZERO` (a pole,
//! DLMF 10.8.1); `Y0(−0)`, `Y0(x<0)`, `Y0(−∞)` are `qNaN` + `INVALID`
//! (`Y` is complex off the positive axis, the Ci/li convention);
//! `Y0(+∞) = +0` (decaying-envelope, ADR-0021/0023). The Normal
//! series is out of Kani scope; the negative-argument arm returns
//! before the evaluator, so it is exercised here.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn y0_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.y0(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn y0_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.y0(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn y0_pos_zero_is_pole() {
    // Y0(+0) = −∞ + DIV_BY_ZERO (a pole, DLMF 10.8.1).
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.y0(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn y0_neg_zero_is_invalid() {
    // Y0(−0): −0 groups with x < 0 (complex), so qNaN + INVALID.
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.y0(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn y0_negative_arg_is_invalid() {
    // x < 0: Y is complex in the reals → qNaN + INVALID.
    let a = BigFloat::try_from_i64_exact(-3, 53).expect("precision >= 1");
    let (r, status) = a.y0(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn y0_pos_inf_is_pos_zero() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, _status) = a.y0(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
}

#[kani::proof]
fn y0_neg_inf_is_invalid() {
    // −∞ groups with x < 0 (complex) → qNaN + INVALID.
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.y0(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}
