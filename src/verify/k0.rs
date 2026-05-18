//! Kani harnesses for [`BigFloat::k0`] (modified Bessel `K0`, real
//! only for `x > 0`).
//!
//! `K0(NaN) = NaN`; `K0(+0) = +∞` raising `DIV_BY_ZERO` (a pole,
//! DLMF 10.30.2/10.30.3; **positive** ∞, the opposite of
//! `Y0(+0) = −∞`); `K0(−0)`, `K0(x<0)`, `K0(−∞)` are `qNaN` +
//! `INVALID` (`K` is complex off the positive axis, the Ci/li
//! convention); `K0(+∞) = +0` (a genuine exponential-decay limit,
//! `Status::OK`, not the decaying-envelope convention). The Normal
//! series is out of Kani scope; the negative-argument arm returns
//! before the evaluator.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn k0_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.k0(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn k0_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.k0(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn k0_pos_zero_is_pole() {
    // K0(+0) = +∞ + DIV_BY_ZERO (a pole; +∞, opposite of Y0).
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.k0(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(!r.is_sign_negative());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn k0_neg_zero_is_invalid() {
    // −0 groups with x < 0 (complex) → qNaN + INVALID.
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.k0(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn k0_negative_arg_is_invalid() {
    let a = BigFloat::try_from_i64_exact(-3, 53).expect("precision >= 1");
    let (r, status) = a.k0(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn k0_pos_inf_is_pos_zero() {
    // K0(+∞) = +0, the genuine exponential-decay limit, Status::OK.
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.k0(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
    assert!(!status.invalid());
}

#[kani::proof]
fn k0_neg_inf_is_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.k0(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}
