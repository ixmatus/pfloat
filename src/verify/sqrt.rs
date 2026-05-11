//! Kani harnesses for [`BigFloat::sqrt`].
//!
//! Coverage: NaN propagation, signaling-NaN INVALID, the
//! `sqrt(negative_finite)` invalid form, the IEEE 754-2019 §5.4.1
//! `sqrt(-0) = -0` rule, and `sqrt(+∞) = +∞`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

/// NaN propagates through `sqrt`.
#[kani::proof]
fn sqrt_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.sqrt(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

/// Signaling NaN raises `INVALID`.
#[kani::proof]
fn sqrt_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.sqrt(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(!r.is_signaling_nan());
    assert!(status.contains(Status::INVALID));
}

/// `sqrt(−1)` raises `INVALID` and returns a quiet NaN.
#[kani::proof]
fn sqrt_neg_finite_is_nan_invalid() {
    let a = BigFloat::try_from_i64_exact(-1, 53).expect("-1 fits");
    let (r, status) = a.sqrt(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.contains(Status::INVALID));
}

/// `sqrt(−∞)` raises `INVALID` and returns a quiet NaN.
#[kani::proof]
fn sqrt_neg_inf_is_nan_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.sqrt(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.contains(Status::INVALID));
}

/// `sqrt(−0) = −0` per IEEE 754-2019 §5.4.1. Sign is preserved
/// without raising any flag.
#[kani::proof]
fn sqrt_neg_zero_is_neg_zero() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.sqrt(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_negative());
    assert!(status.is_empty());
}

/// `sqrt(+0) = +0`.
#[kani::proof]
fn sqrt_pos_zero_is_pos_zero() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.sqrt(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_positive());
    assert!(status.is_empty());
}

/// `sqrt(+∞) = +∞`.
#[kani::proof]
fn sqrt_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.sqrt(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_empty());
}
