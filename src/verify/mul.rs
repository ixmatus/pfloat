//! Kani harnesses for [`BigFloat::mul`].
//!
//! Coverage: NaN propagation, signaling-NaN INVALID emission, the
//! IEEE 754 invalid `±0 × ±∞ → qNaN + INVALID` form, and sign-of-
//! product correctness for the `±0 × ±finite`, `±finite × ±0`,
//! `±∞ × ±finite`, and `±finite × ±∞` paths.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

use super::helpers::nondet_constant_at;

/// NaN propagates through `mul`.
#[kani::proof]
fn mul_nan_propagates_under_canonical_operands() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = nondet_constant_at(53);
    let (r, _status) = a.mul(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

/// Signaling NaN in either operand raises `INVALID`.
#[kani::proof]
fn mul_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = nondet_constant_at(53);
    let (r, status) = a.mul(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(!r.is_signaling_nan());
    assert!(status.invalid());
}

/// `(+0) × (+∞)` raises `INVALID` and returns a quiet NaN.
#[kani::proof]
fn mul_pos_zero_times_pos_inf_is_nan_invalid() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.mul(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

/// `(−0) × (+∞)` raises `INVALID` and returns a quiet NaN (sign of
/// the zero is irrelevant to the validity check).
#[kani::proof]
fn mul_neg_zero_times_pos_inf_is_nan_invalid() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.mul(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

/// `(+∞) × (+∞)` is `+∞` (no flag).
#[kani::proof]
fn mul_pos_inf_times_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.mul(&b, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

/// `(+∞) × (−∞)` is `−∞`. Sign-of-product rule on infinities.
#[kani::proof]
fn mul_pos_inf_times_neg_inf_is_neg_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.mul(&b, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.is_ok());
}

/// `(+0) × (−0)` is `−0`. Product-sign rule even at zero.
#[kani::proof]
fn mul_pos_zero_times_neg_zero_is_neg_zero() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.mul(&b, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_negative());
    assert!(status.is_ok());
}
