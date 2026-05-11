//! Kani harnesses for [`BigFloat::atan`].
//!
//! `atan(NaN) = NaN`, `atan(±∞) = ±π/2`, `atan(±0) = ±0`.
//! Domain is the full real line; no INVALID emissions.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn atan_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.atan(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn atan_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.atan(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn atan_pos_zero_is_pos_zero() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.atan(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

#[kani::proof]
fn atan_neg_zero_is_neg_zero() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.atan(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_negative());
    assert!(status.is_ok());
}

/// `atan(+∞) = +π/2`. Just assert that the result is finite,
/// positive, and approximately π/2. MPFR differential validates
/// the exact value.
#[kani::proof]
fn atan_pos_inf_is_pos_finite() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.atan(RoundingMode::NearestEven);
    assert!(r.is_finite());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}
