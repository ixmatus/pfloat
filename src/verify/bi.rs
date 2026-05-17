//! Kani harnesses for [`BigFloat::bi`] (Airy `Bi`, entire on ℝ).
//!
//! `Bi(NaN) = NaN`; `Bi(±0) = Bi(0) ≈ 0.615` (finite positive
//! normal); `Bi(+∞) = +∞` (exact limit at an infinite argument);
//! `Bi(−∞) = +0` (decaying-envelope convention, ADR-0021).

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn bi_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.bi(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn bi_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.bi(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn bi_pos_zero_is_finite_positive() {
    // Bi(0) = 1/(3^{1/6}Γ(2/3)) ≈ 0.615: finite positive normal.
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, _status) = a.bi(RoundingMode::NearestEven);
    assert!(!r.is_nan());
    assert!(!r.is_infinite());
    assert!(!r.is_zero());
    assert!(!r.is_sign_negative());
}

#[kani::proof]
fn bi_neg_zero_is_finite_positive() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, _status) = a.bi(RoundingMode::NearestEven);
    assert!(!r.is_nan());
    assert!(!r.is_infinite());
    assert!(!r.is_zero());
    assert!(!r.is_sign_negative());
}

#[kani::proof]
fn bi_pos_inf_is_pos_inf() {
    // Bi grows without bound: Bi(+∞) = +∞ (Status::OK, the
    // exp(+∞)/gamma(+∞) convention).
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.bi(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(!r.is_sign_negative());
    assert!(status.is_ok());
}

#[kani::proof]
fn bi_neg_inf_is_pos_zero() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, _status) = a.bi(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
}
