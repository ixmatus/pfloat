//! Kani harnesses for [`BigFloat::ai`] (Airy `Ai`, entire on ℝ).
//!
//! `Ai(NaN) = NaN`; `Ai(±0) = Ai(0) ≈ 0.355` (finite positive
//! normal); `Ai(+∞) = +0`; `Ai(−∞) = +0` (decaying-envelope
//! convention, ADR-0021).

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn ai_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.ai(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn ai_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.ai(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn ai_pos_zero_is_finite_positive() {
    // Ai(0) = 1/(3^{2/3}Γ(2/3)) ≈ 0.355: finite positive normal.
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, _status) = a.ai(RoundingMode::NearestEven);
    assert!(!r.is_nan());
    assert!(!r.is_infinite());
    assert!(!r.is_zero());
    assert!(!r.is_sign_negative());
}

#[kani::proof]
fn ai_neg_zero_is_finite_positive() {
    // Airy is entire: Ai(−0) = Ai(+0) = Ai(0).
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, _status) = a.ai(RoundingMode::NearestEven);
    assert!(!r.is_nan());
    assert!(!r.is_infinite());
    assert!(!r.is_zero());
    assert!(!r.is_sign_negative());
}

#[kani::proof]
fn ai_pos_inf_is_pos_zero() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, _status) = a.ai(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
}

#[kani::proof]
fn ai_neg_inf_is_pos_zero() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, _status) = a.ai(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
}
