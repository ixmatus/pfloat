//! Kani harnesses for [`BigFloat::ai_prime`] (Airy `Ai′`, entire
//! on ℝ).
//!
//! `Ai′(NaN) = NaN`; `Ai′(±0) = Ai′(0) ≈ −0.259` (finite negative
//! normal); `Ai′(+∞) = −0`; `Ai′(−∞) = +0` (decaying-envelope
//! convention, ADR-0021).

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn ai_prime_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.ai_prime(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn ai_prime_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.ai_prime(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn ai_prime_pos_zero_is_finite_negative() {
    // Ai′(0) = −1/(3^{1/3}Γ(1/3)) ≈ −0.259: finite negative normal.
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, _status) = a.ai_prime(RoundingMode::NearestEven);
    assert!(!r.is_nan());
    assert!(!r.is_infinite());
    assert!(!r.is_zero());
    assert!(r.is_sign_negative());
}

#[kani::proof]
fn ai_prime_neg_zero_is_finite_negative() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, _status) = a.ai_prime(RoundingMode::NearestEven);
    assert!(!r.is_nan());
    assert!(!r.is_infinite());
    assert!(!r.is_zero());
    assert!(r.is_sign_negative());
}

#[kani::proof]
fn ai_prime_pos_inf_is_neg_zero() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, _status) = a.ai_prime(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_negative());
}

#[kani::proof]
fn ai_prime_neg_inf_is_pos_zero() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, _status) = a.ai_prime(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
}
