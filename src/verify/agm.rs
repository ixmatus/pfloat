//! Kani harnesses for [`BigFloat::agm`].
//!
//! `agm(NaN, _)` and `agm(_, NaN)` propagate NaN. Signaling NaN
//! raises INVALID. `agm(±0, ±0) = +0`. `agm(neg, _)` raises INVALID
//! and returns qNaN. `agm(+∞, +∞) = +∞`. `agm(+∞, +0)` does not
//! converge: qNaN + INVALID.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn agm_nan_propagates_lhs() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, _status) = a.agm(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn agm_nan_propagates_rhs() {
    let a = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let b = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.agm(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn agm_snan_lhs_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, status) = a.agm(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn agm_snan_rhs_raises_invalid() {
    let a = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let b = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.agm(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn agm_negative_finite_is_invalid() {
    let a = BigFloat::try_from_i64_exact(-1, 53).expect("-1 fits");
    let b = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, status) = a.agm(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn agm_pos_inf_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.agm(&a, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

#[kani::proof]
fn agm_pos_inf_pos_zero_is_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let z = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.agm(&z, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}
