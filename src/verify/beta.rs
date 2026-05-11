//! Kani harnesses for [`BigFloat::beta`].
//!
//! `beta(a, b)` accepts only positive finite `a` and `b` in
//! pfloat's current scope (slice 4c). NaN propagation in either
//! operand, signaling-NaN INVALID, and `INVALID` for any
//! non-positive `a` or `b`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn beta_nan_in_a_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, _status) = a.beta(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn beta_nan_in_b_propagates() {
    let a = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let b = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.beta(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn beta_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, status) = a.beta(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn beta_zero_a_is_nan_invalid() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, status) = a.beta(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn beta_neg_a_is_nan_invalid() {
    let a = BigFloat::try_from_i64_exact(-1, 53).expect("-1 fits");
    let b = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, status) = a.beta(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}
