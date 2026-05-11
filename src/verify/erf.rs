//! Kani harnesses for [`BigFloat::erf`].
//!
//! `erf(NaN) = NaN`, `erf(±∞) = ±1`, `erf(±0) = ±0`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn erf_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.erf(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn erf_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.erf(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.contains(Status::INVALID));
}

#[kani::proof]
fn erf_pos_inf_is_one() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.erf(RoundingMode::NearestEven);
    let one = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (cmp, _) = r.partial_cmp(&one);
    assert_eq!(cmp, Some(core::cmp::Ordering::Equal));
    assert!(status.is_empty());
}

#[kani::proof]
fn erf_neg_inf_is_neg_one() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.erf(RoundingMode::NearestEven);
    let neg_one = BigFloat::try_from_i64_exact(-1, 53).expect("-1 fits");
    let (cmp, _) = r.partial_cmp(&neg_one);
    assert_eq!(cmp, Some(core::cmp::Ordering::Equal));
    assert!(status.is_empty());
}

#[kani::proof]
fn erf_neg_zero_is_neg_zero() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.erf(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_negative());
    assert!(status.is_empty());
}
