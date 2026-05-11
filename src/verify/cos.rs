//! Kani harnesses for [`BigFloat::cos`].
//!
//! `cos(NaN) = NaN`, `cos(±∞) → qNaN + INVALID`,
//! `cos(±0) = +1`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn cos_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.cos(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn cos_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.cos(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn cos_pos_inf_is_nan_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.cos(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn cos_pos_zero_is_one() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.cos(RoundingMode::NearestEven);
    let one = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (cmp, _) = r.partial_cmp(&one);
    assert_eq!(cmp, Some(core::cmp::Ordering::Equal));
    assert!(status.is_ok());
}

#[kani::proof]
fn cos_neg_zero_is_one() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.cos(RoundingMode::NearestEven);
    let one = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (cmp, _) = r.partial_cmp(&one);
    assert_eq!(cmp, Some(core::cmp::Ordering::Equal));
    assert!(status.is_ok());
}
