//! Kani harnesses for [`BigFloat::erfc`].
//!
//! `erfc(NaN) = NaN`, `erfc(±0) = +1`, `erfc(+∞) = +0`,
//! `erfc(−∞) = +2`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn erfc_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.erfc(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn erfc_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.erfc(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn erfc_pos_zero_is_one() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.erfc(RoundingMode::NearestEven);
    let one = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (cmp, _) = r.partial_cmp(&one);
    assert_eq!(cmp, Some(core::cmp::Ordering::Equal));
    assert!(status.is_ok());
}

#[kani::proof]
fn erfc_pos_inf_is_pos_zero() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.erfc(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

#[kani::proof]
fn erfc_neg_inf_is_two() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.erfc(RoundingMode::NearestEven);
    let two = BigFloat::try_from_i64_exact(2, 53).expect("2 fits");
    let (cmp, _) = r.partial_cmp(&two);
    assert_eq!(cmp, Some(core::cmp::Ordering::Equal));
    assert!(status.is_ok());
}
