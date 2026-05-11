//! Kani harnesses for [`BigFloat::cosh`].
//!
//! `cosh(NaN) = NaN`, `cosh(±∞) = +∞`, `cosh(±0) = +1`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn cosh_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.cosh(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn cosh_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.cosh(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn cosh_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.cosh(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

#[kani::proof]
fn cosh_neg_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.cosh(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

#[kani::proof]
fn cosh_pos_zero_is_one() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.cosh(RoundingMode::NearestEven);
    let one = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (cmp, _) = r.partial_cmp(&one);
    assert_eq!(cmp, Some(core::cmp::Ordering::Equal));
    assert!(status.is_ok());
}
