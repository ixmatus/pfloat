//! Kani harnesses for [`BigFloat::asin`].
//!
//! `asin(NaN) = NaN`, `asin(x) → qNaN + INVALID` for `|x| > 1`,
//! `asin(±0) = ±0`. `asin(±1) = ±π/2` is exercised by MPFR
//! differential.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn asin_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.asin(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn asin_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.asin(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn asin_out_of_domain_is_nan_invalid() {
    let a = BigFloat::try_from_i64_exact(2, 53).expect("2 fits");
    let (r, status) = a.asin(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn asin_pos_inf_is_nan_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.asin(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn asin_neg_zero_is_neg_zero() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.asin(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_negative());
    assert!(status.is_ok());
}
