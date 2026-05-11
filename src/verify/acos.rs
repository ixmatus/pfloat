//! Kani harnesses for [`BigFloat::acos`].
//!
//! `acos(NaN) = NaN`, `acos(x) → qNaN + INVALID` for `|x| > 1`,
//! `acos(±∞) → qNaN + INVALID`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn acos_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.acos(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn acos_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.acos(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn acos_out_of_domain_is_nan_invalid() {
    let a = BigFloat::try_from_i64_exact(2, 53).expect("2 fits");
    let (r, status) = a.acos(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn acos_pos_inf_is_nan_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.acos(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}
