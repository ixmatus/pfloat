//! Kani harnesses for [`BigFloat::log1p`].
//!
//! `log1p(x) = ln(1 + x)`. Special cases:
//! `log1p(NaN) = NaN`, `log1p(−1) = −∞ + DIV_BY_ZERO`,
//! `log1p(+∞) = +∞`, `log1p(+0) = +0`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn log1p_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.log1p(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn log1p_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.log1p(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.contains(Status::INVALID));
}

#[kani::proof]
fn log1p_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.log1p(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_empty());
}

#[kani::proof]
fn log1p_neg_inf_is_nan_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.log1p(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.contains(Status::INVALID));
}
