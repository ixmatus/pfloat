//! Kani harnesses for [`BigFloat::digamma`].
//!
//! `digamma(NaN) = NaN`, `digamma(±0) = −∞ + DIV_BY_ZERO`,
//! `digamma(non-positive integer) → −∞ + DIV_BY_ZERO`,
//! `digamma(+∞) → qNaN + INVALID` (per pfloat's current
//! implementation; the asymptotic series is not directly
//! evaluable at +∞).

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn digamma_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.digamma(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn digamma_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.digamma(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.contains(Status::INVALID));
}

#[kani::proof]
fn digamma_pos_zero_is_neg_inf_div_by_zero() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.digamma(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.contains(Status::DIV_BY_ZERO));
}

#[kani::proof]
fn digamma_neg_integer_is_neg_inf_div_by_zero() {
    let a = BigFloat::try_from_i64_exact(-2, 53).expect("-2 fits");
    let (r, status) = a.digamma(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.contains(Status::DIV_BY_ZERO));
}
