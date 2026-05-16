//! Kani harnesses for [`BigFloat::li`].
//!
//! `li(NaN) = NaN`, `li(0) = 0`, `li(1) = −∞ + DIV_BY_ZERO`,
//! `li(+∞) = +∞`, `li(x < 0) = NaN + INVALID`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn li_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.li(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn li_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.li(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn li_zero_is_zero() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.li(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(status.is_ok());
}

#[kani::proof]
fn li_one_is_neg_inf_div_by_zero() {
    let a = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, status) = a.li(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn li_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.li(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(!r.is_sign_negative());
    assert!(status.is_ok());
}

#[kani::proof]
fn li_neg_finite_is_nan_invalid() {
    let a = BigFloat::try_from_i64_exact(-2, 53).expect("-2 fits");
    let (r, status) = a.li(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}
