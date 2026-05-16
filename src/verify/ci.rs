//! Kani harnesses for [`BigFloat::ci`].
//!
//! `Ci(NaN) = NaN`, `Ci(+0) = −∞ + DIV_BY_ZERO`, `Ci(+∞) = +0`,
//! `Ci(x < 0) = NaN + INVALID` (complex in the reals).

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn ci_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.ci(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn ci_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.ci(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn ci_pos_zero_is_neg_inf_div_by_zero() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.ci(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn ci_pos_inf_is_zero() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.ci(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(status.is_ok());
}

#[kani::proof]
fn ci_neg_finite_is_nan_invalid() {
    let a = BigFloat::try_from_i64_exact(-1, 53).expect("-1 fits");
    let (r, status) = a.ci(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn ci_neg_inf_is_nan_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.ci(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}
