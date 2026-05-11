//! Kani harnesses for [`BigFloat::atan2`].
//!
//! The full IEEE 754-2019 §9.2.1 atan2 special-case table is
//! large; this slice covers NaN propagation, signaling-NaN INVALID,
//! and a handful of canonical quadrant edge cases. MPFR
//! differential exercises the dispatch tree on integer pairs.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

#[kani::proof]
fn atan2_nan_in_y_propagates() {
    let y = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let x = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, _status) = y.atan2(&x, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn atan2_signaling_nan_raises_invalid() {
    let y = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let x = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, status) = y.atan2(&x, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

/// `atan2(+0, +1) = +0`.
#[kani::proof]
fn atan2_pos_zero_pos_one_is_pos_zero() {
    let y = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let x = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, status) = y.atan2(&x, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

/// `atan2(-0, +1) = -0`.
#[kani::proof]
fn atan2_neg_zero_pos_one_is_neg_zero() {
    let y = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let x = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, status) = y.atan2(&x, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_negative());
    assert!(status.is_ok());
}
