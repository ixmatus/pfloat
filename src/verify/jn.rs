//! Kani harnesses for [`BigFloat::jn`] (Bessel `Jn`, integer order,
//! entire on ℝ), including negative order.
//!
//! `Jn(NaN) = NaN`; `J0(±0) = 1`, `Jn(±0) = +0` for `n ≠ 0` (exact,
//! DLMF 10.2.2); `Jn(±∞) = +0` for every order; the negative order
//! `J₋ₙ = (−1)ⁿ Jₙ` reduction (DLMF 10.4.1) is exercised on the
//! special-value arms (the Normal series is out of Kani scope).

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn jn_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.jn(2, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn jn_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.jn(2, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn jn_zero_order_zero_is_one() {
    // J0(±0) = 1 via the jn path (m == 0 arm).
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, _status) = a.jn(0, RoundingMode::NearestEven);
    let one = BigFloat::try_from_i64_exact(1, 53).expect("precision >= 1");
    assert!(r.partial_cmp(&one).0 == Some(core::cmp::Ordering::Equal));
}

#[kani::proof]
fn jn_nonzero_order_zero_is_pos_zero() {
    // Jn(±0) = 0 for n ≠ 0 (DLMF 10.2.2).
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, _status) = a.jn(3, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
}

#[kani::proof]
fn jn_negative_order_zero_is_pos_zero() {
    // J₋₃(±0): the m = |n| = 3 reduction still yields +0.
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, _status) = a.jn(-3, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
}

#[kani::proof]
fn jn_negative_order_inf_is_pos_zero() {
    // Jn(±∞) = +0 for every order, negative order included.
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, _status) = a.jn(-5, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
}
