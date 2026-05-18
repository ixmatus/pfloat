//! Kani harnesses for [`BigFloat::in_`] (modified Bessel `Iₙ`,
//! integer order, entire), including negative order.
//!
//! `Iₙ(NaN) = NaN`; `Iₙ(±0) = 0` for `n ≠ 0` (exact); `Iₙ(+∞) =
//! +∞`; the order parity `I₋ₙ(x) = Iₙ(x)` is **even with no sign**
//! (DLMF 10.27.1, unlike `J`/`Y`'s `(−1)ⁿ`), exercised on the
//! special-value arms. A negative argument is finite and never
//! `INVALID` (`I` is entire). The Normal series is out of Kani
//! scope.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn in_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.in_(2, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn in_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.in_(2, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn in_pos_zero_is_zero() {
    // Iₙ(+0) = 0 for n ≠ 0 (DLMF 10.30.1).
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.in_(3, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!status.invalid());
}

#[kani::proof]
fn in_negative_order_pos_zero_is_zero() {
    // I₋₃(+0) = I₃(+0) = 0 (order parity, no sign).
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.in_(-3, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!status.invalid());
}

#[kani::proof]
fn in_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.in_(5, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(!r.is_sign_negative());
    assert!(!status.invalid());
}

#[kani::proof]
fn in_negative_order_pos_inf_is_pos_inf() {
    // I₋₅(+∞) = I₅(+∞) = +∞ (order parity, no sign).
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.in_(-5, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(!r.is_sign_negative());
    assert!(!status.invalid());
}

#[kani::proof]
fn in_negative_arg_is_not_invalid() {
    // I is entire: a negative argument is finite, never INVALID
    // (the argument parity Iₙ(−x) = (−1)ⁿ Iₙ(x), not a domain
    // error — the opposite of K's complex-off-axis convention).
    let a = BigFloat::try_from_i64_exact(-3, 53).expect("precision >= 1");
    let (r, status) = a.in_(2, RoundingMode::NearestEven);
    assert!(!r.is_nan());
    assert!(!status.invalid());
}
