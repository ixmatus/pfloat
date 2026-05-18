//! Kani harnesses for [`BigFloat::kn`] (modified Bessel `Kₙ`,
//! integer order, real only for `x > 0`), including negative order.
//!
//! `Kₙ(NaN) = NaN`; `Kₙ(+0) = +∞` raising `DIV_BY_ZERO` (a pole,
//! DLMF 10.30.2) for every order; `Kₙ(−0)`, `Kₙ(x<0)`, `Kₙ(−∞)` are
//! `qNaN` + `INVALID` (complex off the positive axis); `Kₙ(+∞) =
//! +0` (a genuine exponential-decay limit). The order parity
//! `K₋ₙ(x) = Kₙ(x)` is **even with no sign** (DLMF 10.27.3, unlike
//! `Y`'s `(−1)ⁿ`), exercised on the special-value arms. The Normal
//! series is out of Kani scope.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn kn_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.kn(2, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn kn_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.kn(2, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn kn_pos_zero_is_pole() {
    // Kₙ(+0) = +∞ + DIV_BY_ZERO for every order (DLMF 10.30.2).
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.kn(3, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(!r.is_sign_negative());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn kn_negative_order_pos_zero_is_pole() {
    // K₋₃(+0) = K₃(+0): the m = |n| reduction (no sign) still poles.
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.kn(-3, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(!r.is_sign_negative());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn kn_neg_zero_is_invalid() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.kn(3, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn kn_negative_arg_is_invalid() {
    let a = BigFloat::try_from_i64_exact(-3, 53).expect("precision >= 1");
    let (r, status) = a.kn(2, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn kn_pos_inf_is_pos_zero() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.kn(5, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
    assert!(!status.invalid());
}

#[kani::proof]
fn kn_negative_order_pos_inf_is_pos_zero() {
    // K₋₅(+∞) = K₅(+∞) = +0 (order parity, no sign).
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.kn(-5, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
    assert!(!status.invalid());
}

#[kani::proof]
fn kn_neg_inf_is_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.kn(2, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}
