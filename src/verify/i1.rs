//! Kani harnesses for [`BigFloat::i1`] (modified Bessel `I1`,
//! entire, odd in the argument).
//!
//! `I1(NaN) = NaN`; `I1(±0) = 0` (exact, DLMF 10.30.1); `I1(+∞) =
//! +∞`, `I1(−∞) = −∞` (odd order, the argument parity
//! `I1(−x) = −I1(x)`); a negative argument is finite and never
//! `INVALID` (`I` is entire). The Normal series is out of Kani
//! scope.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn i1_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.i1(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn i1_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.i1(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn i1_pos_zero_is_zero() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.i1(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!status.invalid());
}

#[kani::proof]
fn i1_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.i1(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(!r.is_sign_negative());
    assert!(!status.invalid());
}

#[kani::proof]
fn i1_neg_inf_is_neg_inf() {
    // Odd order: I1(−∞) = −∞ (argument parity I1(−x) = −I1(x)).
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.i1(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(!status.invalid());
}
