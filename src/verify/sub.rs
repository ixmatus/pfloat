//! Kani harnesses for [`BigFloat::sub`].
//!
//! Subtraction routes through `add_with_signs` with the second
//! operand's sign flipped. The harnesses cover NaN propagation,
//! signaling-NaN INVALID emission, ±∞ − ±∞ edge cases, and the
//! sign-of-zero rule for `(±0) − (±0)` (which mirrors `(±0) ± (∓0)`
//! per IEEE 754-2019 §6.3).

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

use super::helpers::{nondet_constant_at, nondet_rounding_mode};

/// `∞ − ∞` raises `INVALID` and returns a quiet NaN.
#[kani::proof]
fn sub_inf_minus_inf_is_nan_with_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.sub(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(!r.is_signaling_nan());
    assert!(status.invalid());
}

/// `∞ − (−∞)` is `+∞` (no flag).
#[kani::proof]
fn sub_inf_minus_neg_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.sub(&b, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

/// NaN propagates through `sub` for any second operand drawn from
/// the canonical set.
#[kani::proof]
fn sub_nan_propagates_under_canonical_operands() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = nondet_constant_at(53);
    let (r, _status) = a.sub(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

/// Signaling NaN in either operand raises `INVALID`.
#[kani::proof]
fn sub_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = nondet_constant_at(53);
    let (r, status) = a.sub(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(!r.is_signaling_nan());
    assert!(status.invalid());
}

/// `(+0) − (+0)` is `+0` under every mode except `TowardNegative`,
/// where it is `−0`. Equivalent to `add(+0, −0)`.
#[kani::proof]
fn sub_pos_zero_minus_pos_zero_sign_rule() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let mode = nondet_rounding_mode();
    let (r, status) = a.sub(&b, mode);
    assert!(r.is_zero());
    assert!(status.is_ok());
    if matches!(mode, RoundingMode::TowardNegative) {
        assert!(r.is_sign_negative());
    } else {
        assert!(r.is_sign_positive());
    }
}

/// `(+0) − (−0)` is `+0` under every rounding mode.
#[kani::proof]
fn sub_pos_zero_minus_neg_zero_is_pos_zero() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let mode = nondet_rounding_mode();
    let (r, status) = a.sub(&b, mode);
    assert!(r.is_zero());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}
