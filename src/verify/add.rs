//! Kani harnesses for [`BigFloat::add`].
//!
//! Slice 6a ships the canonical NaN / infinity / signed-zero
//! harnesses at fixed precision `53`. The proofs target the
//! special-case dispatch in [`crate::ops::addsub`]; finite-finite
//! arithmetic correctness is delegated to the proptest harness in
//! `tests/property_addsub.rs` and the MPFR differential lane in
//! `tests/differential_add.rs`.
//!
//! ## Strategy
//!
//! `BigFloat::add` routes through `add_with_signs`, which dispatches
//! NaN propagation, ±∞ ± ±∞, ±0 ± ±0, and zero + finite before
//! falling through to the alignment + rounding pipeline. The
//! harnesses bound their operands to the eight-constant set in
//! [`super::helpers::nondet_constant_at`], which keeps Kani from
//! exploring the finite-finite pipeline. The bounded-normal
//! generators that exercise the alignment pipeline at small fixed
//! precision land in slice 6b.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

use super::helpers::{nondet_constant_at, nondet_rounding_mode};

/// `∞ + (−∞)` raises `INVALID` and returns a quiet NaN.
#[kani::proof]
fn add_inf_minus_inf_is_nan_with_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.add(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(!r.is_signaling_nan());
    assert!(status.invalid());
}

/// Same-sign infinities produce that signed infinity with no flag.
#[kani::proof]
fn add_inf_plus_inf_same_sign_is_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.add(&b, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

/// NaN propagates through `add` for any second operand drawn from
/// the canonical set.
#[kani::proof]
fn add_nan_propagates_under_canonical_operands() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = nondet_constant_at(53);
    let (r, _status) = a.add(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

/// Signaling NaN in either operand raises `INVALID`, with the result
/// quieted per IEEE 754-2019 §6.2.
#[kani::proof]
fn add_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = nondet_constant_at(53);
    let (r, status) = a.add(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(!r.is_signaling_nan());
    assert!(status.invalid());
}

/// `(+0) + (−0)` is `+0` under every mode except `TowardNegative`,
/// where it is `−0`. Per IEEE 754-2019 §6.3.
#[kani::proof]
fn add_pos_zero_plus_neg_zero_sign_rule() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let mode = nondet_rounding_mode();
    let (r, status) = a.add(&b, mode);
    assert!(r.is_zero());
    assert!(status.is_ok());
    if matches!(mode, RoundingMode::TowardNegative) {
        assert!(r.is_sign_negative());
    } else {
        assert!(r.is_sign_positive());
    }
}

/// `(−0) + (−0)` is `−0` under every rounding mode.
#[kani::proof]
fn add_neg_zero_plus_neg_zero_is_neg_zero() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let b = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let mode = nondet_rounding_mode();
    let (r, status) = a.add(&b, mode);
    assert!(r.is_zero());
    assert!(r.is_sign_negative());
    assert!(status.is_ok());
}
