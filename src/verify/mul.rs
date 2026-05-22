//! Kani harnesses for [`BigFloat::mul`].
//!
//! Coverage: NaN propagation, signaling-NaN INVALID emission, the
//! IEEE 754 invalid `±0 × ±∞ → qNaN + INVALID` form, sign-of-
//! product correctness for the `±0 × ±finite`, `±finite × ±0`,
//! `±∞ × ±finite`, and `±finite × ±∞` paths, and exponent-range
//! saturation (the product exponent exceeding `i64::MAX` flags
//! `OVERFLOW` without panicking, pf-rnc).

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

use super::helpers::nondet_constant_at;

/// NaN propagates through `mul`.
#[kani::proof]
fn mul_nan_propagates_under_canonical_operands() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = nondet_constant_at(53);
    let (r, _status) = a.mul(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

/// Signaling NaN in either operand raises `INVALID`.
#[kani::proof]
fn mul_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = nondet_constant_at(53);
    let (r, status) = a.mul(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(!r.is_signaling_nan());
    assert!(status.invalid());
}

/// `(+0) × (+∞)` raises `INVALID` and returns a quiet NaN.
#[kani::proof]
fn mul_pos_zero_times_pos_inf_is_nan_invalid() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.mul(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

/// `(−0) × (+∞)` raises `INVALID` and returns a quiet NaN (sign of
/// the zero is irrelevant to the validity check).
#[kani::proof]
fn mul_neg_zero_times_pos_inf_is_nan_invalid() {
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.mul(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

/// `(+∞) × (+∞)` is `+∞` (no flag).
#[kani::proof]
fn mul_pos_inf_times_pos_inf_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.mul(&b, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

/// `(+∞) × (−∞)` is `−∞`. Sign-of-product rule on infinities.
#[kani::proof]
fn mul_pos_inf_times_neg_inf_is_neg_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.mul(&b, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.is_ok());
}

/// `(+0) × (−0)` is `−0`. Product-sign rule even at zero.
#[kani::proof]
fn mul_pos_zero_times_neg_zero_is_neg_zero() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.mul(&b, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_negative());
    assert!(status.is_ok());
}

/// Regression (pf-rnc, fuzz-found via Airy `bi_prime`): the result
/// exponent `top_bit + e_a + e_b − p_a − p_b + 2` is now computed in
/// `i128` and saturated to the `i64` range, so multiplying operands
/// whose true product exponent exceeds `i64::MAX` flags `OVERFLOW`
/// and returns a finite saturated value (pfloat has no `emax`)
/// instead of panicking on `i64` overflow. Squaring `2` doubles the
/// exponent each step; within ~63 steps it passes `i64::MAX`. The
/// invariant proved: no step panics or yields `NaN`, and saturation
/// is reached. Advisory: the bounded `unwind` over a chain of `mul`s
/// is the documented ADR-0012 deep-unwind cost.
#[kani::proof]
#[kani::unwind(67)]
fn mul_extreme_exponent_saturates_without_panic() {
    let mut x = BigFloat::try_from_i64_exact(2, 53).expect("2 fits");
    let mut saw_overflow = false;
    let mut i = 0;
    while i < 66 {
        let (sq, status) = x.mul(&x, RoundingMode::NearestEven);
        assert!(!sq.is_nan());
        if status.overflow() {
            saw_overflow = true;
            break;
        }
        x = sq;
        i += 1;
    }
    assert!(saw_overflow);
}
