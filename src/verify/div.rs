//! Kani harnesses for [`BigFloat::div`].
//!
//! Coverage: NaN propagation, signaling-NaN INVALID, the invalid
//! `0/0` and `∞/∞` forms, the `DIV_BY_ZERO` flag for
//! `finite_nonzero / 0`, sign-of-quotient for the
//! `0 / finite_nonzero` and `∞ / finite` cases, and exponent-range
//! saturation (a quotient exponent past `i64::MAX` flags `OVERFLOW`
//! without panicking, pf-rnc).

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

use super::helpers::nondet_constant_at;

/// NaN propagates through `div`.
#[kani::proof]
fn div_nan_propagates_under_canonical_operands() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = nondet_constant_at(53);
    let (r, _status) = a.div(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

/// Signaling NaN in either operand raises `INVALID`.
#[kani::proof]
fn div_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = nondet_constant_at(53);
    let (r, status) = a.div(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(!r.is_signaling_nan());
    assert!(status.invalid());
}

/// `0 / 0` raises `INVALID` and returns a quiet NaN.
#[kani::proof]
fn div_zero_by_zero_is_nan_invalid() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.div(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

/// `∞ / ∞` raises `INVALID` and returns a quiet NaN.
#[kani::proof]
fn div_inf_by_inf_is_nan_invalid() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.div(&b, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

/// `1 / (+0)` is `+∞` and raises `DIV_BY_ZERO` per IEEE 754-2019 §7.3.
#[kani::proof]
fn div_pos_finite_by_pos_zero_is_pos_inf_div_by_zero() {
    let a = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let b = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.div(&b, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.div_by_zero());
}

/// `1 / (−0)` is `−∞` and raises `DIV_BY_ZERO`.
#[kani::proof]
fn div_pos_finite_by_neg_zero_is_neg_inf_div_by_zero() {
    let a = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let b = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.div(&b, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.div_by_zero());
}

/// `0 / 1` is `+0` (no flag); sign comes from the dividend.
#[kani::proof]
fn div_pos_zero_by_pos_finite_is_pos_zero() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, status) = a.div(&b, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

/// `(+∞) / 1` is `+∞` (no flag).
#[kani::proof]
fn div_pos_inf_by_pos_finite_is_pos_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, status) = a.div(&b, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

/// Regression (pf-rnc, fuzz-found via Airy `bi_prime`): the quotient
/// exponent is now computed in `i128` and saturated to the `i64`
/// range, so a quotient whose true exponent exceeds `i64::MAX` flags
/// `OVERFLOW` and returns a finite saturated value instead of
/// panicking on `i64` overflow. Square `2` until the next square
/// would saturate, take the reciprocal, then `big / tiny` has
/// exponent past `i64::MAX`. No step panics or yields `NaN`, and
/// saturation is reached. Advisory: bounded `unwind` over a chain
/// of `mul`s is the documented ADR-0012 deep-unwind cost.
#[kani::proof]
#[kani::unwind(67)]
fn div_extreme_exponent_saturates_without_panic() {
    let mut big = BigFloat::try_from_i64_exact(2, 53).expect("2 fits");
    let mut i = 0;
    while i < 66 {
        let (sq, st) = big.mul(&big, RoundingMode::NearestEven);
        assert!(!sq.is_nan());
        if st.overflow() {
            break;
        }
        big = sq;
        i += 1;
    }
    let one = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (tiny, _) = one.div(&big, RoundingMode::NearestEven);
    let (q, st) = big.div(&tiny, RoundingMode::NearestEven);
    assert!(!q.is_nan());
    assert!(st.overflow());
}
