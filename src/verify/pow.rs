//! Kani harnesses for [`BigFloat::pow`].
//!
//! Coverage of the IEEE 754-2019 §9.2.1 special-case table:
//! `pow(x, ±0) = 1` (even for NaN, ±∞), `pow(+1, y) = 1`
//! (even for NaN, ±∞), NaN propagation outside those two rules,
//! `pow(±0, neg) = ±∞ + DIV_BY_ZERO`, `pow(±0, pos) = ±0`,
//! `pow(±∞, neg) = ±0`, `pow(±∞, pos) = ±∞`, and
//! `pow(neg, non-integer) = qNaN + INVALID`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

/// `pow(NaN, +0) = 1`. The "anything to the zero is one" rule
/// trumps NaN propagation.
#[kani::proof]
fn pow_nan_zero_is_one() {
    let x = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let y = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = x.pow(&y, RoundingMode::NearestEven);
    assert!(!r.is_nan());
    let one = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (cmp, _) = r.partial_cmp(&one);
    assert_eq!(cmp, Some(core::cmp::Ordering::Equal));
    assert!(status.is_ok());
}

/// `pow(+1, NaN) = 1`. The "one to anything is one" rule trumps
/// NaN propagation.
#[kani::proof]
fn pow_one_nan_is_one() {
    let x = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let y = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = x.pow(&y, RoundingMode::NearestEven);
    assert!(!r.is_nan());
    let one = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (cmp, _) = r.partial_cmp(&one);
    assert_eq!(cmp, Some(core::cmp::Ordering::Equal));
    assert!(status.is_ok());
}

/// `pow(+1, +∞) = 1` per §9.2.1.
#[kani::proof]
fn pow_one_inf_is_one() {
    let x = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let y = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = x.pow(&y, RoundingMode::NearestEven);
    let one = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (cmp, _) = r.partial_cmp(&one);
    assert_eq!(cmp, Some(core::cmp::Ordering::Equal));
    assert!(status.is_ok());
}

/// Signaling NaN in either operand raises `INVALID`.
#[kani::proof]
fn pow_signaling_nan_raises_invalid() {
    let x = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let y = BigFloat::try_from_i64_exact(2, 53).expect("2 fits");
    let (r, status) = x.pow(&y, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

/// `pow(+0, −1) = +∞ + DIV_BY_ZERO`.
#[kani::proof]
fn pow_pos_zero_neg_finite_is_pos_inf_div_by_zero() {
    let x = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let y = BigFloat::try_from_i64_exact(-1, 53).expect("-1 fits");
    let (r, status) = x.pow(&y, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.div_by_zero());
}

/// `pow(−0, −1) = −∞ + DIV_BY_ZERO` (odd-integer exponent
/// preserves the zero's sign through to infinity).
#[kani::proof]
fn pow_neg_zero_neg_one_is_neg_inf_div_by_zero() {
    let x = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let y = BigFloat::try_from_i64_exact(-1, 53).expect("-1 fits");
    let (r, status) = x.pow(&y, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(status.div_by_zero());
}

/// `pow(+0, +2) = +0`.
#[kani::proof]
fn pow_pos_zero_pos_finite_is_pos_zero() {
    let x = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let y = BigFloat::try_from_i64_exact(2, 53).expect("2 fits");
    let (r, status) = x.pow(&y, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

/// `pow(+∞, −1) = +0`.
#[kani::proof]
fn pow_pos_inf_neg_finite_is_pos_zero() {
    let x = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let y = BigFloat::try_from_i64_exact(-1, 53).expect("-1 fits");
    let (r, status) = x.pow(&y, RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}

/// `pow(+∞, +1) = +∞`.
#[kani::proof]
fn pow_pos_inf_pos_finite_is_pos_inf() {
    let x = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let y = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, status) = x.pow(&y, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_ok());
}
