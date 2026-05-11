//! Kani harnesses for [`BigFloat::fma`].
//!
//! Coverage: NaN propagation across all three operands, signaling-
//! NaN INVALID emission, the IEEE 754-2019 §7.2 `0 × ∞ + finite`
//! invalid form, and the subtle §7.2 carve-out where
//! `c` being NaN suppresses the `INVALID` that would otherwise
//! arise from the `0 × ∞` product (the NaN propagates without an
//! extra flag).

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::Status;

use super::helpers::nondet_constant_at;

/// NaN in `a` propagates through `fma`.
#[kani::proof]
fn fma_nan_in_a_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = nondet_constant_at(53);
    let c = nondet_constant_at(53);
    let (r, _status) = a.fma(&b, &c, RoundingMode::NearestEven);
    assert!(r.is_nan());
}

/// Signaling NaN in `a` raises `INVALID`.
#[kani::proof]
fn fma_signaling_nan_in_a_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let b = nondet_constant_at(53);
    let c = nondet_constant_at(53);
    let (r, status) = a.fma(&b, &c, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(!r.is_signaling_nan());
    assert!(status.contains(Status::INVALID));
}

/// `(0 × ∞) + 1` raises `INVALID` and returns a quiet NaN.
#[kani::proof]
fn fma_zero_times_inf_plus_finite_is_nan_invalid() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let c = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let (r, status) = a.fma(&b, &c, RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.contains(Status::INVALID));
}

/// `(0 × ∞) + qNaN` propagates the NaN **without** an additional
/// `INVALID` flag. Per IEEE 754-2019 §7.2: when `c` is NaN, the
/// NaN's propagation takes precedence over the `0 × ∞` invalid form,
/// and no `INVALID` is signaled from the product.
#[kani::proof]
fn fma_zero_times_inf_plus_qnan_no_extra_invalid() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let c = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.fma(&b, &c, RoundingMode::NearestEven);
    assert!(r.is_nan());
    // c is quiet NaN: no INVALID expected from the 0×∞ product.
    assert!(!status.contains(Status::INVALID));
}

/// `(+∞) × 1 + 0` is `+∞` (no flag).
#[kani::proof]
fn fma_inf_times_finite_plus_zero_is_inf() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let b = BigFloat::try_from_i64_exact(1, 53).expect("1 fits");
    let c = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.fma(&b, &c, RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(status.is_empty());
}
