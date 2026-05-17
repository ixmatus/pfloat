//! Kani harnesses for [`BigFloat::j0`] (Bessel `J0`, entire on ℝ).
//!
//! `J0(NaN) = NaN`; `J0(±0) = 1` (exact, DLMF 10.2.2); `J0(±∞) = +0`
//! (decaying-envelope convention, ADR-0023).

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn j0_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.j0(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn j0_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.j0(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn j0_pos_zero_is_one() {
    // J0(+0) = 1 exactly (DLMF 10.2.2).
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, _status) = a.j0(RoundingMode::NearestEven);
    let one = BigFloat::try_from_i64_exact(1, 53).expect("precision >= 1");
    assert!(r.partial_cmp(&one).0 == Some(core::cmp::Ordering::Equal));
}

#[kani::proof]
fn j0_neg_zero_is_one() {
    // Bessel J is entire: J0(−0) = J0(+0) = 1.
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, _status) = a.j0(RoundingMode::NearestEven);
    let one = BigFloat::try_from_i64_exact(1, 53).expect("precision >= 1");
    assert!(r.partial_cmp(&one).0 == Some(core::cmp::Ordering::Equal));
}

#[kani::proof]
fn j0_pos_inf_is_pos_zero() {
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, _status) = a.j0(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
}

#[kani::proof]
fn j0_neg_inf_is_pos_zero() {
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, _status) = a.j0(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
}
