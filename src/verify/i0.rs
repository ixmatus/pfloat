//! Kani harnesses for [`BigFloat::i0`] (modified Bessel `I0`,
//! entire on the real line).
//!
//! `I0(NaN) = NaN`; `I0(±0) = 1` (exact, DLMF 10.30.1); `I0(±∞) =
//! +∞` (a genuine infinite limit, `Status::OK`, not the
//! decaying-envelope convention); a negative argument is finite and
//! never `INVALID` (`I` is entire). The Normal series is out of Kani
//! scope; the special-value arms return before the evaluator.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn i0_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, _status) = a.i0(RoundingMode::NearestEven);
    assert!(r.is_nan());
}

#[kani::proof]
fn i0_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.i0(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn i0_pos_zero_is_one() {
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.i0(RoundingMode::NearestEven);
    let one = BigFloat::try_from_i64_exact(1, 53).expect("precision >= 1");
    assert!(matches!(
        r.partial_cmp(&one).0,
        Some(core::cmp::Ordering::Equal)
    ));
    assert!(!status.invalid());
}

#[kani::proof]
fn i0_neg_zero_is_one() {
    // I is entire: −0 behaves like +0.
    let a = BigFloat::try_new_zero(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.i0(RoundingMode::NearestEven);
    let one = BigFloat::try_from_i64_exact(1, 53).expect("precision >= 1");
    assert!(matches!(
        r.partial_cmp(&one).0,
        Some(core::cmp::Ordering::Equal)
    ));
    assert!(!status.invalid());
}

#[kani::proof]
fn i0_pos_inf_is_pos_inf() {
    // I0(+∞) = +∞, a genuine infinite limit, Status::OK.
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.i0(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(!r.is_sign_negative());
    assert!(!status.invalid());
}

#[kani::proof]
fn i0_neg_inf_is_pos_inf() {
    // I0 is even (order 0): I0(−∞) = +∞.
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.i0(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(!r.is_sign_negative());
    assert!(!status.invalid());
}
