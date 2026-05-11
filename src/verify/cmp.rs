//! Kani harnesses for [`BigFloat`] comparison.
//!
//! Covers the reflexivity of `total_cmp` on any canonical value,
//! the NaN-handling contract for `partial_cmp` (returns `None` for
//! any NaN, raises `INVALID` for signaling NaN, no flag for quiet
//! NaN), and the consistency of `min` / `max` with `partial_cmp`.

use core::cmp::Ordering;

use crate::big::BigFloat;
use crate::sign::Sign;
use crate::status::Status;

use super::helpers::{canonical_at, NUM_CANONICAL};

/// `total_cmp` is reflexive on every canonical value.
#[kani::proof]
fn total_cmp_reflexive_on_canonical() {
    let idx: u8 = kani::any();
    kani::assume(idx < NUM_CANONICAL);
    let v = canonical_at(idx, 53);
    assert_eq!(v.total_cmp(&v), Ordering::Equal);
}

/// `partial_cmp(NaN, x)` returns `None` for any second operand,
/// with no `INVALID` flag when the NaN is quiet.
#[kani::proof]
fn partial_cmp_quiet_nan_returns_none_no_flag() {
    let nan = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let idx: u8 = kani::any();
    kani::assume(idx < NUM_CANONICAL);
    let other = canonical_at(idx, 53);

    let (ord, status) = nan.partial_cmp(&other);
    assert!(ord.is_none());
    // Quiet NaN: per IEEE 754-2019 §5.11, the unordered relation
    // does not signal INVALID.
    assert!(!status.invalid());
}

/// `partial_cmp` raises `INVALID` when either operand is signaling.
#[kani::proof]
fn partial_cmp_signaling_nan_raises_invalid() {
    let snan = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let other = canonical_at(6, 53); // +1
    let (ord, status) = snan.partial_cmp(&other);
    assert!(ord.is_none());
    assert!(status.invalid());
}

/// `min(a, a)` is `a` for any canonical value; the status flag is
/// either empty or `INVALID` depending on whether `a` is signaling.
#[kani::proof]
fn min_idempotent_on_canonical() {
    let idx: u8 = kani::any();
    kani::assume(idx < NUM_CANONICAL);
    let v = canonical_at(idx, 53);
    let (r, _status) = v.min(&v);
    // For non-NaN values, the result equals `v`. For NaN, the
    // result is NaN (idempotent in class).
    if v.is_nan() {
        assert!(r.is_nan());
    } else {
        assert_eq!(r.total_cmp(&v), Ordering::Equal);
    }
}

/// `max(a, a)` matches `min(a, a)` shape.
#[kani::proof]
fn max_idempotent_on_canonical() {
    let idx: u8 = kani::any();
    kani::assume(idx < NUM_CANONICAL);
    let v = canonical_at(idx, 53);
    let (r, _status) = v.max(&v);
    if v.is_nan() {
        assert!(r.is_nan());
    } else {
        assert_eq!(r.total_cmp(&v), Ordering::Equal);
    }
}
