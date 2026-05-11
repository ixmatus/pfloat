//! Kani harnesses for [`BigFloat`] classification.
//!
//! Covers the totality invariant (every value is exactly one of
//! NaN, infinite, finite-nonzero, or zero), the sign-of-zero
//! preservation through `abs` / `neg` / `signum`, and the
//! signaling-NaN / quiet-NaN partition.

use crate::big::BigFloat;
use crate::sign::Sign;

use super::helpers::{canonical_at, NUM_CANONICAL};

/// Each canonical value satisfies exactly one of the four
/// classification predicates: NaN, infinite, zero, normal.
#[kani::proof]
fn classify_totality_on_canonical() {
    let idx: u8 = kani::any();
    kani::assume(idx < NUM_CANONICAL);
    let v = canonical_at(idx, 53);

    let count = u32::from(v.is_nan())
        + u32::from(v.is_infinite())
        + u32::from(v.is_zero())
        + u32::from(v.is_normal());
    assert_eq!(count, 1);
}

/// `is_finite ↔ !is_nan ∧ !is_infinite` on the canonical set.
#[kani::proof]
fn classify_finite_partition_on_canonical() {
    let idx: u8 = kani::any();
    kani::assume(idx < NUM_CANONICAL);
    let v = canonical_at(idx, 53);

    let finite = v.is_finite();
    let non_finite = v.is_nan() || v.is_infinite();
    assert_ne!(finite, non_finite);
}

/// Signaling NaN implies NaN; quiet NaN is also NaN.
#[kani::proof]
fn classify_signaling_implies_nan() {
    let snan = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    assert!(snan.is_signaling_nan());
    assert!(snan.is_nan());

    let qnan = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    assert!(!qnan.is_signaling_nan());
    assert!(qnan.is_nan());
}

/// `abs` clears the sign bit of every canonical non-NaN value.
/// (NaN abs is left unspecified for sign-bit interpretation; the
/// classification stays NaN.)
#[kani::proof]
fn abs_clears_sign_on_finite_canonical() {
    let idx: u8 = kani::any();
    kani::assume(idx < NUM_CANONICAL);
    let v = canonical_at(idx, 53);
    kani::assume(!v.is_nan());

    let r = v.abs();
    assert!(r.is_sign_positive());
}

/// `signum` returns +1, -1, +0, -0, or NaN, matching the sign of
/// the input.
#[kani::proof]
fn signum_preserves_sign_on_canonical() {
    let idx: u8 = kani::any();
    kani::assume(idx < NUM_CANONICAL);
    let v = canonical_at(idx, 53);
    kani::assume(!v.is_nan());

    let r = v.signum();
    if v.is_sign_positive() {
        assert!(r.is_sign_positive());
    } else {
        assert!(r.is_sign_negative());
    }
}
