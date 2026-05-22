//! Kani harnesses for [`BigFloat::beta`].
//!
//! `beta(a, b)` is total over the reals (ADR-0030, derived from
//! DLMF 5.12.1 / 5.5.3 / 5.2). These harnesses pin one input per
//! ADR-0030 domain row plus the NaN / signaling-NaN special cases:
//!
//! - **NaN / sNaN** propagation and `INVALID` (unchanged from
//!   slice 4c).
//! - **Row 2** — `a` negative non-integer, finite signed result;
//!   the combined sign is the product of three `gamma_sign_of`
//!   values, so both signs are exercised.
//! - **Row 5** — `a + b` a non-positive integer with `a, b` off the
//!   Γ poles: `+0` (denominator pole), no exception.
//! - **Row 4** — pole/pole cancellation to a finite value via the
//!   `(−1)^m / (m·C(n,m))` closed form.
//! - **Row 3** — a negative integer with no compensating `a + b`
//!   pole: `qNaN + INVALID` (two-sided sign-ambiguous pole).
//! - **Row 0** — a `±0` operand with no cancellation:
//!   `±∞ + DIV_BY_ZERO` (mirrors `gamma(±0)`).
//! - **Row 6** — both operands at Γ poles: `qNaN + INVALID`.
//!
//! The case-2 harnesses drive the `lgamma`/`exp` composition, so
//! they enter Kani's deep-unwind rounding loop (the documented
//! ADR-0012 continue-on-error advisory cost); the classification and
//! closed-form rows are cheap.
//!
//! Note: slice 4c shipped two harnesses asserting the old
//! positive-only reject (`beta(+0, 1)` and `beta(−1, 1)` →
//! `qNaN + INVALID`). ADR-0030 makes `beta(+0, 1) = +∞ +
//! DIV_BY_ZERO` (row 0) and `beta(−1, 1) = −1` (row 4 cancellation),
//! so those two were corrected here to the shipped behavior — the
//! same correction the 8a.2 unit tests received, not a new policy.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

const P: u32 = 53;
const NE: RoundingMode = RoundingMode::NearestEven;

/// `−1/2` as an exact dyadic (negative non-integer, off every Γ
/// pole).
fn neg_half() -> BigFloat {
    let neg_one = BigFloat::try_from_i64_exact(-1, P).expect("-1 fits");
    let two = BigFloat::try_from_i64_exact(2, P).expect("2 fits");
    neg_one.div(&two, NE).0
}

/// `k/4` as an exact dyadic.
fn quarter(k: i64) -> BigFloat {
    let num = BigFloat::try_from_i64_exact(k, P).expect("k fits");
    let four = BigFloat::try_from_i64_exact(4, P).expect("4 fits");
    num.div(&four, NE).0
}

#[kani::proof]
fn beta_nan_in_a_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, P, &[]).expect("precision >= 1");
    let b = BigFloat::try_from_i64_exact(1, P).expect("1 fits");
    let (r, _status) = a.beta(&b, NE);
    assert!(r.is_nan());
}

#[kani::proof]
fn beta_nan_in_b_propagates() {
    let a = BigFloat::try_from_i64_exact(1, P).expect("1 fits");
    let b = BigFloat::try_new_quiet_nan(Sign::Positive, P, &[]).expect("precision >= 1");
    let (r, _status) = a.beta(&b, NE);
    assert!(r.is_nan());
}

#[kani::proof]
fn beta_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, P, &[]).expect("precision >= 1");
    let b = BigFloat::try_from_i64_exact(1, P).expect("1 fits");
    let (r, status) = a.beta(&b, NE);
    assert!(r.is_nan());
    assert!(status.invalid());
}

/// ADR-0030 row 2: `B(−1/2, 1/4) > 0`, finite. Sign is
/// `sign Γ(−½)·sign Γ(¼)·sign Γ(−¼) = (−)(+)(−) = +`.
#[kani::proof]
fn beta_negative_non_integer_positive_sign() {
    let (r, status) = neg_half().beta(&quarter(1), NE);
    assert!(!r.is_nan() && !r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(!status.invalid() && !status.div_by_zero());
}

/// ADR-0030 row 2: `B(−1/2, 3/4) < 0`, finite. Sign is
/// `sign Γ(−½)·sign Γ(¾)·sign Γ(¼) = (−)(+)(+) = −`.
#[kani::proof]
fn beta_negative_non_integer_negative_sign() {
    let (r, status) = neg_half().beta(&quarter(3), NE);
    assert!(!r.is_nan() && !r.is_infinite());
    assert!(r.is_sign_negative());
    assert!(!status.invalid() && !status.div_by_zero());
}

/// ADR-0030 row 5: `a + b ∈ {0, −1, …}` with `a, b` off the Γ
/// poles. `B(−1/2, 1/2) = +0` exactly (denominator pole), no
/// exception.
#[kani::proof]
fn beta_sum_nonpos_integer_is_zero() {
    let (r, status) = neg_half().beta(&quarter(2), NE); // 1/4·2 = 1/2
    assert!(r.is_zero() && r.is_sign_positive());
    assert!(!status.invalid() && !status.div_by_zero());
}

/// ADR-0030 row 4: pole/pole cancellation. `a = −3` (negative
/// integer), `b = 2` (positive integer), `a + b = −1` a non-positive
/// integer, so `B(−3, 2) = (−1)^2 / (2·C(3,2)) = 1/6`, finite, no
/// exception. Closed-form path (no `lgamma`).
#[kani::proof]
fn beta_pole_cancellation_finite() {
    let a = BigFloat::try_from_i64_exact(-3, P).expect("-3 fits");
    let b = BigFloat::try_from_i64_exact(2, P).expect("2 fits");
    let (r, status) = a.beta(&b, NE);
    assert!(!r.is_nan() && !r.is_infinite());
    assert!(r.is_sign_positive());
    assert!(!status.invalid() && !status.div_by_zero());
}

/// ADR-0030 row 3: `a = −1` a negative integer, `b = 2`, `a + b = 1`
/// positive so no pole cancellation. Two-sided sign-ambiguous pole
/// → `qNaN + INVALID` (corrects the slice-4c `beta(−1, 1)` harness,
/// which is now row 4 cancellation).
#[kani::proof]
fn beta_negative_integer_uncancelled_invalid() {
    let a = BigFloat::try_from_i64_exact(-1, P).expect("-1 fits");
    let b = BigFloat::try_from_i64_exact(2, P).expect("2 fits");
    let (r, status) = a.beta(&b, NE);
    assert!(r.is_nan());
    assert!(status.invalid());
}

/// ADR-0030 row 0: `a = +0`, `b = 2`, no cancellation. Mirrors
/// `gamma(+0)`: `+∞ + DIV_BY_ZERO` (corrects the slice-4c
/// `beta(+0, 1)` harness, which asserted `qNaN + INVALID`).
#[kani::proof]
fn beta_zero_operand_signed_divzero() {
    let a = BigFloat::try_new_zero(Sign::Positive, P).expect("precision >= 1");
    let b = BigFloat::try_from_i64_exact(2, P).expect("2 fits");
    let (r, status) = a.beta(&b, NE);
    assert!(r.is_infinite() && r.is_sign_positive());
    assert!(status.div_by_zero());
}

/// ADR-0030 row 6: both operands at Γ poles. `B(−2, −3)` is a net
/// pole → `qNaN + INVALID`.
#[kani::proof]
fn beta_both_nonpos_integers_invalid() {
    let a = BigFloat::try_from_i64_exact(-2, P).expect("-2 fits");
    let b = BigFloat::try_from_i64_exact(-3, P).expect("-3 fits");
    let (r, status) = a.beta(&b, NE);
    assert!(r.is_nan());
    assert!(status.invalid());
}
