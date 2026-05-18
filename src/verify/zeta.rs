//! Kani harnesses for [`BigFloat::zeta`] (Riemann zeta, real
//! argument; slice 6r).
//!
//! The harnesses cover the special-value dispatch, which returns
//! before the heavy Borwein / functional-equation evaluator:
//! `ζ(NaN) = NaN`; `sNaN` raises `INVALID`; the simple pole
//! `ζ(1) = +∞` raising `DIV_BY_ZERO` (residue +1, `+∞` the
//! `s → 1⁺` side, DLMF 25.2); `ζ(0) = −1/2` exact (DLMF 25.6.1);
//! the trivial zero `ζ(−2) = +0` (DLMF 25.6.4); `ζ(+∞) = 1` (the
//! genuine Dirichlet limit); `ζ(−∞) = NaN` + `INVALID` (the
//! unbounded-non-converging convention, the K-vs-Y `+∞`
//! distinction, ADR-0025/0026). The Normal accelerator/FE path is
//! out of Kani scope (the [`super::k0`] precedent: the iterative
//! kernel is covered by the differential and property suites).

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[kani::proof]
fn zeta_nan_propagates() {
    let a = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.zeta(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(!status.invalid());
}

#[kani::proof]
fn zeta_signaling_nan_raises_invalid() {
    let a = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).expect("precision >= 1");
    let (r, status) = a.zeta(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}

#[kani::proof]
fn zeta_pole_at_one_is_div_by_zero() {
    // ζ(1) = +∞ + DIV_BY_ZERO (simple pole, residue +1; +∞ the
    // s → 1⁺ side, the Ci/li/K pole convention).
    let a = BigFloat::try_from_i64_exact(1, 53).expect("precision >= 1");
    let (r, status) = a.zeta(RoundingMode::NearestEven);
    assert!(r.is_infinite());
    assert!(!r.is_sign_negative());
    assert!(status.div_by_zero());
}

#[kani::proof]
fn zeta_at_zero_is_minus_half() {
    // ζ(0) = −1/2 exact: finite, negative, non-zero, OK.
    let a = BigFloat::try_new_zero(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.zeta(RoundingMode::NearestEven);
    assert!(!r.is_nan());
    assert!(!r.is_infinite());
    assert!(!r.is_zero());
    assert!(r.is_sign_negative());
    assert!(!status.invalid());
    assert!(!status.div_by_zero());
}

#[kani::proof]
fn zeta_trivial_zero_at_minus_two() {
    // ζ(−2) = +0 exact (a trivial zero), Status::OK.
    let a = BigFloat::try_from_i64_exact(-2, 53).expect("precision >= 1");
    let (r, status) = a.zeta(RoundingMode::NearestEven);
    assert!(r.is_zero());
    assert!(!r.is_sign_negative());
    assert!(!status.invalid());
}

#[kani::proof]
fn zeta_pos_inf_is_one() {
    // ζ(+∞) = 1, the genuine Dirichlet limit: finite, positive,
    // non-zero, OK.
    let a = BigFloat::try_new_infinity(Sign::Positive, 53).expect("precision >= 1");
    let (r, status) = a.zeta(RoundingMode::NearestEven);
    assert!(!r.is_nan());
    assert!(!r.is_infinite());
    assert!(!r.is_zero());
    assert!(!r.is_sign_negative());
    assert!(!status.invalid());
}

#[kani::proof]
fn zeta_neg_inf_is_invalid() {
    // ζ(−∞) = NaN + INVALID (unbounded non-converging oscillation;
    // not the decaying-envelope convention).
    let a = BigFloat::try_new_infinity(Sign::Negative, 53).expect("precision >= 1");
    let (r, status) = a.zeta(RoundingMode::NearestEven);
    assert!(r.is_nan());
    assert!(status.invalid());
}
