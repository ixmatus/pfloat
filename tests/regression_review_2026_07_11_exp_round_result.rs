//! Regression guard for pf-291u (review finding `elementary/type/4`),
//! adjudicated CONFIRMED and fixed under ADR-0129.
//!
//! `exp_round` was the only elementary `*_round` kernel that returned
//! `(Self, Status)` and panicked in release on `target_precision == 0`
//! (the `debug_assert` was compiled out, then `exp_kernel` hit an
//! `.expect`). It now returns `Result<(Self, Status), BuildError>` and
//! yields `Err(BuildError::PrecisionZero)` for target 0, exactly like
//! `ln_round`/`sin_round`/etc.
//!
//! Run: `cargo test --release --features std,big,exp-log
//! --test regression_review_2026_07_11_exp_round_result`.

#![cfg(all(feature = "big", feature = "exp-log"))]

use pfloat::{BigFloat, BuildError, RoundingMode};

const NE: RoundingMode = RoundingMode::NearestEven;

#[test]
fn exp_round_precision_zero_is_err_not_panic() {
    // The defect: this panicked in release. It must now be an `Err`.
    let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
    assert!(matches!(
        two.exp_round(0, NE),
        Err(BuildError::PrecisionZero)
    ));
}

#[test]
fn exp_round_nonzero_precision_is_ok() {
    let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
    let (v, s) = two.exp_round(53, NE).expect("precision 53 is valid");
    assert_eq!(v.precision(), 53);
    assert!(s.inexact(), "e^2 is irrational");
}

#[test]
fn exp_convenience_still_infallible() {
    // `exp` unwraps internally on the always-nonzero `self.precision`;
    // it keeps its `(Self, Status)` shape and never panics.
    let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
    let (v, _s) = two.exp(NE);
    assert_eq!(v.precision(), 53);
}
