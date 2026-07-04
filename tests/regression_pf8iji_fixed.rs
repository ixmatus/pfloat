//! Regression tests for epic pf-8iji review remediation, defect
//! pf-sn3n (ADR-0119).
//!
//! Before the fix, `FixedFloat<0>` was instantiable: a const-generic
//! fixed-precision float with zero significand bits. That is a
//! nonsensical illegal state. Its methods panicked while claiming a
//! `PREC >= 1` bound that did not exist, and `to_big()` leaked a
//! precision-0 `BigFloat`, violating `BigFloat`'s own `precision >= 1`
//! invariant.
//!
//! The remedy makes `FixedFloat<0>` uninstantiable at compile time.
//! Each value-birthing constructor carries a method-level where-clause
//! bound (the internal `require_nonzero_precision` guard) that fails
//! const-evaluation for `PREC == 0`. The type name stays inert: it may
//! appear in a signature (`Option<FixedFloat<0>>` is a well-formed
//! type), but no value of it can be born in safe code, so the illegal
//! state is gone.
//!
//! # The exclusion is a compile-time property
//!
//! The guarantee lives in the type system, so it holds without a
//! runtime assertion; a test can only pin the positive side. Every
//! constructor call below is a COMPILE ERROR and must stay one —
//! verified manually with `rustc` in a consumer that enables
//! `generic_const_exprs`. Each fails const-evaluation of the
//! `require_nonzero_precision(PREC)` bound, e.g.:
//!
//! ```text
//! error[E0080]: evaluation panicked: FixedFloat<0> is not a valid
//! type: a fixed-precision float needs at least one significand bit
//! (PREC >= 1)
//! ```
//!
//! for each of:
//!
//! ```text
//! use pfloat::FixedFloat;
//! let _ = FixedFloat::<0>::zero();                 // does not compile
//! let _ = FixedFloat::<0>::infinity();             // does not compile
//! let _ = FixedFloat::<0>::nan(Sign::Positive);    // does not compile
//! let _ = FixedFloat::<0>::try_from_i64_round(0, RoundingMode::NearestEven);
//! ```
//!
//! (A machine-checked `compile_fail` doctest is deliberately avoided:
//! it would pass spuriously in any build that omits the `fixed`
//! feature, where `FixedFloat` is simply absent, or that forgets the
//! `generic_const_exprs` consumer feature, where every use is an
//! unrelated E0275. The guarantee is structural; these tests pin that
//! every legal precision still works and that `to_big()` never yields
//! a precision-0 `BigFloat`.)

// A consumer crate that names `FixedFloat<PREC>` must enable the same
// nightly feature the library declares, so its `[(); limbs_for(PREC)]`
// bound can be evaluated (otherwise every use is an E0275 well-formed
// overflow, unrelated to this fix).
#![feature(generic_const_exprs)]
#![allow(incomplete_features)]

use core::cmp::Ordering;

use pfloat::{BigFloat, FixedFloat, RoundingMode, Sign};

/// `PREC == 1` is the smallest legal precision (one significand bit).
/// It must construct, classify, and round-trip through `BigFloat`
/// with a precision that honors the `precision >= 1` invariant.
#[test]
fn smallest_legal_precision_one_works() {
    // The integer 1 has exactly one significant bit, so it is exact
    // at one-bit precision.
    let one = FixedFloat::<1>::try_from_i64_exact(1).unwrap();
    assert!(one.is_normal());
    assert!(one.is_sign_positive());
    assert_eq!(one.precision(), 1);
    assert_eq!(FixedFloat::<1>::PRECISION, 1);
    assert_eq!(FixedFloat::<1>::LIMBS, 1);

    // to_big() cannot yield a precision-0 BigFloat.
    let big = one.to_big();
    assert_eq!(big.precision(), 1);
    assert!(big.precision() >= 1);
}

/// Special-value constructors and sign handling at the floor
/// precision, and their `BigFloat` projection.
#[test]
fn precision_one_specials_project_soundly() {
    let z = FixedFloat::<1>::zero();
    assert!(z.is_zero() && z.is_sign_positive());
    assert!(z.to_big().precision() >= 1);

    let neg_inf = FixedFloat::<1>::neg_infinity();
    assert!(neg_inf.is_infinite() && neg_inf.is_sign_negative());
    assert!(neg_inf.to_big().precision() >= 1);

    let nan = FixedFloat::<1>::nan(Sign::Positive);
    assert!(nan.is_quiet_nan());
    assert_eq!(nan.to_big().precision(), 1);
}

/// The common IEEE 754 significand widths are unaffected: binary32
/// (24), binary64 (53), binary128 (113) still construct and compute.
#[test]
fn common_ieee_widths_still_work() {
    let a = FixedFloat::<24>::try_from_i64_exact(3).unwrap();
    let b = FixedFloat::<24>::try_from_i64_exact(4).unwrap();
    let (sum, status) = a.add(&b, RoundingMode::NearestEven);
    assert!(status.is_ok());
    let seven = FixedFloat::<24>::try_from_i64_exact(7).unwrap();
    assert_eq!(sum.partial_cmp(&seven).0, Some(Ordering::Equal));
    assert_eq!(sum.precision(), 24);

    assert_eq!(FixedFloat::<53>::PRECISION, 53);
    assert_eq!(FixedFloat::<113>::PRECISION, 113);
    assert_eq!(FixedFloat::<256>::PRECISION, 256);
}

/// Across several legal precisions, the `BigFloat` projection always
/// preserves the fixed precision, which is `>= 1` by the type-level
/// bound. This is the exact leak the fix closes.
#[test]
fn to_big_never_yields_precision_zero() {
    let one_big: BigFloat = FixedFloat::<1>::try_from_i64_exact(1).unwrap().into();
    assert!(one_big.precision() >= 1);

    let two_big = FixedFloat::<2>::zero().to_big();
    assert!(two_big.precision() >= 1);

    let wide_big = FixedFloat::<113>::try_from_i64_exact(5).unwrap().to_big();
    assert_eq!(wide_big.precision(), 113);
    assert!(wide_big.precision() >= 1);
}
