//! Correctly-rounded `f64` elementary and reciprocal functions.
//!
//! The `f64` companion to [`crate::f32`]; the same widen, evaluate, and
//! enclose-then-round contract holds at the wider format. The `f64`
//! surface rests on differential testing plus worst-case vectors rather
//! than an exhaustive sweep, since the 2^64 input space cannot be
//! enumerated.

use pfloat::BigFloat;

use crate::round::{drive, unary, F64Shell};
use crate::{RoundingMode, Status};

// --- exponentials (saturation fast-path: see `crate::saturate`) ---
unary!(
    f64,
    F64Shell,
    exp,
    exp_round,
    direct_sat,
    crate::saturate::sat_exp,
    "`exp`"
);
unary!(
    f64,
    F64Shell,
    exp2,
    exp2_round,
    result_sat,
    crate::saturate::sat_exp2,
    "`exp2`"
);
unary!(
    f64,
    F64Shell,
    exp10,
    exp10_round,
    result_sat,
    crate::saturate::sat_exp10,
    "`exp10`"
);
unary!(
    f64,
    F64Shell,
    expm1,
    expm1_round,
    result_sat,
    crate::saturate::sat_expm1,
    "`expm1`"
);

// --- logarithms ---
unary!(f64, F64Shell, ln, ln_round, result, "`ln`");
unary!(f64, F64Shell, log2, log2_round, result, "`log2`");
unary!(f64, F64Shell, log10, log10_round, result, "`log10`");
unary!(f64, F64Shell, log1p, log1p_round, result, "`log1p`");

// --- roots ---
unary!(f64, F64Shell, sqrt, sqrt_round, result, "`sqrt`");
unary!(f64, F64Shell, cbrt, cbrt_round, result, "`cbrt`");

// --- circular ---
unary!(f64, F64Shell, sin, sin_round, result, "`sin`");
unary!(f64, F64Shell, cos, cos_round, result, "`cos`");
unary!(f64, F64Shell, tan, tan_round, result, "`tan`");
unary!(f64, F64Shell, cot, cot_round, result, "`cot`");
unary!(f64, F64Shell, sec, sec_round, result, "`sec`");
unary!(f64, F64Shell, csc, csc_round, result, "`csc`");

// --- inverse circular ---
unary!(f64, F64Shell, asin, asin_round, result, "`asin`");
unary!(f64, F64Shell, acos, acos_round, result, "`acos`");
unary!(f64, F64Shell, atan, atan_round, result, "`atan`");

// --- hyperbolic (saturation fast-path: see `crate::saturate`) ---
unary!(
    f64,
    F64Shell,
    sinh,
    sinh_round,
    result_sat,
    crate::saturate::sat_sinh,
    "`sinh`"
);
unary!(
    f64,
    F64Shell,
    cosh,
    cosh_round,
    result_sat,
    crate::saturate::sat_cosh,
    "`cosh`"
);
unary!(
    f64,
    F64Shell,
    tanh,
    tanh_round,
    result_sat,
    crate::saturate::sat_tanh,
    "`tanh`"
);

// --- inverse hyperbolic ---
unary!(f64, F64Shell, asinh, asinh_round, result, "`asinh`");
unary!(f64, F64Shell, acosh, acosh_round, result, "`acosh`");
unary!(f64, F64Shell, atanh, atanh_round, result, "`atanh`");

/// `hypot(x, y) = sqrt(x² + y²)`, correctly rounded to `f64` under
/// round-to-nearest-even.
#[must_use]
pub fn hypot(x: f64, y: f64) -> f64 {
    hypot_round(x, y, RoundingMode::NearestEven).0
}

/// `hypot(x, y)` correctly rounded to `f64` under `mode`, with the IEEE
/// 754 status flags raised. IEEE 754-2019 §9.2.1 special cases are
/// inherited from pfloat's kernel.
#[must_use]
pub fn hypot_round(x: f64, y: f64, mode: RoundingMode) -> (f64, Status) {
    let yb = BigFloat::from_f64(y);
    drive::<F64Shell>(x, mode, move |xb, w, dir| {
        xb.hypot_round(&yb, w, dir)
            .expect("w = PREC + guard >= 1: BuildError only on precision 0")
    })
}

/// `rootn(x, n)`, the real `n`-th root of `x`, correctly rounded to `f64`
/// under round-to-nearest-even.
#[must_use]
pub fn rootn(x: f64, n: i32) -> f64 {
    rootn_round(x, n, RoundingMode::NearestEven).0
}

/// `rootn(x, n)` correctly rounded to `f64` under `mode`, with the IEEE
/// 754 status flags raised. IEEE 754-2019 §9.2 sign and domain rules are
/// inherited from pfloat's kernel.
#[must_use]
pub fn rootn_round(x: f64, n: i32, mode: RoundingMode) -> (f64, Status) {
    drive::<F64Shell>(x, mode, move |xb, w, dir| {
        xb.rootn_round(n, w, dir)
            .expect("w = PREC + guard >= 1: BuildError only on precision 0")
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const MODES: [RoundingMode; 5] = [
        RoundingMode::NearestEven,
        RoundingMode::NearestAway,
        RoundingMode::TowardZero,
        RoundingMode::TowardPositive,
        RoundingMode::TowardNegative,
    ];

    /// Kernel-exact: every mode agrees and INEXACT is clear.
    fn exact_clear(f: fn(f64, RoundingMode) -> (f64, Status), x: f64, want: f64) {
        for &mode in &MODES {
            let (got, st) = f(x, mode);
            assert_eq!(got.to_bits(), want.to_bits(), "value x={x} mode={mode:?}");
            assert!(st.is_ok() && !st.inexact(), "status x={x} mode={mode:?}");
        }
    }

    /// Mathematically exact via a composed transcendental: value exact,
    /// INEXACT may be conservatively set (ADR-0057).
    fn exact_value(f: fn(f64, RoundingMode) -> (f64, Status), x: f64, want: f64) {
        for &mode in &MODES {
            let (got, _st) = f(x, mode);
            assert_eq!(got.to_bits(), want.to_bits(), "value x={x} mode={mode:?}");
        }
    }

    #[test]
    fn kernel_exact_values_clear_inexact() {
        exact_clear(exp_round, 0.0, 1.0);
        exact_clear(ln_round, 1.0, 0.0);
        exact_clear(log2_round, 8.0, 3.0);
        exact_clear(sqrt_round, 16.0, 4.0);
        exact_clear(cbrt_round, 27.0, 3.0);
        exact_clear(cbrt_round, -27.0, -3.0);
        exact_clear(cos_round, 0.0, 1.0);
        exact_clear(sec_round, -0.0, 1.0);
    }

    #[test]
    fn composed_exact_values_round_correctly() {
        exact_value(exp10_round, 3.0, 1000.0);
        exact_value(log10_round, 1000.0, 3.0);
    }

    #[test]
    fn inexact_is_set_for_transcendentals() {
        assert!(exp_round(1.0, RoundingMode::NearestEven).1.inexact());
        assert!(ln_round(2.0, RoundingMode::NearestEven).1.inexact());
    }

    #[test]
    fn specials() {
        assert!(exp(f64::NAN).is_nan());
        assert!(exp(f64::INFINITY).is_infinite() && exp(f64::INFINITY).is_sign_positive());
        // exp(-inf) = +0 exactly: no UNDERFLOW.
        let (v, st) = exp_round(f64::NEG_INFINITY, RoundingMode::NearestEven);
        assert_eq!(v.to_bits(), 0.0f64.to_bits());
        assert!(st.is_ok() && !st.underflow());
        // exp overflow / underflow of the format.
        let (v, st) = exp_round(1000.0, RoundingMode::NearestEven);
        assert!(v.is_infinite() && st.overflow() && st.inexact());
        let (v, st) = exp_round(-1000.0, RoundingMode::NearestEven);
        assert_eq!(v.to_bits(), 0.0f64.to_bits());
        assert!(st.underflow() && st.inexact());
        // log poles / domain.
        let (v, st) = ln_round(0.0, RoundingMode::NearestEven);
        assert!(v.is_infinite() && v.is_sign_negative() && st.div_by_zero());
        let (v, st) = ln_round(-1.0, RoundingMode::NearestEven);
        assert!(v.is_nan() && st.invalid());
        // sqrt domain / sign.
        assert!(sqrt(-1.0).is_nan());
        assert_eq!(sqrt(-0.0).to_bits(), (-0.0f64).to_bits());
        // inverse domain.
        assert!(asin(2.0).is_nan());
        assert!(acosh(0.5).is_nan());
        let (v, st) = atanh_round(1.0, RoundingMode::NearestEven);
        assert!(v.is_infinite() && st.div_by_zero());
        // hyperbolic saturation.
        assert_eq!(tanh(f64::INFINITY).to_bits(), 1.0f64.to_bits());
    }

    #[test]
    fn reciprocal_trig_poles() {
        let (v, st) = csc_round(0.0, RoundingMode::NearestEven);
        assert!(v.is_infinite() && v.is_sign_positive() && st.div_by_zero());
        let (v, st) = cot_round(-0.0, RoundingMode::NearestEven);
        assert!(v.is_infinite() && v.is_sign_negative() && st.div_by_zero());
        assert!(sec(f64::INFINITY).is_nan());
    }

    #[test]
    fn sin_of_huge_input_is_bounded() {
        let (v, st) = sin_round(f64::MAX, RoundingMode::NearestEven);
        assert!(v.abs() <= 1.0 && !st.invalid() && !v.is_nan());
    }

    #[test]
    fn hypot_and_rootn() {
        assert_eq!(hypot(3.0, 4.0).to_bits(), 5.0f64.to_bits());
        assert!(hypot(f64::INFINITY, f64::NAN).is_infinite());
        assert!(hypot(f64::NAN, 2.0).is_nan());
        assert_eq!(rootn(27.0, 3).to_bits(), 3.0f64.to_bits());
        assert_eq!(rootn(-8.0, 3).to_bits(), (-2.0f64).to_bits());
        let (v, st) = rootn_round(-8.0, 2, RoundingMode::NearestEven);
        assert!(v.is_nan() && st.invalid());
        let (v, st) = rootn_round(8.0, 0, RoundingMode::NearestEven);
        assert!(v.is_nan() && st.invalid());
    }
}
