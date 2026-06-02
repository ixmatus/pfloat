//! Correctly-rounded `f32` elementary and reciprocal functions.
//!
//! Each entry widens its `f32` argument to an exact `BigFloat`, evaluates
//! pfloat's correctly-rounded kernel through the outer Ziv loop in
//! [`crate::round`], and rounds back to `f32` straight onto the format
//! grid. The bare form rounds to nearest-even; the `_round` form takes an
//! explicit [`RoundingMode`] and returns the IEEE 754 [`Status`] flags the
//! call raised.

use pfloat::BigFloat;

use crate::round::{drive, unary, F32Shell};
use crate::{RoundingMode, Status};

// --- exponentials (saturation fast-path: see `crate::saturate`) ---
unary!(
    f32,
    F32Shell,
    exp,
    exp_round,
    direct_sat,
    crate::saturate::sat_exp,
    "`exp`"
);
unary!(
    f32,
    F32Shell,
    exp2,
    exp2_round,
    result_sat,
    crate::saturate::sat_exp2,
    "`exp2`"
);
unary!(
    f32,
    F32Shell,
    exp10,
    exp10_round,
    result_sat,
    crate::saturate::sat_exp10,
    "`exp10`"
);
unary!(
    f32,
    F32Shell,
    expm1,
    expm1_round,
    result_sat,
    crate::saturate::sat_expm1,
    "`expm1`"
);

// --- logarithms ---
unary!(f32, F32Shell, ln, ln_round, result, "`ln`");
unary!(f32, F32Shell, log2, log2_round, result, "`log2`");
unary!(f32, F32Shell, log10, log10_round, result, "`log10`");
unary!(f32, F32Shell, log1p, log1p_round, result, "`log1p`");

// --- root ---
unary!(f32, F32Shell, sqrt, sqrt_round, result, "`sqrt`");
unary!(f32, F32Shell, cbrt, cbrt_round, result, "`cbrt`");

// --- circular ---
unary!(f32, F32Shell, sin, sin_round, result, "`sin`");
unary!(f32, F32Shell, cos, cos_round, result, "`cos`");
unary!(f32, F32Shell, tan, tan_round, result, "`tan`");
unary!(f32, F32Shell, cot, cot_round, result, "`cot`");
unary!(f32, F32Shell, sec, sec_round, result, "`sec`");
unary!(f32, F32Shell, csc, csc_round, result, "`csc`");

// --- inverse circular ---
unary!(f32, F32Shell, asin, asin_round, result, "`asin`");
unary!(f32, F32Shell, acos, acos_round, result, "`acos`");
unary!(f32, F32Shell, atan, atan_round, result, "`atan`");

// --- hyperbolic (saturation fast-path: see `crate::saturate`) ---
unary!(
    f32,
    F32Shell,
    sinh,
    sinh_round,
    result_sat,
    crate::saturate::sat_sinh,
    "`sinh`"
);
unary!(
    f32,
    F32Shell,
    cosh,
    cosh_round,
    result_sat,
    crate::saturate::sat_cosh,
    "`cosh`"
);
unary!(
    f32,
    F32Shell,
    tanh,
    tanh_round,
    result_sat,
    crate::saturate::sat_tanh,
    "`tanh`"
);

// --- inverse hyperbolic ---
unary!(f32, F32Shell, asinh, asinh_round, result, "`asinh`");
unary!(f32, F32Shell, acosh, acosh_round, result, "`acosh`");
unary!(f32, F32Shell, atanh, atanh_round, result, "`atanh`");

/// `hypot(x, y) = sqrt(x² + y²)`, correctly rounded to `f32` under
/// round-to-nearest-even.
#[must_use]
pub fn hypot(x: f32, y: f32) -> f32 {
    hypot_round(x, y, RoundingMode::NearestEven).0
}

/// `hypot(x, y)` correctly rounded to `f32` under `mode`, with the IEEE
/// 754 status flags raised. IEEE 754-2019 §9.2.1 special cases (an
/// infinity argument yields `+∞` even if the other is NaN) are inherited
/// from pfloat's kernel.
#[must_use]
pub fn hypot_round(x: f32, y: f32, mode: RoundingMode) -> (f32, Status) {
    let yb = BigFloat::from_f32(y);
    drive::<F32Shell>(x, mode, move |xb, w, dir| {
        xb.hypot_round(&yb, w, dir)
            .expect("w = PREC + guard >= 1: BuildError only on precision 0")
    })
}

/// `rootn(x, n)`, the real `n`-th root of `x`, correctly rounded to `f32`
/// under round-to-nearest-even.
#[must_use]
pub fn rootn(x: f32, n: i32) -> f32 {
    rootn_round(x, n, RoundingMode::NearestEven).0
}

/// `rootn(x, n)` correctly rounded to `f32` under `mode`, with the IEEE
/// 754 status flags raised. IEEE 754-2019 §9.2 sign and domain rules
/// (even root of a negative is `qNaN` + `INVALID`, odd root preserves the
/// sign) are inherited from pfloat's kernel.
#[must_use]
pub fn rootn_round(x: f32, n: i32, mode: RoundingMode) -> (f32, Status) {
    drive::<F32Shell>(x, mode, move |xb, w, dir| {
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

    /// A value the kernel computes exactly: every mode agrees and INEXACT
    /// is clear (the exact short-circuit in the driver fires).
    fn exact_clear(f: fn(f32, RoundingMode) -> (f32, Status), x: f32, want: f32) {
        for &mode in &MODES {
            let (got, st) = f(x, mode);
            assert_eq!(got.to_bits(), want.to_bits(), "value x={x} mode={mode:?}");
            assert!(st.is_ok() && !st.inexact(), "status x={x} mode={mode:?}");
        }
    }

    /// A value that is mathematically exact and correctly returned in every
    /// mode, but reached through a composed-transcendental kernel that may
    /// conservatively report INEXACT (the kernel cannot cheaply prove the
    /// result is exact). The correctly-rounded VALUE is asserted here; the
    /// flag conservatism is documented in ADR-0057.
    fn exact_value(f: fn(f32, RoundingMode) -> (f32, Status), x: f32, want: f32) {
        for &mode in &MODES {
            let (got, _st) = f(x, mode);
            assert_eq!(got.to_bits(), want.to_bits(), "value x={x} mode={mode:?}");
        }
    }

    #[test]
    fn kernel_exact_values_clear_inexact() {
        exact_clear(exp_round, 0.0, 1.0);
        exact_clear(ln_round, 1.0, 0.0);
        exact_clear(log2_round, 8.0, 3.0); // log2 has an exact power-of-two path
        exact_clear(sqrt_round, 4.0, 2.0);
        exact_clear(cbrt_round, 8.0, 2.0);
        exact_clear(cbrt_round, -8.0, -2.0);
        exact_clear(sin_round, 0.0, 0.0);
        exact_clear(cos_round, 0.0, 1.0);
        exact_clear(cos_round, -0.0, 1.0);
        exact_clear(sec_round, 0.0, 1.0);
        exact_clear(tanh_round, 0.0, 0.0);
        exact_clear(atan_round, 0.0, 0.0);
    }

    #[test]
    fn composed_exact_values_round_correctly() {
        // log10 = ln/ln(10), exp10 = exp(x·ln 10), exp2 = exp(x·ln 2): all
        // compose transcendentals, so INEXACT may be conservatively set even
        // though the result is exact. The value must still be exact.
        exact_value(log10_round, 1000.0, 3.0);
        exact_value(exp10_round, 2.0, 100.0);
        exact_value(exp2_round, 10.0, 1024.0);
    }

    #[test]
    fn cbrt_signed_zero() {
        assert_eq!(cbrt(0.0).to_bits(), 0.0f32.to_bits());
        assert_eq!(cbrt(-0.0).to_bits(), (-0.0f32).to_bits());
    }

    #[test]
    fn inexact_is_set_for_transcendentals() {
        let (_v, st) = exp_round(1.0, RoundingMode::NearestEven);
        assert!(st.inexact());
        let (_v, st) = sqrt_round(2.0, RoundingMode::NearestEven);
        assert!(st.inexact());
    }

    #[test]
    fn nan_propagates() {
        let x = f32::NAN;
        assert!(exp(x).is_nan());
        assert!(ln(x).is_nan());
        assert!(sin(x).is_nan());
        assert!(sqrt(x).is_nan());
        assert!(cbrt(x).is_nan());
    }

    #[test]
    fn exp_infinities_and_saturation() {
        assert!(exp(f32::INFINITY).is_infinite() && exp(f32::INFINITY).is_sign_positive());
        // exp(-inf) = +0 EXACTLY: no UNDERFLOW (the value is exact, not tiny).
        let (v, st) = exp_round(f32::NEG_INFINITY, RoundingMode::NearestEven);
        assert_eq!(v.to_bits(), 0.0f32.to_bits());
        assert!(st.is_ok() && !st.underflow());
        // Far overflow / underflow of the format.
        let (v, st) = exp_round(1000.0, RoundingMode::NearestEven);
        assert!(v.is_infinite() && v.is_sign_positive() && st.overflow() && st.inexact());
        let (v, st) = exp_round(-1000.0, RoundingMode::NearestEven);
        assert_eq!(v.to_bits(), 0.0f32.to_bits());
        assert!(st.underflow() && st.inexact());
    }

    #[test]
    fn log_domain_and_poles() {
        for z in [0.0f32, -0.0] {
            let (v, st) = ln_round(z, RoundingMode::NearestEven);
            assert!(v.is_infinite() && v.is_sign_negative() && st.div_by_zero());
        }
        let (v, st) = ln_round(-1.0, RoundingMode::NearestEven);
        assert!(v.is_nan() && st.invalid());
    }

    #[test]
    fn sqrt_domain_and_sign() {
        let (v, st) = sqrt_round(-1.0, RoundingMode::NearestEven);
        assert!(v.is_nan() && st.invalid());
        // sqrt(-0) = -0.
        assert_eq!(sqrt(-0.0).to_bits(), (-0.0f32).to_bits());
    }

    #[test]
    fn inverse_domain_errors() {
        let (v, st) = asin_round(2.0, RoundingMode::NearestEven);
        assert!(v.is_nan() && st.invalid());
        let (v, st) = acos_round(2.0, RoundingMode::NearestEven);
        assert!(v.is_nan() && st.invalid());
        let (v, st) = acosh_round(0.5, RoundingMode::NearestEven);
        assert!(v.is_nan() && st.invalid());
        let (v, st) = atanh_round(2.0, RoundingMode::NearestEven);
        assert!(v.is_nan() && st.invalid());
        // atanh(1) = +inf + DIV_BY_ZERO (DLMF §4.37).
        let (v, st) = atanh_round(1.0, RoundingMode::NearestEven);
        assert!(v.is_infinite() && v.is_sign_positive() && st.div_by_zero());
    }

    #[test]
    fn hyperbolic_saturation() {
        assert!(cosh(f32::INFINITY).is_infinite());
        assert_eq!(tanh(f32::INFINITY).to_bits(), 1.0f32.to_bits());
        assert_eq!(tanh(f32::NEG_INFINITY).to_bits(), (-1.0f32).to_bits());
    }

    #[test]
    fn reciprocal_trig_poles() {
        // cot/csc are odd: poles at ±0 take the sign of the zero.
        let (v, st) = cot_round(0.0, RoundingMode::NearestEven);
        assert!(v.is_infinite() && v.is_sign_positive() && st.div_by_zero());
        let (v, st) = cot_round(-0.0, RoundingMode::NearestEven);
        assert!(v.is_infinite() && v.is_sign_negative() && st.div_by_zero());
        let (v, st) = csc_round(0.0, RoundingMode::NearestEven);
        assert!(v.is_infinite() && v.is_sign_positive() && st.div_by_zero());
        // cot/sec/csc of an infinity are qNaN + INVALID.
        for f in [cot as fn(f32) -> f32, sec, csc] {
            assert!(f(f32::INFINITY).is_nan());
        }
    }

    #[test]
    fn sin_of_huge_input_is_bounded() {
        // f32::MAX ≈ 2^128, far inside pfloat's Payne-Hanek table budget,
        // so the result is a bounded finite value with NO INVALID.
        let (v, st) = sin_round(f32::MAX, RoundingMode::NearestEven);
        assert!(v.abs() <= 1.0 && !st.invalid() && !v.is_nan());
    }

    #[test]
    fn hypot_known_and_specials() {
        assert_eq!(hypot(3.0, 4.0).to_bits(), 5.0f32.to_bits());
        assert_eq!(hypot(5.0, 12.0).to_bits(), 13.0f32.to_bits());
        let (_v, st) = hypot_round(1.0, 1.0, RoundingMode::NearestEven);
        assert!(st.inexact());
        // §9.2.1: hypot of an infinity is +inf even if the other is NaN.
        assert!(hypot(f32::INFINITY, f32::NAN).is_infinite());
        assert!(hypot(f32::NAN, f32::INFINITY).is_infinite());
        assert!(hypot(f32::NAN, 2.0).is_nan());
    }

    #[test]
    fn rootn_known_and_specials() {
        assert_eq!(rootn(27.0, 3).to_bits(), 3.0f32.to_bits());
        assert_eq!(rootn(16.0, -2).to_bits(), 0.25f32.to_bits());
        assert_eq!(rootn(-8.0, 3).to_bits(), (-2.0f32).to_bits());
        // Even root of a negative, and n = 0, are domain errors.
        let (v, st) = rootn_round(-8.0, 2, RoundingMode::NearestEven);
        assert!(v.is_nan() && st.invalid());
        let (v, st) = rootn_round(8.0, 0, RoundingMode::NearestEven);
        assert!(v.is_nan() && st.invalid());
        // rootn(±0, negative) is a pole.
        let (v, st) = rootn_round(0.0, -2, RoundingMode::NearestEven);
        assert!(v.is_infinite() && st.div_by_zero());
    }
}
