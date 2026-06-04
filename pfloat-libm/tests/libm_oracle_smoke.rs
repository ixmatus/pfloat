//! Plumbing smoke for the MPFR verification harness.
//!
//! Asserts the oracle bridge wires end to end: known values certify and
//! match the shell; domain errors and poles certify the right
//! `INVALID`/`DIV_BY_ZERO`; composed-exact results clear `INEXACT` and
//! the enclosure-derived `INEXACT` gate agrees (pf-njs5, ADR-0060); both
//! widths work. Sub-second; not the coverage gate (that is
//! `libm_smoke_gate`).

#![cfg(all(unix, feature = "differential-mpfr"))]

mod harness;

use harness::{
    certified_round, enclose, verify_input, Hw, LibmArg, LibmFnId, StatusGate, Verdict, ALL_MODES,
};
use pfloat_libm::RoundingMode;

const GATE: StatusGate = StatusGate::ValueAndDomainHard;

/// Every mode certifies and matches the shell for this unary input.
fn ok_unary_f32(f: LibmFnId, x: f32) {
    for &mode in ALL_MODES {
        let v = verify_input::<f32>(f, x.to_bits(), LibmArg::None, mode, GATE);
        assert!(
            matches!(v, Verdict::Ok),
            "f32 {f:?}({x}) mode={mode:?}: {v:?}"
        );
    }
}

fn ok_unary_f64(f: LibmFnId, x: f64) {
    for &mode in ALL_MODES {
        let v = verify_input::<f64>(f, x.to_bits(), LibmArg::None, mode, GATE);
        assert!(
            matches!(v, Verdict::Ok),
            "f64 {f:?}({x}) mode={mode:?}: {v:?}"
        );
    }
}

#[test]
fn enclosure_brackets_known_values() {
    // exp(0) = 1 exactly; the bracket's NE rounding is 1.0.
    let xf = <f32 as Hw>::lift(0.0f32.to_bits(), 64);
    let enc = enclose(LibmFnId::Exp, &xf, None, 64);
    assert_eq!(
        certified_round::<f32>(&enc, RoundingMode::NearestEven).map(f32::to_bits),
        Some(1.0f32.to_bits())
    );
    // sqrt(2) certified f32.
    let two = <f32 as Hw>::lift(2.0f32.to_bits(), 64);
    let enc = enclose(LibmFnId::Sqrt, &two, None, 64);
    let got = certified_round::<f32>(&enc, RoundingMode::NearestEven).unwrap();
    assert_eq!(got, 2.0f32.sqrt());
}

#[test]
fn normal_range_unary_clean_both_widths() {
    for &f in LibmFnId::UNARY {
        // An input inside every function's domain (0.5 is in all; trig
        // reciprocals are away from their poles there).
        ok_unary_f32(f, 0.5);
        ok_unary_f64(f, 0.5);
    }
    // A second anchor for the x >= 1 functions.
    ok_unary_f32(LibmFnId::Acosh, 1.5);
    ok_unary_f32(LibmFnId::Ln, 2.0);
    ok_unary_f32(LibmFnId::Atanh, 0.5);
    ok_unary_f64(LibmFnId::Acosh, 1.5);
}

#[test]
fn domain_errors_certify_invalid() {
    // True value NaN from a finite input: INVALID expected, and the
    // shell must agree (hard gate).
    ok_unary_f32(LibmFnId::Ln, -1.0);
    ok_unary_f32(LibmFnId::Log2, -1.0);
    ok_unary_f32(LibmFnId::Sqrt, -1.0);
    ok_unary_f32(LibmFnId::Asin, 2.0);
    ok_unary_f32(LibmFnId::Acos, 2.0);
    ok_unary_f32(LibmFnId::Acosh, 0.5);
    ok_unary_f32(LibmFnId::Atanh, 2.0);
    ok_unary_f32(LibmFnId::Log1p, -2.0);
    ok_unary_f64(LibmFnId::Ln, -1.0);
    ok_unary_f64(LibmFnId::Atanh, 2.0);
}

#[test]
fn poles_certify_div_by_zero() {
    ok_unary_f32(LibmFnId::Ln, 0.0);
    ok_unary_f32(LibmFnId::Ln, -0.0);
    ok_unary_f32(LibmFnId::Log10, 0.0);
    ok_unary_f32(LibmFnId::Log1p, -1.0);
    ok_unary_f32(LibmFnId::Atanh, 1.0);
    ok_unary_f32(LibmFnId::Atanh, -1.0);
    ok_unary_f32(LibmFnId::Cot, 0.0);
    ok_unary_f32(LibmFnId::Cot, -0.0);
    ok_unary_f32(LibmFnId::Csc, 0.0);
    ok_unary_f64(LibmFnId::Ln, 0.0);
    ok_unary_f64(LibmFnId::Cot, 0.0);
}

#[test]
fn composed_exact_values_clear_inexact_and_gate_agrees() {
    // log10(1000)=3, exp10(2)=100, exp2(10)=1024 are exact. The
    // exact-input dispatch (pf-njs5, ADR-0060) clears INEXACT, and the
    // enclosure-derived INEXACT gate confirms it: value matches AND the
    // now-hard INEXACT flag agrees.
    ok_unary_f32(LibmFnId::Log10, 1000.0);
    ok_unary_f32(LibmFnId::Exp10, 2.0);
    ok_unary_f32(LibmFnId::Exp2, 10.0);
    ok_unary_f64(LibmFnId::Log10, 1000.0);
}

#[test]
fn non_finite_inputs_clean() {
    // NaN propagates (no INVALID); inf maps per IEEE.
    for &f in LibmFnId::UNARY {
        ok_unary_f32(f, f32::NAN);
    }
    ok_unary_f32(LibmFnId::Exp, f32::INFINITY);
    ok_unary_f32(LibmFnId::Exp, f32::NEG_INFINITY);
    ok_unary_f32(LibmFnId::Atan, f32::INFINITY);
    ok_unary_f32(LibmFnId::Tanh, f32::INFINITY);
    // Trig of an infinity is qNaN + INVALID (the enclosure-derived
    // expectation must match the shell).
    ok_unary_f32(LibmFnId::Sin, f32::INFINITY);
    ok_unary_f32(LibmFnId::Cos, f32::INFINITY);
    ok_unary_f32(LibmFnId::Cot, f32::INFINITY);
    ok_unary_f32(LibmFnId::Sec, f32::INFINITY);
}

#[test]
fn binary_functions_clean() {
    // hypot(3, 4) = 5.
    for &mode in ALL_MODES {
        let v = verify_input::<f32>(
            LibmFnId::Hypot,
            3.0f32.to_bits(),
            LibmArg::HypotY(u64::from(4.0f32.to_bits())),
            mode,
            GATE,
        );
        assert!(matches!(v, Verdict::Ok), "hypot(3,4) mode={mode:?}: {v:?}");
    }
    // rootn(27, 3) = 3; rootn(-8, 3) = -2; even root of negative and
    // n = 0 are INVALID; rootn(0, -2) is a pole.
    let cases: &[(i32, f32, &str)] = &[
        (3, 27.0, "cube root"),
        (3, -8.0, "odd root of negative"),
        (2, -8.0, "even root of negative -> INVALID"),
        (0, 8.0, "n=0 -> INVALID"),
        (-2, 0.0, "rootn(0,-2) -> DIV_BY_ZERO"),
    ];
    for &(n, x, label) in cases {
        for &mode in ALL_MODES {
            let v = verify_input::<f32>(LibmFnId::Rootn(n), x.to_bits(), LibmArg::None, mode, GATE);
            assert!(
                matches!(v, Verdict::Ok),
                "rootn({x}, {n}) [{label}] mode={mode:?}: {v:?}"
            );
        }
    }
}
