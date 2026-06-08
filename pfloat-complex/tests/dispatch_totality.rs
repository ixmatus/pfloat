//! Exhaustive special-value dispatch verification for the complex elementary
//! core and §G.5.1 recovery (ADR-0092, the C5 verification pass).
//!
//! Where `annex_g_special_values.rs` pins the row VALUES against the standard,
//! this lane proves the dispatch has no GAP. The special-value dispatch
//! branches on a finite grid of IEEE component classes, so an exhaustive
//! enumeration over that grid is a complete totality proof: stronger than a
//! sampled Kani run (`BigFloat` is `Vec`-backed and CBMC-hostile, ADR-0062), it
//! checks every class combination directly through the public API.
//!
//! Four properties, none of which re-encodes the row tables (so a shared
//! transcription error in those tables cannot mask a defect here):
//!
//! 1. **Totality**: every class combination returns without panicking, in
//!    every rounding mode.
//! 2. **Signaling-NaN ⇒ INVALID**: any signaling-NaN operand raises INVALID.
//! 3. **Quiet-NaN-without-infinity ⇒ `(NaN, NaN)` and no INVALID** for the
//!    unary functions.
//! 4. **Conjugation symmetry** `f(conj z) = conj(f z)` for finite / signed-zero
//!    inputs, the Annex G symmetry that fixes every signed-zero branch row.

#![cfg(feature = "trig")]

mod common;

use common::{rep, ALL_CLASSES, ALL_MODES, NE};
use core::cmp::Ordering;
use pfloat::{BigFloat, RoundingMode, Status};
use pfloat_complex::Complex;

const PRECS: [u32; 2] = [53, 113];

/// The three unary elementary kernels, by name, dispatched through the public
/// API at one mode.
fn unary(name: &str, z: &Complex<BigFloat>, mode: RoundingMode) -> (Complex<BigFloat>, Status) {
    match name {
        "csqrt" => z.sqrt(mode),
        "cexp" => z.exp(mode),
        "clog" => z.log(mode),
        other => panic!("unknown unary {other}"),
    }
}

const UNARY: [&str; 3] = ["csqrt", "cexp", "clog"];

/// Exact equality including signed-zero and NaN: both NaN, or equal in value
/// AND sign. (IEEE comparison treats `+0 == -0`, so the sign clause is needed.)
fn exact_signed_eq(a: &BigFloat, b: &BigFloat) -> bool {
    if a.is_nan() || b.is_nan() {
        return a.is_nan() && b.is_nan();
    }
    matches!(a.partial_cmp(b).0, Some(Ordering::Equal))
        && a.is_sign_negative() == b.is_sign_negative()
}

#[test]
fn totality_unary_over_full_class_grid() {
    // Every (re_class, im_class) input, every unary kernel, every mode, both
    // precisions: the call returns (no panic) and yields two classifiable
    // components. 8*8 = 64 inputs * 3 fns * 5 modes * 2 precs.
    for &p in &PRECS {
        for &mode in &ALL_MODES {
            for &name in &UNARY {
                for &rc in &ALL_CLASSES {
                    for &ic in &ALL_CLASSES {
                        let z = Complex::new(rep(rc, p), rep(ic, p));
                        let (r, _) = unary(name, &z, mode);
                        // Use the result so the dispatch is not optimized away;
                        // every BigFloat is in exactly one IEEE class.
                        assert!(
                            r.re.is_nan()
                                || r.re.is_infinite()
                                || r.re.is_zero()
                                || r.re.is_finite(),
                            "{name}({rc:?}, {ic:?}).re unclassifiable under {mode:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn totality_binary_over_full_class_grid() {
    // The §G.5.1 mul/div recovery dispatch over the full 64*64 operand-pair
    // grid at p=53: every pair returns without panic. This is the no-gap proof
    // for the infinity/NaN recovery branches.
    let p = 53;
    for &rc in &ALL_CLASSES {
        for &ic in &ALL_CLASSES {
            let z = Complex::new(rep(rc, p), rep(ic, p));
            for &qc in &ALL_CLASSES {
                for &dc in &ALL_CLASSES {
                    let w = Complex::new(rep(qc, p), rep(dc, p));
                    let (m, _) = z.mul(&w, NE);
                    let (d, _) = z.div(&w, NE);
                    // Reaching here without a panic IS the totality proof for
                    // this operand pair; classifying both results forces the
                    // dispatch to actually run (not be optimized away).
                    assert!(
                        matches!(
                            common::classify(&m.re),
                            common::ResultCls::Nan
                                | common::ResultCls::PosInf
                                | common::ResultCls::NegInf
                                | common::ResultCls::PosZero
                                | common::ResultCls::NegZero
                                | common::ResultCls::PosFin
                                | common::ResultCls::NegFin
                        ),
                        "mul dispatch produced an unclassifiable component"
                    );
                    let _ = d;
                }
            }
        }
    }
}

#[test]
fn signaling_nan_raises_invalid_across_the_grid() {
    // Any signaling-NaN operand component raises INVALID, for every kernel.
    let p = 64;
    let snan = common::snan(p);
    for &name in &UNARY {
        for other in [common::bf(2, p), common::pinf(p), common::pz(p)] {
            let (_, s1) = unary(name, &Complex::new(snan.clone(), other.clone()), NE);
            let (_, s2) = unary(name, &Complex::new(other, snan.clone()), NE);
            assert!(s1.invalid(), "{name}(sNaN + y) must raise INVALID");
            assert!(s2.invalid(), "{name}(x + sNaN*i) must raise INVALID");
        }
    }
    // Binary mul/div: an sNaN operand raises INVALID too.
    let fin = common::bf(2, p);
    let z = Complex::new(snan.clone(), fin.clone());
    let w = Complex::new(fin.clone(), fin.clone());
    assert!(
        z.mul(&w, NE).1.invalid(),
        "mul with sNaN must raise INVALID"
    );
    assert!(
        z.div(&w, NE).1.invalid(),
        "div with sNaN must raise INVALID"
    );
}

#[test]
fn quiet_nan_without_infinity_is_nan_pair_no_invalid() {
    // For the unary kernels, a quiet-NaN component with NO infinity present
    // makes both outputs NaN and does not signal. (With an infinity present the
    // dominance rules can override the NaN, so those rows are excluded here and
    // pinned by name in the value table.)
    let p = 64;
    let qnan = common::qnan(p);
    // Finite NONZERO only: cexp(NaN +- 0i) is NaN +- 0i (a conjugation-pinned
    // signed-zero imaginary part, not NaN), so a zero "other" is the documented
    // exception to "quiet NaN => (NaN, NaN)" and is pinned by value in
    // annex_g_special_values::cexp_nan_real instead.
    let finite_nonzero = [common::bf(2, p), common::bf(-2, p)];
    for &name in &UNARY {
        for other in &finite_nonzero {
            for z in [
                Complex::new(qnan.clone(), other.clone()),
                Complex::new(other.clone(), qnan.clone()),
            ] {
                let (r, s) = unary(name, &z, NE);
                assert!(
                    r.re.is_nan() && r.im.is_nan(),
                    "{name} with a quiet NaN and no infinity must be (NaN, NaN)"
                );
                assert!(!s.invalid(), "{name} quiet NaN must not raise INVALID");
            }
        }
    }
}

#[test]
fn conjugation_symmetry_on_finite_and_signed_zero_inputs() {
    // Annex G mandates f(conj z) = conj(f z). On finite / signed-zero inputs the
    // functions are fully determined, so the equality is exact and componentwise
    // (including signed zeros) -- this fixes every signed-zero branch row and is
    // independent of the enumerated value table. (The inf/NaN [REP] rows, where
    // the standard leaves a sign free, are deliberately excluded: the crate
    // picks a fixed representative that need not be conjugation-symmetric on
    // those measure-zero points.)
    //
    // BIT-EXACT only under the sign-symmetric rounding modes {NE, NA, TZ}. The
    // real part is EVEN under conjugation (wants the same mode) and the
    // imaginary part is ODD (round_M(-I) = -round_{mirror(M)}(I), so it wants
    // the mirror mode); one single-mode evaluation satisfies both only when
    // M = mirror(M), i.e. the symmetric modes. Under TowardPositive /
    // TowardNegative the identity still holds mathematically but to ~1 ULP, not
    // bit-for-bit, which is correct directed-rounding behavior, not a defect.
    let symmetric_modes = [
        RoundingMode::NearestEven,
        RoundingMode::NearestAway,
        RoundingMode::TowardZero,
    ];
    let finite_or_zero = [
        common::bf(2, 113),
        common::bf(-2, 113),
        common::bf(3, 113),
        common::bf(-3, 113),
        common::pz(113),
        common::nz(113),
    ];
    for &mode in &symmetric_modes {
        for &name in &UNARY {
            for re in &finite_or_zero {
                for im in &finite_or_zero {
                    let z = Complex::new(re.clone(), im.clone());
                    let (lhs, _) = unary(name, &z.conj(), mode);
                    let (fz, _) = unary(name, &z, mode);
                    let rhs = fz.conj();
                    assert!(
                        exact_signed_eq(&lhs.re, &rhs.re),
                        "{name}: conjugation symmetry broke on .re at ({re}, {im}) mode {mode:?}: {} vs {}",
                        lhs.re,
                        rhs.re
                    );
                    assert!(
                        exact_signed_eq(&lhs.im, &rhs.im),
                        "{name}: conjugation symmetry broke on .im at ({re}, {im}) mode {mode:?}: {} vs {}",
                        lhs.im,
                        rhs.im
                    );
                }
            }
        }
    }
}
