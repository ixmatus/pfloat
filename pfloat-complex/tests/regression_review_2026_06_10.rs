//! Regression guard for the 2026-06-10 workspace deep review,
//! pfloat-complex slice (epic pf-8iji: pf-qm8a clog real-part band
//! collapse). Each test began red and lands with its fix (ADR-0100).

#![cfg(feature = "trig")]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode};
use pfloat_complex::Complex;

const NE: RoundingMode = RoundingMode::NearestEven;

fn assert_bit_exact(label: &str, got: &BigFloat, reference: &str, p: u32) {
    let expected = BigFloat::parse_str(reference, p, NE).unwrap().0;
    assert_eq!(
        got.total_cmp(&expected),
        Ordering::Equal,
        "{label}: got {got}, want {expected}"
    );
}

/// pf-qm8a: clog(1 + 2^-545 i) at p64. The hypot bracket straddles
/// |z| = 1 to a depth (2·545 bits) the static GUARDS schedule
/// (≤ p + 1024) cannot resolve, and the exhausted loop silently
/// returned the unconverged end: re = +0 with INEXACT where the
/// truth ≈ 2^-1091 is representable. References: mpmath 1.4.1
/// @4400 bits at the exact dyadic input.
#[test]
fn clog_near_unit_circle_resolves_the_band() {
    let re_in = BigFloat::try_from_i64_exact(1, 64).unwrap();
    let (im_in, s) = BigFloat::try_from_i64_exact(1, 64)
        .unwrap()
        .scale_by_pow2(-545);
    assert!(s.is_ok());
    let z = Complex::new(re_in, im_in);
    let (w, st) = z.log(NE);
    assert!(
        !w.re().is_zero(),
        "clog real part collapsed to a signed zero"
    );
    assert_bit_exact(
        "re(clog(1 + 2^-545 i))",
        &w.re().clone(),
        "3.76942173645970568982367548269822214023635714e-329",
        64,
    );
    assert_bit_exact(
        "im(clog(1 + 2^-545 i))",
        &w.im().clone(),
        "8.68265136517608391969809927449437537856454907e-165",
        64,
    );
    assert!(st.inexact());
}

/// Control: the measure-zero-adjacent shallow case the old schedule
/// already handled (depth well inside p + 1024).
#[test]
fn clog_near_unit_circle_shallow_control() {
    let re_in = BigFloat::try_from_i64_exact(1, 64).unwrap();
    let (im_in, s) = BigFloat::try_from_i64_exact(1, 64)
        .unwrap()
        .scale_by_pow2(-100);
    assert!(s.is_ok());
    let z = Complex::new(re_in, im_in);
    let (w, _) = z.log(NE);
    // re = 0.5·ln(1 + 2^-200) ≈ 2^-201: nonzero and tiny.
    assert!(!w.re().is_zero());
    assert!(w.re().is_finite());
}

/// The adversarial verifier's refutation of this slice's first
/// draft: at depth >= 576 the inner scalar hypot exhausts ITS OWN
/// Ziv cap and returns a falsely-exact 1, so the outer bracket
/// "converged exactly" on [0,0] at the first guard and the
/// depth-scaled growth never ran (component status OK — worse than
/// the original defect). The exponent carries the depth here, so
/// the input is cheap to build at p64. No nontrivial dyadic point
/// lies on the unit circle, so a both-ends-1 bracket with nonzero
/// components is always a lie; the fixed loop treats it as
/// unresolved. mpmath 1.4.1 @9000 bits.
#[test]
fn clog_exponent_encoded_depth_resolves() {
    let re_in = BigFloat::try_from_i64_exact(1, 64).unwrap();
    let (im_in, s) = BigFloat::try_from_i64_exact(1, 64)
        .unwrap()
        .scale_by_pow2(-2000);
    assert!(s.is_ok());
    let z = Complex::new(re_in, im_in);
    let (w, st) = z.log(NE);
    assert!(!w.re().is_zero(), "deep band collapsed again");
    assert_bit_exact(
        "re(clog(1 + 2^-2000 i))",
        &w.re().clone(),
        "3.79303935173368928611560265184341857540759251e-1205",
        64,
    );
    assert_bit_exact(
        "im(clog(1 + 2^-2000 i))",
        &w.im().clone(),
        "8.70980981621721667557619549477887229585910374e-603",
        64,
    );
    assert!(st.inexact());
}
