//! Regression guards for the 2026-06-10 workspace deep review
//! (epic pf-8iji, remediation arc R1: the certified-wrong-answer
//! family).
//!
//! Each test encodes the *correct* behaviour for one confirmed
//! defect. Every test began red against the defect it guards and
//! lands in the same commit as its fix, so the lane records one
//! expected count per defect bucket rather than an aggregate floor.
//!
//! Oracle strategy mirrors `regression_review_2026_05_29.rs`:
//! external references are computed with `mpmath` 1.4.1 at 4000 bits
//! on *exactly representable* inputs (single-bit or few-bit
//! mantissas, so pfloat's input and the oracle's input are
//! bit-identical) and quoted inline; where the high-precision path
//! is correct once fixed, precision-refinement self-consistency
//! (`f(x)@target == round(f(x)@HIGH, target)`) backs it up.
//!
//! Run: `cargo test --test regression_review_2026_06_10 \
//!        --features std,fmt,big,agm,trig,specials,zeta`

#![cfg(all(
    feature = "big",
    feature = "agm",
    feature = "trig",
    feature = "specials",
    feature = "zeta"
))]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode};

const NE: RoundingMode = RoundingMode::NearestEven;

/// `|got - ref| / |ref| < 1e-12` (about 40 bits): far tighter than
/// the defects this lane guards (relative errors from 5e-3 up to
/// "wrong class entirely") and far looser than correct rounding at
/// 53 bits, so it cleanly separates broken from fixed without being
/// brittle.
fn assert_close(label: &str, got: &BigFloat, reference: &str) {
    let r = BigFloat::parse_str(reference, 200, NE).unwrap().0;
    assert_eq!(
        got.is_sign_negative(),
        r.is_sign_negative(),
        "{label}: sign mismatch (got {got}, want {reference})"
    );
    let (diff, _) = got.sub(&r, NE);
    let (rel, _) = diff.abs().div(&r.abs(), NE);
    let bound = BigFloat::parse_str("1e-12", 200, NE).unwrap().0;
    assert_eq!(
        rel.partial_cmp(&bound).0,
        Some(Ordering::Less),
        "{label}: relative error {rel} exceeds 1e-12 (got {got}, want {reference})"
    );
}

/// `value * 2^k`, exact (scaling never touches the mantissa).
fn scaled(value: i64, prec: u32, k: i64) -> BigFloat {
    let (x, status) = BigFloat::try_from_i64_exact(value, prec)
        .unwrap()
        .scale_by_pow2(k);
    assert!(status.is_ok(), "scaled({value}, {prec}, {k}) not exact");
    x
}

// ---------------------------------------------------------------
// pf-ddfl: agm convergence floor was absolute (-w - 4), not
// relative to the operand magnitude. Small operands tripped it
// before the first Gauss iteration and the kernel returned the
// arithmetic mean (0.5% relative error here) with Status OK.
// ---------------------------------------------------------------

/// agm(2^-300, 3*2^-302): mpmath 1.4.1 @4000 bits, inputs exact.
/// The broken kernel returned exactly (a + b) / 2 = 0.875 * 2^-300.
#[test]
fn agm_small_operands_iterates_to_the_agm() {
    let a = scaled(1, 53, -300);
    let b = scaled(3, 53, -302);
    let (r, status) = a.agm(&b, NE);
    assert_close(
        "agm(2^-300, 3*2^-302)",
        &r,
        "4.273399828000648542805471530695713670719e-91",
    );
    // The true AGM of unequal operands is not representable at 53
    // bits; OK was part of the defect.
    assert!(status.inexact(), "agm small-operand status must be INEXACT");
    // And it must not be the arithmetic mean bitwise.
    let (sum, _) = a.add(&b, NE);
    let (am, _) = sum.scale_by_pow2(-1);
    assert_ne!(
        r.total_cmp(&am),
        Ordering::Equal,
        "agm returned the arithmetic mean"
    );
}

/// Same defect family at the opposite scale: for large operands the
/// absolute floor 2^(-w-4) was unreachable, so the loop always ran
/// all 64 iterations (wasted work, result still correct). The
/// relative criterion must keep this case right.
/// agm(2^300, 3*2^298): mpmath 1.4.1 @4000 bits, inputs exact.
#[test]
fn agm_large_operands_unchanged_by_relative_floor() {
    let a = scaled(1, 53, 300);
    let b = scaled(3, 53, 298);
    let (r, status) = a.agm(&b, NE);
    assert_close(
        "agm(2^300, 3*2^298)",
        &r,
        "1.773253911834204859984452477441122355622e+90",
    );
    assert!(status.inexact());
}

/// Precision-refinement self-consistency on the defect input:
/// correct rounding implies agm@53 == round(agm@2000 -> 53).
#[test]
fn agm_small_operands_refinement_consistency() {
    let a = scaled(1, 53, -300);
    let b = scaled(3, 53, -302);
    let (r53, _) = a.agm(&b, NE);
    let a_hi = scaled(1, 2000, -300);
    let b_hi = scaled(3, 2000, -302);
    let (r_hi, _) = a_hi.agm_round(&b_hi, 2000, NE).unwrap();
    let (r_hi_53, _) = r_hi.round_to_precision(53, NE).unwrap();
    assert_eq!(
        r53.total_cmp(&r_hi_53),
        Ordering::Equal,
        "agm@53 disagrees with round(agm@2000 -> 53): {r53} vs {r_hi_53}"
    );
}
