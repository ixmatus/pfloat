//! Regression tests for epic pf-8iji review remediation R4 (ADR-0115):
//! status and sign fidelity in complex divide, multiply, and csqrt.
//!
//! - pf-pz9r: a zero-dividend part's sign follows the inputs, not the mode.
//! - pf-yprp: csqrt's negative-real-axis imaginary part rounds on the correct
//!   side under directed modes (negate after the mirrored directed round).
//! - pf-hdq1: a quiet NaN propagates without INVALID; a signaling NaN raises it.
//! - pf-bv2i: an exact quotient reports OK; an exponent-saturating enclosure
//!   carries the saturation flag rather than a wrong value with a clean flag.

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Sign};
use pfloat_complex::Complex;

fn bf(n: i64, p: u32) -> BigFloat {
    BigFloat::try_from_i64_exact(n, p).unwrap()
}
fn nz(p: u32) -> BigFloat {
    BigFloat::try_new_zero(Sign::Negative, p).unwrap()
}
fn pow2(k: i64, p: u32) -> BigFloat {
    bf(1, p).scale_by_pow2(k).0
}
fn c(re: BigFloat, im: BigFloat) -> Complex<BigFloat> {
    Complex::new(re, im)
}
fn eq(v: &BigFloat, w: &BigFloat) -> bool {
    matches!(v.partial_cmp(w).0, Some(Ordering::Equal))
}

const NE: RoundingMode = RoundingMode::NearestEven;

// pf-pz9r: (-0 - 0i)/(3 + 0i) has a zero real dividend; componentwise IEEE
// division gives re = -0 (sign from the input, not the rounding mode).
#[test]
fn pz9r_zero_dividend_sign_follows_inputs() {
    let z = c(nz(64), nz(64));
    let w = c(bf(3, 64), bf(0, 64));
    let (q, _) = z.div(&w, NE);
    assert!(
        q.re.is_zero() && q.re.is_sign_negative(),
        "re must be -0, got {} (neg={})",
        q.re,
        q.re.is_sign_negative()
    );
    // The imaginary part (a cancelling difference) stays +0 under NE.
    assert!(q.im.is_zero() && q.im.is_sign_positive());
}

// A positive zero dividend keeps +0, and a negative divisor still leaves the
// non-negative denominator, so the quotient sign follows the numerator alone.
#[test]
fn pz9r_positive_zero_dividend_stays_positive() {
    let z = c(bf(0, 64), bf(0, 64));
    let w = c(bf(3, 64), bf(0, 64));
    let (q, _) = z.div(&w, NE);
    assert!(q.re.is_zero() && q.re.is_sign_positive());
}

// z/z with z = (1 + 2^-199) + i at p=200 is exactly 1 + 0i; ADR-0090 requires
// an exact quotient to report OK, not INEXACT (pf-bv2i part a).
#[test]
fn bv2i_a_exact_self_quotient_is_ok() {
    let p = 200u32;
    let re = bf(1, p + 4).add(&pow2(-199, p + 4), NE).0;
    let z = c(re, bf(1, p + 4));
    let (q, s) = z.div(&z, NE);
    assert!(eq(&q.re, &bf(1, p)), "re must be exactly 1, got {}", q.re);
    assert!(q.im.is_zero(), "im must be exactly 0, got {}", q.im);
    assert!(
        s.is_ok(),
        "an exact quotient must report OK (no INEXACT), got {s:?}"
    );
}

// (1 + 1i)/(c + ci) with c = 2^-(2^62 + 100): forming c² + d² underflows the
// i64 exponent, so the enclosure is unsound. The status must carry the
// saturation flag rather than returning a wrong value with a clean OK flag
// (pf-bv2i part b).
#[test]
fn bv2i_b_saturation_carries_status() {
    let p = 64u32;
    let k = -(4_611_686_018_427_387_904_i64 + 100); // -(2^62 + 100)
    let cval = pow2(k, p);
    let z = c(bf(1, p), bf(1, p));
    let w = c(cval.clone(), cval);
    let (_, s) = z.div(&w, NE);
    assert!(
        !s.is_ok() && (s.underflow() || s.overflow()),
        "an exponent-saturating enclosure must not report OK, got {s:?}"
    );
}

// (qNaN + 1i)/(2 + 0i): a quiet NaN propagates SILENTLY -- no INVALID
// (pf-hdq1 part a).
#[test]
fn hdq1_a_quiet_nan_propagates_without_invalid() {
    let p = 64u32;
    let qnan = BigFloat::try_new_quiet_nan(Sign::Positive, p, &[]).unwrap();
    let z = c(qnan, bf(1, p));
    let w = c(bf(2, p), bf(0, p));
    let (q, s) = z.div(&w, NE);
    assert!(q.re.is_nan() && q.im.is_nan());
    assert!(
        !s.invalid(),
        "a quiet NaN operand must not raise INVALID, got {s:?}"
    );
}

// A signaling-NaN divisor part MUST raise INVALID even where a quiet NaN would
// not (the IEEE 754 signaling rule).
#[test]
fn hdq1_a_signaling_nan_divisor_raises_invalid() {
    let p = 64u32;
    let snan = BigFloat::try_new_signaling_nan(Sign::Positive, p, &[]).unwrap();
    let z = c(bf(1, p), bf(1, p));
    let w = c(snan, bf(2, p));
    let (_, s) = z.div(&w, NE);
    assert!(s.invalid(), "a signaling NaN operand must raise INVALID");
}

// mul((inf + sNaN i), (2 + 3i)): the §G.5.1 recovery must not swallow the
// INVALID a signaling NaN operand raises (pf-hdq1 part b).
#[test]
fn hdq1_b_signaling_nan_in_mul_recovery_raises_invalid() {
    let p = 64u32;
    let inf = BigFloat::try_new_infinity(Sign::Positive, p).unwrap();
    let snan = BigFloat::try_new_signaling_nan(Sign::Positive, p, &[]).unwrap();
    let z = c(inf, snan);
    let w = c(bf(2, p), bf(3, p));
    let (_, s) = z.mul(&w, NE);
    assert!(
        s.invalid(),
        "a signaling NaN operand must raise INVALID through mul recovery"
    );
}

// A quiet NaN through the same recovery must NOT raise INVALID.
#[test]
fn hdq1_b_quiet_nan_in_mul_recovery_no_invalid() {
    let p = 64u32;
    let inf = BigFloat::try_new_infinity(Sign::Positive, p).unwrap();
    let qnan = BigFloat::try_new_quiet_nan(Sign::Positive, p, &[]).unwrap();
    let z = c(inf, qnan);
    let w = c(bf(2, p), bf(3, p));
    let (_, s) = z.mul(&w, NE);
    assert!(
        !s.invalid(),
        "a quiet NaN operand must not raise INVALID, got {s:?}"
    );
}

// Adversarial (pf-bv2i): a genuinely inexact quotient must STILL report
// INEXACT. (1 + 2i)/(3 + 4i) = 11/25 + (2/25)i, non-dyadic, so both parts
// round; the new exactness certificate must not clear INEXACT here.
#[test]
fn bv2i_inexact_quotient_still_inexact() {
    let z = c(bf(1, 64), bf(2, 64));
    let w = c(bf(3, 64), bf(4, 64));
    let (_, s) = z.div(&w, NE);
    assert!(s.inexact(), "a non-dyadic quotient must report INEXACT");
}

// Adversarial (pf-hdq1): a signaling NaN in multiply with NO infinity present
// (so §G.5.1 recovery returns None) must still raise INVALID via the naive
// fused product.
#[test]
fn hdq1_signaling_nan_mul_without_infinity_invalid() {
    let p = 64u32;
    let snan = BigFloat::try_new_signaling_nan(Sign::Positive, p, &[]).unwrap();
    let z = c(snan, bf(1, p));
    let w = c(bf(2, p), bf(3, p));
    let (_, s) = z.mul(&w, NE);
    assert!(s.invalid(), "a signaling NaN operand must raise INVALID");
}

// pf-yprp: csqrt(-2 - 0i) under a directed mode must round the imaginary part
// on the side the mode asks for. The imaginary part is -sqrt(2); under
// TowardPositive it must be >= -sqrt(2) (the value nearer +inf), i.e. the
// magnitude rounded TowardNegative then negated.
#[cfg(feature = "exp-log")]
#[test]
fn yprp_csqrt_negative_axis_directed_imaginary_correct_side() {
    use RoundingMode::{TowardNegative as TN, TowardPositive as TP};
    let p = 64u32;
    let z = c(bf(-2, p), nz(p));
    // TowardPositive: expected im = -(sqrt(2) rounded TowardNegative).
    let (q_tp, _) = z.sqrt(TP);
    let expected_tp = bf(2, p).sqrt_round(p, TN).unwrap().0.negated();
    assert!(
        eq(&q_tp.im, &expected_tp),
        "TP im must be -(sqrt2 @TN) = {expected_tp}, got {}",
        q_tp.im
    );
    // TowardNegative: expected im = -(sqrt(2) rounded TowardPositive).
    let (q_tn, _) = z.sqrt(TN);
    let expected_tn = bf(2, p).sqrt_round(p, TP).unwrap().0.negated();
    assert!(
        eq(&q_tn.im, &expected_tn),
        "TN im must be -(sqrt2 @TP) = {expected_tn}, got {}",
        q_tn.im
    );
    // Correct-side ordering: TP result is >= TN result (TP nearer +inf).
    assert!(matches!(
        q_tp.im.partial_cmp(&q_tn.im).0,
        Some(Ordering::Greater)
    ));
}

// The positive imaginary branch (y = +0) is not negated, so a directed mode
// applies to the magnitude directly.
#[cfg(feature = "exp-log")]
#[test]
fn yprp_csqrt_negative_axis_positive_imaginary_directed() {
    use RoundingMode::TowardPositive as TP;
    let p = 64u32;
    let z = c(bf(-2, p), bf(0, p));
    let (q, _) = z.sqrt(TP);
    let expected = bf(2, p).sqrt_round(p, TP).unwrap().0;
    assert!(
        eq(&q.im, &expected),
        "y=+0 TP im must be +sqrt2 @TP = {expected}, got {}",
        q.im
    );
}
