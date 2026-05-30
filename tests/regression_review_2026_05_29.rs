//! Failing regression tests for the 2026-05-29 correctness review.
//!
//! Each test encodes the *correct* behaviour, so it fails on the
//! current tree and goes green when the underlying defect is fixed.
//! The defects cluster into two systemic root causes plus a handful
//! of conformance / panic nits; the tests are grouped accordingly.
//!
//! Two oracle strategies are used, both reproducible:
//!
//! - **Phase-reduction findings (root cause 1).** pfloat is wrong at
//!   every reachable precision for these (the working precision does
//!   not scale with the argument magnitude), so self-reference is not
//!   a valid oracle. These assert closeness to an external reference
//!   computed with `mpmath` 1.4.1 at 4000 bits on the *exactly
//!   representable* `2^k` input (so pfloat's input and the oracle's
//!   input are bit-identical). Reference values are quoted inline.
//!
//! - **Cancellation findings (root cause 2).** Here the high-precision
//!   path *is* correct, so the test needs no external oracle: it
//!   asserts the precision-refinement self-consistency that correct
//!   rounding guarantees, namely `f(x)@target == round(f(x)@HIGH,
//!   target)`. The inputs sit a few hundred bits from each function's
//!   zero (verified near-zero with mpmath), where the relative
//!   half-width under-bounds the absolute cancellation error.
//!
//! Run: `cargo test --test regression_review_2026_05_29 \
//!        --features std,fmt,big,integrals,airy,bessel`

#![cfg(all(
    feature = "big",
    feature = "integrals",
    feature = "airy",
    feature = "bessel"
))]

use pfloat::{BigFloat, RoundingMode, Sign};

const NE: RoundingMode = RoundingMode::NearestEven;

/// Exact `2^k` (a single-bit mantissa, representable at any precision).
fn pow2(k: u32, prec: u32) -> BigFloat {
    let two = BigFloat::try_from_i64_exact(2, prec).unwrap();
    let mut x = BigFloat::try_from_i64_exact(1, prec).unwrap();
    for _ in 0..k {
        x = x.mul(&two, NE).0;
    }
    x
}

/// `|got - ref| / |ref| < 1e-12` (about 40 bits): far tighter than any
/// of the defects (which give relative error near 1, often wrong sign)
/// and far looser than correct rounding to 53 bits (~1e-16), so it
/// cleanly separates broken from fixed without being brittle.
fn assert_close(label: &str, got: &BigFloat, reference: &str) {
    let r = BigFloat::parse_str(reference, 200, NE).unwrap().0;
    assert_eq!(
        got.is_sign_negative(),
        r.is_sign_negative(),
        "{label}: sign mismatch (got {got}, expected {reference})"
    );
    let diff = got.sub(&r, NE).0.abs();
    let rel = diff.div(&r.abs(), NE).0;
    let tol = BigFloat::parse_str("1e-12", 200, NE).unwrap().0;
    assert!(
        matches!(
            rel.partial_cmp(&tol).0,
            Some(core::cmp::Ordering::Less | core::cmp::Ordering::Equal)
        ),
        "{label}: got {got}, expected ~{reference} (relative error {rel})"
    );
}

/// Correct rounding is precision-refinement stable: the result at
/// `target` must equal the result computed at `HIGH` then rounded to
/// `target`. The high-precision path is correct for these inputs, so a
/// disagreement is a misround at `target`.
fn assert_refinement_stable(label: &str, x: &BigFloat, f: impl Fn(&BigFloat, u32) -> BigFloat) {
    // Comfortably above the test inputs' cancellation depth (~262 bits
    // for the li case, the deepest) so f@HIGH is correct well past 53
    // bits, while staying in the kernels' fast regime.
    const HIGH: u32 = 320;
    let direct = f(x, 53);
    let refined = f(x, HIGH).round_to_precision(53, NE).unwrap().0;
    assert_eq!(
        direct.partial_cmp(&refined).0,
        Some(core::cmp::Ordering::Equal),
        "{label}: direct@53 = {direct} but round(@{HIGH}->53) = {refined} (misround at p=53)"
    );
}

// =====================================================================
// Root cause 1: working / phase precision does not scale with |x|.
// The Ziv driver certifies garbage on iteration 1 because its
// relative half-width is taken around the already-wrong value.
// =====================================================================

// trig_reduce::reduce (src/math/trig_reduce.rs:102-105): mul_prec omits
// e_x, so x*(2/pi) is rounded to <= 2048 significant bits regardless of
// the exponent. For 2^2048 <= |x| the fractional part vanishes and the
// reduced argument collapses to 0.
#[test]
fn rc1_trig_large_argument_is_correctly_reduced() {
    let x = pow2(3000, 53);
    let (s, sst) = x.sin(NE);
    let (c, cst) = x.cos(NE);
    // Current (buggy) behaviour: sin = 0, cos = 1, status OK.
    assert_close(
        "sin(2^3000)",
        &s,
        "-0.9023960377032291054512045038415469221809",
    );
    assert_close(
        "cos(2^3000)",
        &c,
        "-0.4309076364344362774210502315402372443805",
    );
    // Whatever the reduction does, the result of a transcendental on a
    // representable input is inexact, never an exact 0/1 with OK.
    assert!(
        sst.inexact() && cst.inexact(),
        "sin/cos of 2^3000 cannot be exact"
    );
}

// airy_asymptotic_neg (src/math/airy.rs:571): working = target + 32, so
// zeta = (2/3) t^(3/2) ~ 2^300 loses all fractional bits at large |x|;
// the phase phi = zeta - pi/4 reaching sin/cos is already garbage.
#[test]
fn rc1_airy_large_negative_argument_has_correct_sign() {
    let x = pow2(200, 300).negated();
    let (ai, st) = x.ai_round(53, NE).unwrap();
    // Current (buggy) behaviour: +4.24e-16 (wrong sign), status OK-ish.
    assert_close(
        "Ai(-2^200)",
        &ai,
        "-3.004443173224147928442938076132835769249e-16",
    );
    assert!(st.inexact(), "Ai of a representable input is inexact");
}

// bessel_j_asymptotic: the Hankel phase omega = x - m*pi/2 - pi/4 is
// computed at working capped at target+512, and the Ziv loop certifies
// on iteration 1 (working = target+64) before it can grow. Both the
// hard cap and the magnitude-independent start are wrong.
#[test]
fn rc1_bessel_large_argument_has_correct_sign() {
    let (j4, _) = pow2(400, 53).j0_round(53, NE).unwrap();
    assert_close(
        "J0(2^400)",
        &j4,
        "-3.033405919027553113526656815748938467038e-61",
    );
    // 2^800 is beyond the target+512 cap entirely: cannot be recovered
    // without removing the cap.
    let (j8, _) = pow2(800, 53).j0_round(53, NE).unwrap();
    assert_close(
        "J0(2^800)",
        &j8,
        "-7.131908349607747958111886972283760654196e-122",
    );
}

// =====================================================================
// Root cause 2: Ziv's relative half-width d = |y|*2^(guard-working)
// under-bounds the ABSOLUTE error of a near-zero result formed by
// cancellation of O(1)-magnitude operands. error_guard = 24 covers
// only 24 bits of cancellation; the reflection branches and near-zero
// integrals have unbounded cancellation depth.
// =====================================================================

// li(x) = Ei(ln x); near the Ramanujan-Soldner constant the inner Ei
// cancels. Double composition makes this the worst of the family
// (wrong sign, ~3800x at p=53).
#[test]
fn rc2_li_near_its_zero_is_correctly_rounded() {
    // ~1e-79 from the li zero; deep in the cancellation band at p=53.
    let x = BigFloat::parse_str(
        "1.4513692348833810502839684858920274494930322836480158630930045576624255957545178",
        700,
        NE,
    )
    .unwrap()
    .0;
    assert_refinement_stable("li", &x, |x, p| x.li_round(p, NE).unwrap().0);
}

// lgamma reflection ln|Gamma(x)| = ln(pi) - ln|sin(pi x)| - lgamma(1-x)
// cancels near the negative-axis roots of |Gamma| = 1.
#[test]
fn rc2_lgamma_near_negative_root_is_correctly_rounded() {
    let x = BigFloat::parse_str(
        "-2.45702473822080062303945414765117954323659789544059009815495",
        700,
        NE,
    )
    .unwrap()
    .0;
    assert_refinement_stable("lgamma", &x, |x, p| x.lgamma_round(p, NE).unwrap().0);
}

// digamma reflection psi(x) = psi(1-x) - pi*cot(pi x) cancels near the
// negative-axis zeros of psi.
#[test]
fn rc2_digamma_near_negative_zero_is_correctly_rounded() {
    let x = BigFloat::parse_str(
        "-0.504083008264455409258269304533302498955385439742794143169891",
        700,
        NE,
    )
    .unwrap()
    .0;
    assert_refinement_stable("digamma", &x, |x, p| x.digamma_round(p, NE).unwrap().0);
}

// Ci(x) = γ + ln|x| + Σ … cancels near its real zero (~0.61650549),
// the same relative-half-width defect as li (single composition rather
// than li's Ei(ln x) double composition).
#[test]
fn rc2_ci_near_its_zero_is_correctly_rounded() {
    let x = BigFloat::parse_str(
        "0.6165054856207162337971104041001727475394958981816653305721207211247532",
        700,
        NE,
    )
    .unwrap()
    .0;
    assert_refinement_stable("ci", &x, |x, p| x.ci_round(p, NE).unwrap().0);
}

// expm1/log1p cancellation boost is capped (inner_w = min(w+cancel,
// w+1024), src/math/expm1.rs:152). For x below ~2^-(target+1088) the
// 1-subtraction collapses to exactly 0 and half_width(0)=0 certifies
// it on iteration 1. expm1(2^-2000) must round to ~2^-2000, not 0.
#[test]
fn rc2_expm1_tiny_x_does_not_collapse_to_zero() {
    let x = BigFloat::try_from_i64_exact(1, 53)
        .unwrap()
        .div(&pow2(2000, 53), NE)
        .0; // 2^-2000, exact
    let (em, _) = x.expm1(NE);
    assert!(
        !em.is_zero(),
        "expm1(2^-2000) collapsed to 0 (true ~8.71e-603)"
    );
    assert!(
        em.is_sign_positive(),
        "expm1 of a positive value is positive"
    );
}

// The collapse re-breaks tanh, which ADR-0050 believed it had fixed via
// expm1: tanh(positive tiny) returns -0 (wrong value AND wrong sign).
#[test]
fn rc2_tanh_tiny_x_keeps_sign_and_magnitude() {
    let x = BigFloat::try_from_i64_exact(1, 53)
        .unwrap()
        .div(&pow2(2000, 53), NE)
        .0;
    let (t, _) = x.tanh(NE);
    assert!(
        !t.is_zero(),
        "tanh(2^-2000) collapsed to 0 (true ~8.71e-603)"
    );
    assert!(t.is_sign_positive(), "tanh(+x) must be positive, got {t}");
}

// =====================================================================
// Arithmetic core and IEEE conformance.
// =====================================================================

// huge_gap_short_circuit (src/ops/addsub.rs:433-470) routes the larger
// operand through the pipeline with pre_sticky only, ignoring the
// residue *direction*. For opposite-sign subtraction the true value is
// just below `large`, so TowardZero must reach the predecessor.
#[test]
fn arith_huge_gap_directed_rounding_is_correct() {
    let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
    let tiny = one.div(&pow2(200, 53), NE).0; // 2^-200
                                              // 1 - 2^-200 is strictly in (1 - 2^-53, 1); TowardZero floors it.
    let (r, _) = one.sub(&tiny, RoundingMode::TowardZero);
    assert_eq!(
        r.partial_cmp(&one).0,
        Some(core::cmp::Ordering::Less),
        "TowardZero(1 - 2^-200) must be < 1 (the predecessor), got {r}"
    );
}

// div.rs:118 matches (_, Zero) before any (Infinity, _) arm, so Inf/0
// raises DIV_BY_ZERO. IEEE 754-2019 §7.3 requires a finite dividend;
// Inf/0 is exact +Inf with no flag.
#[test]
fn conformance_inf_div_zero_raises_no_flag() {
    let inf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
    let zero = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
    let (q, st) = inf.div(&zero, NE);
    assert!(
        q.is_infinite() && q.is_sign_positive(),
        "Inf/0 value is +Inf"
    );
    assert!(!st.div_by_zero(), "Inf/0 must not raise DIV_BY_ZERO (§7.3)");
    assert!(st.is_ok(), "Inf/0 raises no exception flag");
}

// cmp.rs:116-121: max(x, sNaN) returns Status::INVALID but never calls
// auto_raise, so the std thread-local flag bag disagrees with the
// returned status (min does raise it).
#[cfg(feature = "std")]
#[test]
fn conformance_max_snan_raises_threadlocal_invalid() {
    let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
    let snan = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
    pfloat::flags::clear();
    let (_m, st) = one.max(&snan);
    assert!(
        st.invalid(),
        "max(x, sNaN) returns INVALID in the explicit status"
    );
    assert!(
        pfloat::flags::test().invalid(),
        "max(x, sNaN) must also raise the thread-local INVALID (min does)"
    );
    pfloat::flags::clear();
}

// =====================================================================
// Panic-safety: i64/u32 exponent arithmetic is not hardened the way
// mul/div are (pf-rnc i128 saturation). These panic in debug / wrap to
// a silently wrong result in release on reachable saturated-exponent
// values. cargo test runs debug, so each currently panics (= fails).
// =====================================================================

// fmt.rs:256: `bits + num_p2` overflows u32 for a value with a large
// binary exponent (2^(2^40) is reachable by 40 squarings). Display /
// to_decimal_string must not panic on any finite value.
#[test]
fn panic_fmt_large_exponent_does_not_overflow() {
    let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
    let mut v = two;
    for _ in 0..40 {
        v = v.mul(&v, NE).0; // v = 2^(2^40), binary exponent 1099511627776
    }
    let s = v.to_decimal_string(17, NE);
    assert!(
        s.bytes().any(|b| b.is_ascii_digit()),
        "to_decimal_string of a large-exponent value must return digits, not panic"
    );
}

// addsub.rs:297 / sqrt.rs:148: `e - p + 1` underflows i64 for an operand
// whose exponent saturated to ~i64::MIN (the reciprocal of a saturated-
// exponent value). mul/div were i128-hardened (pf-rnc); add/sub/sqrt
// were not. Must not panic on any finite operands.
#[test]
fn panic_addsub_sqrt_saturated_exponent_is_bounded() {
    let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
    let mut big = two;
    for _ in 0..200 {
        let (sq, st) = big.mul(&big, NE); // exponent saturates to i64::MAX
        big = sq;
        if st.overflow() {
            break;
        }
    }
    let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
    let (tiny, _) = one.div(&big, NE); // exponent ~ i64::MIN + 1
    let (s, _) = tiny.add(&tiny, NE);
    assert!(
        !s.is_nan(),
        "add of a near-i64::MIN-exponent operand must not panic"
    );
    let (q, _) = tiny.sqrt(NE);
    assert!(
        !q.is_nan(),
        "sqrt of a near-i64::MIN-exponent operand must not panic"
    );
}

// trig_reduce.rs:73: `e_x + slack` overflows i64 for a value whose
// exponent saturated to i64::MAX (reachable by repeated squaring). The
// module doc promises qNaN + INVALID for out-of-range |x|, not a panic.
#[test]
fn panic_trig_saturated_exponent_returns_invalid() {
    let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
    let mut big = two;
    // Square until the exponent saturates (OVERFLOW flag) at i64::MAX.
    for _ in 0..200 {
        let (sq, st) = big.mul(&big, NE);
        big = sq;
        if st.overflow() {
            break;
        }
    }
    let (r, st) = big.sin(NE);
    assert!(
        r.is_nan() && st.invalid(),
        "sin of an out-of-range |x| must be qNaN + INVALID, not a panic"
    );
}
