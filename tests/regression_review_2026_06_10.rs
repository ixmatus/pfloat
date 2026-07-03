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
//!        --features std,fmt,big,agm,trig,integrals,zeta`
//! (`integrals` implies `specials`; arc R3 adds Ei/li rows that need
//! it. Both CI jobs that exercise this lane — the full-feature-union
//! `cargo test` and the release MPFR job — already enable it.)

#![cfg(all(
    feature = "big",
    feature = "agm",
    feature = "trig",
    feature = "integrals",
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

// ---------------------------------------------------------------
// pf-7z66: exp at the exponent ceiling, three confirmed failures.
// (a) deep underflow returned a garbage Normal near i64::MIN with
// INEXACT only; (b) k = round(x/ln2) wrapped `2^63 as i64` to
// i64::MIN inside the reduction, returning a tiny garbage Normal
// where the truth is near 2^(i64::MAX); (c) just below the wrap
// window the result must be a representable finite at exponent
// i64::MAX, not +inf. mpmath 1.4.1 @400 bits pins every window
// classification and mantissa quoted below; the inputs parse
// exactly at the stated precisions.
// ---------------------------------------------------------------

const ALL_MODES: [RoundingMode; 5] = [
    RoundingMode::NearestEven,
    RoundingMode::NearestAway,
    RoundingMode::TowardZero,
    RoundingMode::TowardPositive,
    RoundingMode::TowardNegative,
];

/// Smallest positive value at precision `p` (`2^(i64::MIN)`).
fn min_pos(p: u32) -> BigFloat {
    use pfloat::Sign;
    BigFloat::try_new_zero(Sign::Positive, p)
        .unwrap()
        .next_up()
        .0
}

/// Largest finite value at precision `p` (all-ones mantissa at
/// exponent `i64::MAX`).
fn max_finite(p: u32) -> BigFloat {
    use pfloat::Sign;
    BigFloat::try_new_infinity(Sign::Positive, p)
        .unwrap()
        .next_down()
        .0
}

/// `m_str * 2^k` at precision 53 under `mode`: the 40-digit decimal
/// mantissa rounds to 53 bits, then the power-of-two scale is exact.
fn mantissa_at_pow2(m_str: &str, k: i64, mode: RoundingMode) -> BigFloat {
    let (m, _) = BigFloat::parse_str(m_str, 53, mode).unwrap();
    let (v, s) = m.scale_by_pow2(k);
    assert!(s.is_ok(), "scale by 2^{k} must be exact");
    v
}

/// (a) exp(-1e300): certain deep underflow. The truth is below half
/// of `MinPos`, so every mode except `TowardPositive` rounds to +0;
/// `TowardPositive` must round up to `MinPos`. `UNDERFLOW|INEXACT` either
/// way. The broken kernel returned a garbage Normal near `i64::MIN`
/// with INEXACT only.
#[test]
fn exp_deep_underflow_is_mode_aware_zero() {
    // "-1e300" is not dyadic; the 53-bit parse is INEXACT and the
    // parsed neighbour is equally deep in the certain-underflow
    // region, so exactness is irrelevant here.
    let (x, _) = BigFloat::parse_str("-1e300", 53, NE).unwrap();
    for mode in ALL_MODES {
        let (r, st) = x.exp_round(53, mode);
        assert!(st.underflow(), "exp(-1e300) {mode:?}: UNDERFLOW missing");
        assert!(st.inexact(), "exp(-1e300) {mode:?}: INEXACT missing");
        if matches!(mode, RoundingMode::TowardPositive) {
            assert_eq!(
                r.total_cmp(&min_pos(53)),
                Ordering::Equal,
                "exp(-1e300) TP must be MinPos, got {r}"
            );
        } else {
            assert!(
                r.is_zero() && !r.is_sign_negative(),
                "exp(-1e300) {mode:?} must be +0, got {r}"
            );
        }
    }
}

/// Certain overflow at exp(1e300) and just past the top window
/// (x = 6393154322601327830 > 2^63·ln2 ≈ ...829.894): +inf for the
/// to-nearest and upward modes, `MaxFinite` for the inward ones,
/// `OVERFLOW|INEXACT` always.
#[test]
fn exp_overflow_is_mode_aware_inf() {
    // "1e300" parses inexactly (any neighbour is equally deep past
    // the rim); the integer probe just past 2^63·ln2 is exact.
    for s in ["1e300", "6393154322601327830"] {
        let (x, _) = BigFloat::parse_str(s, 64, NE).unwrap();
        for mode in ALL_MODES {
            let (r, st) = x.exp_round(53, mode);
            assert!(st.overflow(), "exp({s}) {mode:?}: OVERFLOW missing");
            assert!(st.inexact(), "exp({s}) {mode:?}: INEXACT missing");
            if matches!(
                mode,
                RoundingMode::TowardZero | RoundingMode::TowardNegative
            ) {
                assert_eq!(
                    r.total_cmp(&max_finite(53)),
                    Ordering::Equal,
                    "exp({s}) {mode:?} must be MaxFinite, got {r}"
                );
            } else {
                assert!(r.is_infinite(), "exp({s}) {mode:?} must be +inf, got {r}");
            }
        }
    }
}

/// (b)+(c) the top window: both probes have `floor(x/ln2)` = `i64::MAX`,
/// so the truth is a representable finite `m·2^(i64::MAX)`. mpmath
/// 1.4.1 @400 bits, inputs exact at 80 bits.
#[test]
fn exp_top_window_returns_the_representable_finite() {
    let cases = [
        (
            // (c) k = i64::MAX without wrap; the broken kernel
            // returned +inf. The .5 literal is dyadic: exact at 80.
            "6393154322601327829.5",
            "1.348282864605751640420850372405533755632",
        ),
        (
            // (b) round(x/ln2) = 2^63 wrapped to i64::MIN; truth is
            // at i64::MAX. The .7 literal is NOT dyadic, so the
            // oracle is evaluated at the 80-bit parsed value
            // (mpmath at mp.prec=80 for the input, 400 for exp) —
            // the bit-identical-input rule.
            "6393154322601327829.7",
            "1.646791383993419741410683242114440875528",
        ),
    ];
    for (xs, ms) in cases {
        let (x, _) = BigFloat::parse_str(xs, 80, NE).unwrap();
        for mode in [NE, RoundingMode::TowardZero] {
            let (r, st) = x.exp_round(53, mode);
            let expected = mantissa_at_pow2(ms, i64::MAX, mode);
            assert_eq!(
                r.total_cmp(&expected),
                Ordering::Equal,
                "exp({xs}) {mode:?}: got {r}, want {expected}"
            );
            assert!(st.inexact(), "exp({xs}) {mode:?}: INEXACT missing");
            assert!(
                !st.overflow() && !st.underflow(),
                "exp({xs}) {mode:?}: representable finite must not flag range, got {st:?}"
            );
        }
    }
}

/// Bottom window: `floor(x/ln2)` = `i64::MIN` exactly; the truth is the
/// representable normal `m·2^(i64::MIN)` with no UNDERFLOW.
#[test]
fn exp_bottom_window_returns_the_representable_finite() {
    let (x, sp) = BigFloat::parse_str("-6393154322601327829.5", 80, NE).unwrap();
    assert!(sp.is_ok());
    let (r, st) = x.exp_round(53, NE);
    let expected = mantissa_at_pow2("1.483368254913493617047065800911430569035", i64::MIN, NE);
    assert_eq!(
        r.total_cmp(&expected),
        Ordering::Equal,
        "exp bottom window: got {r}, want {expected}"
    );
    assert!(st.inexact());
    assert!(
        !st.underflow(),
        "representable normal at the floor must not flag UNDERFLOW"
    );
}

/// The underflow sliver: `floor(x/ln2)` = `i64::MIN` − 1, so the truth
/// lies in [`MinPos`/2, `MinPos`) — above the to-nearest midpoint
/// (mantissa 2^(t−floor t) ∈ (1,2) strictly). Nearest and upward
/// modes give `MinPos`; inward modes give +0. `UNDERFLOW|INEXACT`.
#[test]
fn exp_underflow_sliver_is_mode_aware() {
    let (x, _) = BigFloat::parse_str("-6393154322601327830.2", 80, NE).unwrap();
    for mode in ALL_MODES {
        let (r, st) = x.exp_round(53, mode);
        assert!(st.underflow(), "sliver {mode:?}: UNDERFLOW missing");
        assert!(st.inexact(), "sliver {mode:?}: INEXACT missing");
        if matches!(
            mode,
            RoundingMode::TowardZero | RoundingMode::TowardNegative
        ) {
            assert!(r.is_zero(), "sliver {mode:?} must be +0, got {r}");
        } else {
            assert_eq!(
                r.total_cmp(&min_pos(53)),
                Ordering::Equal,
                "sliver {mode:?} must be MinPos, got {r}"
            );
        }
    }
}

/// One step deeper (floor = `i64::MIN` − 2): below half of `MinPos`, so
/// only `TowardPositive` rounds away from +0.
#[test]
fn exp_just_below_sliver_rounds_to_zero() {
    let (x, sp) = BigFloat::parse_str("-6393154322601327831", 64, NE).unwrap();
    assert!(sp.is_ok());
    let (r_ne, st_ne) = x.exp_round(53, NE);
    assert!(r_ne.is_zero(), "NE must be +0, got {r_ne}");
    assert!(st_ne.underflow() && st_ne.inexact());
    let (r_tp, _) = x.exp_round(53, RoundingMode::TowardPositive);
    assert_eq!(r_tp.total_cmp(&min_pos(53)), Ordering::Equal);
}

/// Adversarial verification finding D1: the certified-floor bracket
/// cap must scale with the input precision. x = −`RD_1100(ln2)·2^63`
/// is an exact 1100-bit dyadic one part in ~2^1037 below
/// `i64::MIN`·ln2 in magnitude, so `floor(x/ln2)` = `i64::MIN` and the
/// truth is the representable normal (1 + ~2^-1037)·`2^(i64::MIN)`.
/// A fixed q = 1024 bracket cannot separate the floor from
/// `i64::MIN` − 1 and the fall-through dispatched the SLIVER: +0
/// under `TowardZero` (wrong value), `MinPos` under `TowardPositive`
/// (1 ulp low), spurious UNDERFLOW under `NearestEven`.
#[test]
fn exp_certified_floor_cap_scales_with_input_precision() {
    let two = BigFloat::try_from_i64_exact(2, 1100).unwrap();
    let (ln2_rd, _) = two.ln_round(1100, RoundingMode::TowardNegative).unwrap();
    let (mag, s) = ln2_rd.scale_by_pow2(63);
    assert!(s.is_ok());
    let x = mag.negated();

    let (r_ne, st_ne) = x.exp_round(53, NE);
    assert_eq!(
        r_ne.total_cmp(&min_pos(53)),
        Ordering::Equal,
        "NE must be MinPos (truth = (1+eps)·2^MIN), got {r_ne}"
    );
    assert!(
        !st_ne.underflow(),
        "representable normal at the floor must not flag UNDERFLOW"
    );
    assert!(st_ne.inexact());

    let (r_tz, _) = x.exp_round(53, RoundingMode::TowardZero);
    assert_eq!(
        r_tz.total_cmp(&min_pos(53)),
        Ordering::Equal,
        "TZ must be MinPos (truth strictly above it), got {r_tz}"
    );

    let (r_tp, _) = x.exp_round(53, RoundingMode::TowardPositive);
    let next_up_minpos = min_pos(53).next_up().0;
    assert_eq!(
        r_tp.total_cmp(&next_up_minpos),
        Ordering::Equal,
        "TP must be nextUp(MinPos), got {r_tp}"
    );
}

/// Adversarial verification finding D2: a certified-rounding carry
/// in the top window must dispatch as overflow, not reach
/// `scale_by_pow2`'s clamp after the carry replaced the mantissa.
/// x = `RD_130(ln2)·2^63` (exact 130-bit dyadic): the truth is
/// (2 − ~2^-72)·`2^(i64::MAX)`, strictly above `MaxFinite` at p53, so
/// the to-nearest and upward modes overflow to +inf while
/// `TowardZero` gives `MaxFinite`. The broken compose returned a
/// non-monotone 1.0·`2^(i64::MAX)`.
#[test]
fn exp_top_window_carry_overflows_mode_aware() {
    // floor(ln2·2^130), verified against mpmath 1.4.1 @300 bits.
    let (m, sp) = BigFloat::parse_str("943463052902053176551776571056617937597", 130, NE).unwrap();
    assert!(sp.is_ok(), "the 130-bit integer must parse exactly");
    let (x, s) = m.scale_by_pow2(-67);
    assert!(s.is_ok());

    for mode in [NE, RoundingMode::NearestAway, RoundingMode::TowardPositive] {
        let (r, st) = x.exp_round(53, mode);
        assert!(r.is_infinite(), "carry {mode:?} must be +inf, got {r}");
        assert!(st.overflow() && st.inexact(), "carry {mode:?}: {st:?}");
    }
    // TowardZero cannot carry, and per IEEE 754-2019 §7.4 it does
    // NOT overflow here: the unbounded-exponent TZ rounding of the
    // truth (2 − 2^-72)·2^MAX is exactly MaxFinite, which does not
    // exceed MaxFinite. INEXACT only — unlike the certain-overflow
    // dispatch (truth ≥ 2^(MAX+1)), where even TZ's unbounded
    // rounding exceeds MaxFinite and the flag is due.
    let (r_tz, st_tz) = x.exp_round(53, RoundingMode::TowardZero);
    assert_eq!(
        r_tz.total_cmp(&max_finite(53)),
        Ordering::Equal,
        "carry TZ must be MaxFinite, got {r_tz}"
    );
    assert!(
        st_tz.inexact() && !st_tz.overflow(),
        "carry TZ is INEXACT only (§7.4), got {st_tz:?}"
    );

    // Monotonicity guard: exp at the lower window probe must not
    // exceed this input's NE result class (the broken compose gave
    // 1.0·2^MAX here, BELOW the 1.348·2^MAX of the smaller input).
    let (x5, _) = BigFloat::parse_str("6393154322601327829.5", 80, NE).unwrap();
    let (r5, _) = x5.exp_round(53, NE);
    let (r_big, _) = x.exp_round(53, NE);
    assert!(
        matches!(r5.partial_cmp(&r_big).0, Some(Ordering::Less)),
        "monotonicity: exp(smaller) {r5} must be < exp(larger) {r_big}"
    );
}

// ---------------------------------------------------------------
// pf-smcb / pf-rylv / pf-wmv7: one mechanism, three kernels. The
// Ziv first-iteration error model (relative half-width
// 2^-(w-guard)) is violated when the evaluation cancels against
// structure carried in the input, and the interval test certifies
// a wrong value. ln just below 1 cancelled ln(2x) against ln2
// (2^15 ulps wrong); asin near 1 amplified the x^2 rounding error
// through 1-x^2 (1.3e29 ulps@400 wrong); lgamma/digamma had the
// RC2 cancellation boost on the negative reflection branch only,
// leaving the positive-branch roots (1, 2, and digamma's 1.46...)
// unprotected. All references mpmath 1.4.1 @4000 bits at the exact
// dyadic inputs; expected values are bit-exact (parse the quoted
// decimal at the target precision under the same mode).
// ---------------------------------------------------------------

/// Bit-exact comparison against a quoted high-precision decimal:
/// parsing it at `p` under `mode` yields the correctly rounded
/// target value (the quoted digits carry ~3x the target bits).
fn assert_bit_exact(label: &str, got: &BigFloat, reference: &str, p: u32, mode: RoundingMode) {
    let expected = BigFloat::parse_str(reference, p, mode).unwrap().0;
    assert_eq!(
        got.total_cmp(&expected),
        Ordering::Equal,
        "{label}: got {got}, want {expected}"
    );
}

/// pf-smcb: ln(1 - 2^-80) at p100 -> 53. The broken kernel
/// returned -8.2718061255904621e-25 (2^15 ulps wrong, certified).
#[test]
fn ln_just_below_one_resolves_the_cancellation() {
    let one = BigFloat::try_from_i64_exact(1, 100).unwrap();
    let (t, _) = scaled(1, 100, -80).round_to_precision(100, NE).unwrap();
    let (x, sx) = one.sub(&t, NE);
    assert!(sx.is_ok(), "1 - 2^-80 must be exact at p100");
    for mode in [NE, RoundingMode::TowardZero] {
        let (r, st) = x.ln_round(53, mode).unwrap();
        assert_bit_exact(
            "ln(1-2^-80)",
            &r,
            "-8.27180612553027674871409034183845745366854817e-25",
            53,
            mode,
        );
        assert!(st.inexact());
    }
    // Control on the other side of 1 (the e = 0 path, already
    // correct before the fix): ln(1 + 2^-80).
    let (xp, _) = one.add(&t, NE);
    let (rp, _) = xp.ln_round(53, NE).unwrap();
    assert_bit_exact(
        "ln(1+2^-80)",
        &rp,
        "8.27180612553027674871408349956079961764769405e-25",
        53,
        NE,
    );
}

/// pf-rylv, the REAL reproducer: the 1-x^2 amplification needs a
/// DENSE delta. x = 1 − RN150(π)·2^-202 (exact at p400, span 350
/// bits) makes x² span ~701 bits, past every early Ziv working
/// precision; the x² rounding error, amplified ~2^100 through the
/// 1−x² cancellation, certified a value ~2^18 ulps@400 wrong with
/// INEXACT (verified by run pre-fix). The review's named input
/// 1 − 2^-200 was a misadjudication: a single-bit delta gives the
/// sparse x² = 1 − 2^-199 + 2^-400, exactly representable at the
/// first working precision, and the recorded review output was in
/// fact correctly rounded (limb-exact vs mpmath RN400) — the
/// reviewer compared Display-truncated digits. That input stays
/// below as a control. References: mpmath 1.4.1 @4000 bits at the
/// exact dyadic inputs, quoted to 140 digits.
#[test]
fn asin_near_one_resolves_the_amplification() {
    // RN150(pi) = 1120957716564506572603712206968581818470252692 · 2^-148.
    let (pi150, sp) =
        BigFloat::parse_str("1120957716564506572603712206968581818470252692", 150, NE).unwrap();
    assert!(sp.is_ok());
    let (delta, sd) = pi150.scale_by_pow2(-148 - 202);
    assert!(sd.is_ok());
    let one = BigFloat::try_from_i64_exact(1, 400).unwrap();
    let (x, sx) = one.sub(&delta, NE);
    assert!(sx.is_ok(), "1 - RN150(pi)*2^-202 must be exact at p400");
    let (r, st) = x.asin_round(400, NE).unwrap();
    assert_bit_exact(
        "asin(1 - RN150(pi)*2^-202)",
        &r,
        "1.5707963267948966192313216916387627515736957026686211910464475327745659276006879479084314382208390440733644072398300969833046812798815176663",
        400,
        NE,
    );
    assert!(st.inexact());
    // And the negative side mirrors through the sign flip.
    let (rn, _) = x.negated().asin_round(400, NE).unwrap();
    assert_eq!(rn.total_cmp(&r.negated()), Ordering::Equal);
}

/// Control (the review's named input): a single-bit delta is saved
/// by sparseness and was never broken; it must stay correct.
#[test]
fn asin_near_one_sparse_delta_control() {
    let one = BigFloat::try_from_i64_exact(1, 400).unwrap();
    let (t, _) = scaled(1, 400, -200).round_to_precision(400, NE).unwrap();
    let (x, sx) = one.sub(&t, NE);
    assert!(sx.is_ok());
    let (r, _) = x.asin_round(400, NE).unwrap();
    assert_bit_exact(
        "asin(1-2^-200)",
        &r,
        "1.5707963267948966192313216916386358243075952280870463612137523814655773123323791142996378956929584997632418679662861247156808759054053569442",
        400,
        NE,
    );
}

/// pf-wmv7 (lgamma at the root x = 2): lgamma(2 + 2^-100) at
/// p120 -> 53 returned a value with relative error 2.5e-3.
#[test]
fn lgamma_positive_root_at_two_is_boosted() {
    let two = BigFloat::try_from_i64_exact(2, 120).unwrap();
    let (t, _) = scaled(1, 120, -100).round_to_precision(120, NE).unwrap();
    let (x, sx) = two.add(&t, NE);
    assert!(sx.is_ok());
    let (r, st) = x.lgamma_round(53, NE).unwrap();
    assert_bit_exact(
        "lgamma(2+2^-100)",
        &r,
        "3.33518033299040380894617486649025130707532963e-31",
        53,
        NE,
    );
    assert!(st.inexact());
}

/// pf-wmv7 (lgamma at the root x = 1): lgamma(1 + 2^-90) at
/// p100 -> 53, same mechanism at the other positive root.
#[test]
fn lgamma_positive_root_at_one_is_boosted() {
    let one = BigFloat::try_from_i64_exact(1, 100).unwrap();
    let (t, _) = scaled(1, 100, -90).round_to_precision(100, NE).unwrap();
    let (x, sx) = one.add(&t, NE);
    assert!(sx.is_ok());
    let (r, _) = x.lgamma_round(53, NE).unwrap();
    assert_bit_exact(
        "lgamma(1+2^-90)",
        &r,
        "-4.66271100848098738705521743984492072391485789e-28",
        53,
        NE,
    );
}

/// pf-rlrb (lgamma's Spouge non-window path): for `working_prec >
/// STIRLING_REACH_THRESHOLD = 600` with `x` outside the positive
/// root windows, `lgamma_at_w` fell through to
/// `lgamma_positive_at_w(..).0`, which dispatches to
/// `spouge_lgamma_scaled` and DISCARDED the returned operand scale —
/// the Spouge S-sum's internal alternating cancellation (~0.1·working
/// at z ≈ 2.5). The Ziv half-width `|y|·2^-(working − guard)` was
/// then violated by ~0.1·working bits and a wrong VALUE certified:
/// `lgamma(5/2)@1024` carried only ~73 accurate bits pre-fix (rel
/// ~1.4e-22). The input 5/2 is exactly representable, so it carries
/// no proximity depth — the cancellation is purely internal to the
/// Spouge sum, the mechanism the `spouge_lgamma_scaled` docstring
/// says the caller MUST re-drive through `cancellation_boosted`. The
/// digamma sibling (pf-0r1l, ADR-0110) routes the whole Spouge regime
/// through `cancellation_boosted`; lgamma now does too. The
/// `differential_gamma` lane caps at p256 (working ≤ 320 < 600), so
/// only these high-target rows exercise the Spouge path. References:
/// mpmath 1.4.1 @6000 bits at the exact input 5/2.
#[test]
fn lgamma_high_precision_spouge_path_is_boosted() {
    let (x, sx) = BigFloat::parse_str("2.5", 1024, NE).unwrap();
    assert!(sx.is_ok(), "5/2 is exactly representable");
    const LGAMMA_5_2: &str = "0.2846828704729191596324946696827019243201376955598947292501458503867759342216325755537007359586395675549719731391654888545210665183759870006270009069022688993799045846598114482500020082760125818560974299014792940076733089815139218871348731373151368550000610073883851700695759798293745014574546300382974729188985063926850201707466788942995896750724437363955678851450165762272223341797710883195233662922241861538669886526106809705331913486076194482374424865365347274797942082567290252635013323104413157405368732329470075198720164365473094559324340985495166190893855031207388878088248615189751211402258880348932758766301836841377252204809766770961217276363882243802785949496796610068685421040097330148704";
    for mode in [NE, RoundingMode::TowardZero, RoundingMode::TowardPositive] {
        let (r, st) = x.lgamma_round(1024, mode).unwrap();
        assert_bit_exact("lgamma(5/2)@1024", &r, LGAMMA_5_2, 1024, mode);
        assert!(st.inexact());
    }
    // Control on the Stirling path (target 256 → working ≤ 320 < 600):
    // untouched by the fix, correct before and after.
    let (x256, _) = BigFloat::parse_str("2.5", 256, NE).unwrap();
    let (r256, _) = x256.lgamma_round(256, NE).unwrap();
    assert_bit_exact("lgamma(5/2)@256", &r256, LGAMMA_5_2, 256, NE);
}

/// pf-rlrb blast radius: gamma composes `exp(lgamma(x))` and beta
/// composes `lgamma(a) + lgamma(b) − lgamma(a+b)`, each calling
/// `lgamma_round(w, …)` at its own Ziv working precision `w`. Past
/// the Spouge threshold those inner calls were the lying kernel, so
/// gamma/beta certified wrong values transitively. gamma(5/2) =
/// 3√π/4 and beta(5/2, 3/2) = π/16 at target 1024, mpmath 1.4.1
/// @6000 bits.
#[test]
fn gamma_beta_high_precision_inherit_the_lgamma_boost() {
    let (x52, sx) = BigFloat::parse_str("2.5", 1024, NE).unwrap();
    assert!(sx.is_ok());
    let (g, sg) = x52.gamma_round(1024, NE).unwrap();
    assert_bit_exact(
        "gamma(5/2)@1024",
        &g,
        "1.329340388179137020473625612505858887098162092091790346160355842389683463443274136031212992553908499062170117718211927999677114649293316951893820282202090301346528273989828842137443879771713119671699071534450972100130979261513609790387525142638925513939085230871184480235441331644429662304064499375679798805710300108106365075250992342024388877306596588373871",
        1024,
        NE,
    );
    assert!(sg.inexact());
    // beta routes a + b = 4 (a positive integer) through the same
    // lgamma compositions at high working precision.
    let (x32, _) = BigFloat::parse_str("1.5", 1024, NE).unwrap();
    let (b, sb) = x52.beta_round(&x32, 1024, NE).unwrap();
    assert_bit_exact(
        "beta(5/2,3/2)@1024",
        &b,
        "0.196349540849362077403915211454968930262323087460944113810934037019238525392888062414252176583882316748884255407080144165443365288096911389482834963008030069840642756418871157569097477788934309683148872776800685979120840383029728014611673947829450119321603035432716271788153395415513337100453765571329607786687912894724260930095057560176828380732210272993287",
        1024,
        NE,
    );
    assert!(sb.inexact());
}

/// pf-wmv7 (digamma at its positive root 1.46163...): the p100
/// parse of the root's 45-digit decimal sits ~2^-100.6 from the
/// root, a ~103-bit cancellation against the O(1) composition
/// terms — past what the first Ziv iteration's error model absorbs
/// at target 53, so the broken kernel certified garbage.
#[test]
fn digamma_positive_root_is_boosted() {
    let (x, _) =
        BigFloat::parse_str("1.4616321449683623412626595423257213284681962", 100, NE).unwrap();
    let (r, st) = x.digamma_round(53, NE).unwrap();
    assert_bit_exact(
        "digamma(RN100(root))",
        &r,
        "-4.91827109164659661351319016418437012176874373e-31",
        53,
        NE,
    );
    assert!(st.inexact());
}

/// Shallow control at the same root: a 53-bit dyadic ~2^-53.3 away
/// (psi = -9.24e-17). The ~57-bit cancellation stays inside the
/// first iteration's headroom at target 53, so this was correct
/// before the fix and must stay correct.
#[test]
fn digamma_positive_root_shallow_control() {
    let (x, sx) = BigFloat::try_from_i64_exact(6582605983432255, 53)
        .unwrap()
        .scale_by_pow2(-52);
    assert!(sx.is_ok());
    let (r, _) = x.digamma_round(53, NE).unwrap();
    assert_bit_exact(
        "digamma(6582605983432255*2^-52)",
        &r,
        "-9.24126552172942751679235141515988768650772057e-17",
        53,
        NE,
    );
}

/// pf-2thy (the digamma Spouge probe, resolved by pf-0r1l / ADR-0110
/// and locked in here): digamma's positive branch now dispatches to
/// `spouge_digamma_scaled` past `working_prec > 600` and re-drives
/// through `cancellation_boosted`, so a high-precision evaluation off
/// the root windows carries full accuracy (the probe's "920-bit
/// ceiling + 2^28-iteration shift" is gone — Spouge's cost is linear
/// in `a ∝ working` with no shift loop). digamma(5/2) at target 1024,
/// mpmath 1.4.1 @6000 bits at the exact input 5/2.
#[test]
fn digamma_high_precision_spouge_path_is_boosted() {
    let (x, sx) = BigFloat::parse_str("2.5", 1024, NE).unwrap();
    assert!(sx.is_ok());
    let (r, st) = x.digamma_round(1024, NE).unwrap();
    assert_bit_exact(
        "digamma(5/2)@1024",
        &r,
        "0.703156640645243187225690333667911099473507062006232559619539412795011695949612564517992949382082542068032257375718212718335122180663862716377151165751591206551711191266986458548378646626288061477325050427202466808022754863088151740191007734417696879630107020674946735984012026572311494321432437157666747120129795446926129310661956727232590773985887649806139170",
        1024,
        NE,
    );
    assert!(st.inexact());
}

/// Residual pf-smcb family found by the slice's adversarial
/// verification: when the input's proximity-to-1 depth exceeds the
/// working precision, `round_to_precision(x, w)` collapses `x_w` to
/// exactly 1, the series returns exact 0, and `half_width(0) = 0`
/// certifies +0 on the first iteration. Both sides of 1.
#[test]
fn ln_deep_proximity_to_one_does_not_collapse() {
    let one = BigFloat::try_from_i64_exact(1, 400).unwrap();
    let t = scaled(1, 400, -200);
    let (below, sb) = one.sub(&t, NE);
    assert!(sb.is_ok());
    let (r_below, st_b) = below.ln_round(53, NE).unwrap();
    assert_bit_exact(
        "ln(1-2^-200 @p400 -> 53)",
        &r_below,
        "-6.22301527786114170714406405378012424059025217e-61",
        53,
        NE,
    );
    assert!(st_b.inexact());
    let (above, sa) = one.add(&t, NE);
    assert!(sa.is_ok());
    let (r_above, _) = above.ln_round(53, NE).unwrap();
    assert_bit_exact(
        "ln(1+2^-200 @p400 -> 53)",
        &r_above,
        "6.22301527786114170714406405378012424059025217e-61",
        53,
        NE,
    );
}

/// Residual pf-wmv7 family found by the slice's adversarial
/// verification: the Spouge sum's internal alternating cancellation
/// (~0.1·w bits, hidden behind the ln) was not charged into the
/// reported operand scale, so the boost stopped short and
/// lgamma(2 + 2^-500 @p520) -> 400 certified a value ~2^21 ulps
/// wrong. Reference: mpmath 1.4.1 @4000 bits, 140 digits.
#[test]
fn lgamma_deep_root_charges_spouge_sum_cancellation() {
    let two = BigFloat::try_from_i64_exact(2, 520).unwrap();
    let (t, _) = scaled(1, 520, -500).round_to_precision(520, NE).unwrap();
    let (x, sx) = two.add(&t, NE);
    assert!(sx.is_ok());
    let (r, st) = x.lgamma_round(400, NE).unwrap();
    assert_bit_exact(
        "lgamma(2+2^-500 @p520 -> 400)",
        &r,
        "1.2915792392103094830071831684709114814217799294928330931850290029682202696010604343104894315468756662115269681485459002743903561637456534562e-151",
        400,
        NE,
    );
    assert!(st.inexact());
}

// ---------------------------------------------------------------
// pf-gg96 / pf-k68i / pf-pdda: input-structure-aware dispatch. The
// input encodes proximity (to the zeta pole, to a multiple of
// pi/2, to a gamma pole) beyond the working resolution; the kernel
// must resolve it exactly or grow resolution to the input's
// precision, and half_width(non-Normal) = 0 must never certify a
// collapsed special.
// ---------------------------------------------------------------

/// pf-gg96: zeta(1 + 2^-5000) at p5001 -> 53. The conditioning
/// probe rounded s to target+8 bits, collapsing s - 1 to zero; the
/// working round of s then made 1 - 2^(1-s) exactly 0, the
/// `DIV_BY_ZERO` from `eta/0` was discarded, and `half_width(inf)` = 0
/// certified +Inf with Status OK. The truth is ~2^5000 + gamma
/// (mpmath 1.4.1 @5400 bits), which rounds at 53 bits to exactly
/// 2^5000 (the Euler-gamma correction sits ~4946 bits below the
/// ulp), INEXACT.
#[test]
fn zeta_near_one_resolves_input_encoded_proximity() {
    let one = BigFloat::try_from_i64_exact(1, 5001).unwrap();
    let (t, _) = scaled(1, 5001, -5000).round_to_precision(5001, NE).unwrap();
    let (s, ss) = one.add(&t, NE);
    assert!(ss.is_ok(), "1 + 2^-5000 must be exact at p5001");
    let (r, st) = s.zeta_round(53, NE).unwrap();
    assert!(!r.is_infinite(), "zeta collapsed to +Inf, got {r}");
    let expected = scaled(1, 53, 5000);
    assert_eq!(
        r.total_cmp(&expected),
        Ordering::Equal,
        "zeta(1+2^-5000) must be 2^5000 at 53 bits, got {r}"
    );
    assert!(st.inexact(), "INEXACT missing (OK was the defect)");
    assert!(!st.div_by_zero(), "no pole was hit");
}

/// pf-k68i: `sin(RN(pi, 2048))` -> 53. `reduce()`'s `mul_prec` clamped to
/// [2048, 4096]; y = x*(2/pi) rounded to exactly 2.0, the residual
/// subtracted to exact zero, and `half_width(0)` = 0 certified -0
/// with Status OK. True value at the exact 2048-bit input: mpmath
/// 1.4.1 @4000 bits.
#[test]
fn sin_near_pi_resolves_the_reduction_residual() {
    let (x, st_pi) = pfloat::constants::pi(2048, NE);
    assert!(st_pi.inexact(), "pi at 2048 bits is inexact");
    let (r, st) = x.sin_round(53, NE).unwrap();
    assert!(!r.is_zero(), "sin collapsed to a signed zero");
    assert_bit_exact(
        "sin(RN2048(pi))",
        &r,
        "-9.19480169066345569190858554773515892681967123e-618",
        53,
        NE,
    );
    assert!(st.inexact());
}

/// pf-pdda: beta(0.5 + 2^-60 @p61, -3.5 @p3). The ADR-0030 case-5
/// dispatch classified the pole on a+b rounded to max(operand
/// precisions) = 61 bits, where -3 + 2^-60 ties-evens to exactly
/// -3: it returned +0 with Status OK where the truth is NEGATIVE
/// (-2.4913e-18, mpmath 1.4.1 @5400 bits at the exact inputs).
#[test]
fn beta_pole_dispatch_classifies_on_the_exact_sum() {
    let half = scaled(1, 61, -1);
    let (t, _) = scaled(1, 61, -60).round_to_precision(61, NE).unwrap();
    let (a, sa) = half.add(&t, NE);
    assert!(sa.is_ok(), "0.5 + 2^-60 must be exact at p61");
    let (b, sb) = BigFloat::try_from_i64_exact(-7, 3)
        .unwrap()
        .scale_by_pow2(-1);
    assert!(sb.is_ok(), "-3.5 must be exact at p3");
    let (r, st) = a.beta(&b, NE);
    assert!(
        r.is_sign_negative() && !r.is_zero(),
        "beta must be a negative normal, got {r} (status {st:?})"
    );
    assert_bit_exact(
        "beta(0.5+2^-60, -3.5)",
        &r,
        "-2.49133464143473706409934468754659056398247682e-18",
        61,
        NE,
    );
    assert!(st.inexact());
}

/// The missing fourth member of the ADR-0098 family, found by the
/// slice's adversarial verification (pre-existing): the lgamma
/// reflection collapses input-encoded proximity to the
/// negative-axis POLES inside `pi*x` before sin sees it, and the
/// realised-cancellation probe never fires because the result is
/// O(depth), not near zero. lgamma(-3+2^-80 @p84) -> 53 certified a
/// value wrong from bit ~41. mpmath 1.4.1 @4000 bits.
#[test]
fn lgamma_reflection_resolves_pole_proximity() {
    let three = BigFloat::try_from_i64_exact(-3, 84).unwrap();
    let (t, _) = scaled(1, 84, -80).round_to_precision(84, NE).unwrap();
    let (x, sx) = three.add(&t, NE);
    assert!(sx.is_ok());
    let (r, st) = x.lgamma_round(53, NE).unwrap();
    assert_bit_exact(
        "lgamma(-3+2^-80)",
        &r,
        "53.6600149755675697525660933973096055854146489",
        53,
        NE,
    );
    assert!(st.inexact());
    // digamma's reflection has the same pole structure.
    let (rd, _) = x.digamma_round(53, NE).unwrap();
    assert_bit_exact(
        "digamma(-3+2^-80)",
        &rd,
        "-1208925819614629174706174.74388233156819952727",
        53,
        NE,
    );
}

/// The deep-beta consumer of the same defect: the exact-sum handoff
/// (this slice) routes B(0.5+2^-120, -3.5)'s near-pole sum into the
/// lgamma reflection, which returned the SIGN of Gamma flipped
/// pre-fix (+2.16e-36 where the truth is negative).
#[test]
fn beta_deep_near_pole_sum_keeps_the_sign() {
    let half = scaled(1, 121, -1);
    let (t, _) = scaled(1, 121, -120).round_to_precision(121, NE).unwrap();
    let (a, sa) = half.add(&t, NE);
    assert!(sa.is_ok());
    let (b, sb) = BigFloat::try_from_i64_exact(-7, 3)
        .unwrap()
        .scale_by_pow2(-1);
    assert!(sb.is_ok());
    let (r, st) = a.beta(&b, NE);
    assert!(
        r.is_sign_negative(),
        "beta sign must be negative, got {r} ({st:?})"
    );
    assert_bit_exact(
        "beta(0.5+2^-120, -3.5)",
        &r,
        "-2.160888344505549714961133716867520191111181e-36",
        121,
        NE,
    );
}

// ---------------------------------------------------------------
// pf-qm0h / pf-6nn5 / pf-lkno hotfix (ADR-0101): the exp-composing
// kernels discarded exp's rim Status and certified bare specials —
// exp2/exp10 with INEXACT only, expm1/sinh/cosh with Status::OK on
// transcendental results (the Class::Normal INEXACT force skips
// non-Normals). Each now forwards through exp's ADR-0096 rim
// machinery, mode-aware values and flags included. Found by the
// R1-merge CI failure in pfloat-libm's consistency lane.
// ---------------------------------------------------------------

/// exp2 past the rim: 2^(2^90) exceeds every representable.
#[test]
fn exp2_rim_is_mode_aware_and_flagged() {
    let x = scaled(1, 53, 90);
    let (r_ne, st_ne) = x.exp2_round(53, NE).unwrap();
    assert!(r_ne.is_infinite(), "NE must be +inf, got {r_ne}");
    assert!(st_ne.overflow() && st_ne.inexact(), "got {st_ne:?}");
    let (r_tz, st_tz) = x.exp2_round(53, RoundingMode::TowardZero).unwrap();
    assert_eq!(r_tz.total_cmp(&max_finite(53)), Ordering::Equal);
    assert!(st_tz.overflow() && st_tz.inexact());
    // The negative mirror underflows.
    let (r_neg, st_neg) = x.negated().exp2_round(53, NE).unwrap();
    assert!(r_neg.is_zero(), "NE must be +0, got {r_neg}");
    assert!(st_neg.underflow() && st_neg.inexact());
    let (r_tp, _) = x
        .negated()
        .exp2_round(53, RoundingMode::TowardPositive)
        .unwrap();
    assert_eq!(r_tp.total_cmp(&min_pos(53)), Ordering::Equal);
}

/// exp2's representable window at the rim band: floor(x) is exact
/// integer arithmetic, so 2^(2^62 + 0.5) = sqrt(2)·2^(2^62) must
/// come back finite with INEXACT only. mpmath 1.4.1: sqrt(2)
/// mantissa quoted.
#[test]
fn exp2_rim_window_returns_the_representable_finite() {
    let half = scaled(1, 64, -1);
    let (x, sx) = scaled(1, 64, 62).add(&half, NE);
    assert!(sx.is_ok());
    let (r, st) = x.exp2_round(53, NE).unwrap();
    let expected = mantissa_at_pow2(
        "1.414213562373095048801688724209698078570",
        4_611_686_018_427_387_904,
        NE,
    );
    assert_eq!(r.total_cmp(&expected), Ordering::Equal, "got {r}");
    assert!(st.inexact() && !st.overflow() && !st.underflow());
}

/// exp10/expm1/cosh/sinh forward exp's rim flags instead of
/// certifying bare specials with OK.
#[test]
fn exp_composers_forward_rim_flags() {
    let (x, _) = BigFloat::parse_str("1e300", 64, NE).unwrap();
    let (r10, st10) = x.exp10_round(53, NE).unwrap();
    assert!(r10.is_infinite() && st10.overflow() && st10.inexact());
    let (rm1, stm1) = x.expm1_round(53, NE).unwrap();
    assert!(rm1.is_infinite() && stm1.overflow() && stm1.inexact());
    let (rc, stc) = x.cosh_round(53, NE).unwrap();
    assert!(
        rc.is_infinite() && stc.overflow() && stc.inexact(),
        "cosh(1e300) = ({rc}, {stc:?}); Status::OK was the defect"
    );
    let (rs, sts) = x.sinh_round(53, NE).unwrap();
    assert!(rs.is_infinite() && sts.overflow() && sts.inexact());
    // sinh's negative side mirrors the directed modes through the
    // negation: under TowardPositive the inward magnitude gives
    // -MaxFinite, not -inf.
    let (rsn, stsn) = x
        .negated()
        .sinh_round(53, RoundingMode::TowardPositive)
        .unwrap();
    assert_eq!(
        rsn.total_cmp(&max_finite(53).negated()),
        Ordering::Equal,
        "sinh(-1e300) TP must be -MaxFinite, got {rsn}"
    );
    assert!(stsn.overflow() && stsn.inexact());
    // cosh under TowardZero takes the inward finite.
    let (rcz, _) = x.cosh_round(53, RoundingMode::TowardZero).unwrap();
    assert_eq!(rcz.total_cmp(&max_finite(53)), Ordering::Equal);
}

/// ADR-0101 verifier round 2: the cosh/sinh rim forward initially
/// landed |x| in [2^62, 2^62+ln2) on exp's LEGACY path, whose
/// reduction cancelled ~`e_x` bits against a flat 24-bit guard —
/// `cosh(2^62+0.5)` at 53 NE certified the wrong NE side (mantissa
/// ...c800 where the triple-oracled truth is 0.0227 ulp above the
/// midpoint: ...d000). The legacy reduction now carries `e_x` extra
/// bits (closing pf-t6ht's probed band for every exp caller).
#[test]
fn cosh_band_above_rim_trigger_rounds_the_ne_side() {
    let (h, _) = BigFloat::try_from_i64_exact(1, 64)
        .unwrap()
        .scale_by_pow2(-1);
    let b = scaled(1, 64, 62);
    let (x, sx) = b.add(&h, NE);
    assert!(sx.is_ok());
    for f in [BigFloat::cosh_round, BigFloat::sinh_round] {
        let (r, st) = f(&x, 53, NE).unwrap();
        match r.parts() {
            pfloat::Parts::Normal {
                exponent, mantissa, ..
            } => {
                assert_eq!(exponent, 6_653_256_548_922_161_245);
                assert_eq!(
                    mantissa,
                    &[13_916_968_917_729_267_712_u64],
                    "wrong NE side (the ...c800 defect)"
                );
            }
            other => panic!("expected a Normal, got {other:?}"),
        }
        assert!(st.inexact() && !st.overflow());
    }
}

// ---------------------------------------------------------------
// pf-71u2 / pf-e2ow (arc R2, slice R2.1, ADR-0102): the scalar
// Ziv-cap lies on input-encoded depth. hypot's eval absorbs the
// small operand's square whenever 2·gap ≥ working, collapses onto
// |big| (on-grid), the interval test never converges, and the
// exhausted fall-through certified the collapsed value — falsely
// EXACT whenever |big| rounds exactly at the target. atan (and
// atan2 through it) has the same shape at depth 2·|e|: the cubic
// correction sits below every working precision the driver
// reaches, so directed modes returned the argument itself. The
// depth is input-encoded and exactly computable up front, so both
// kernels now dispatch through round_with_infinitesimal past the
// representable band (truth strictly inside the boundary-free gap
// above/below the base). Directions pinned with mpmath 1.4.1:
// sqrt(1+2^-4000) − 1 ≈ 2^-4001 ∈ (0, 2^-127); atan(ε) < ε with
// ε − atan(ε) ≈ ε³/3 (both signs).
// ---------------------------------------------------------------

/// pf-71u2, the named reproducer: hypot(1, 2^-2000) at target 128.
/// The broken kernel returned (exactly 1, `Status::OK`) — wrong
/// value under `TowardPositive` and falsely exact in every mode.
#[test]
fn hypot_deep_gap_is_mode_aware_and_inexact() {
    let one = BigFloat::try_from_i64_exact(1, 128).unwrap();
    let tiny = scaled(1, 128, -2000);
    for mode in ALL_MODES {
        let (r, st) = one.hypot_round(&tiny, 128, mode).unwrap();
        let expected = if matches!(mode, RoundingMode::TowardPositive) {
            one.next_up().0
        } else {
            one.clone()
        };
        assert_eq!(
            r.total_cmp(&expected),
            Ordering::Equal,
            "hypot(1, 2^-2000) {mode:?}: got {r}, want {expected}"
        );
        assert!(
            st.inexact(),
            "hypot(1, 2^-2000) {mode:?}: INEXACT missing (OK was the defect)"
        );
    }
    // Symmetric in operand order and sign-independent.
    let (r_sw, st_sw) = tiny
        .negated()
        .hypot_round(&one.negated(), 128, RoundingMode::TowardPositive)
        .unwrap();
    assert_eq!(r_sw.total_cmp(&one.next_up().0), Ordering::Equal);
    assert!(st_sw.inexact());
}

/// Same defect at a high target: the gap (2·1600) clears the old
/// driver cap (2000 + 1024) and the representable band both.
#[test]
fn hypot_deep_gap_high_target() {
    let one = BigFloat::try_from_i64_exact(1, 64).unwrap();
    let tiny = scaled(1, 64, -1600);
    let (r_tp, st_tp) = one
        .hypot_round(&tiny, 2000, RoundingMode::TowardPositive)
        .unwrap();
    let one_2000 = BigFloat::try_from_i64_exact(1, 2000).unwrap();
    assert_eq!(
        r_tp.total_cmp(&one_2000.next_up().0),
        Ordering::Equal,
        "TP must be nextUp(1@2000), got {r_tp}"
    );
    assert!(st_tp.inexact());
    let (r_ne, st_ne) = one.hypot_round(&tiny, 2000, NE).unwrap();
    assert_eq!(r_ne.total_cmp(&one_2000), Ordering::Equal);
    assert!(st_ne.inexact());
}

/// Control: a deep-ish-gap exact Pythagorean triple stays exact
/// with Status OK. (2^62 − 1)² + (2^32)² = (2^62 + 1)²; the gap
/// (29) sits inside the Ziv-certifiable band, which is untouched.
#[test]
fn hypot_pythagorean_band_control() {
    let a = BigFloat::try_from_i64_exact((1_i64 << 62) - 1, 63).unwrap();
    let b = scaled(1, 63, 32);
    let (r, st) = a.hypot_round(&b, 63, NE).unwrap();
    let expected = BigFloat::try_from_i64_exact((1_i64 << 62) + 1, 63).unwrap();
    assert_eq!(r.total_cmp(&expected), Ordering::Equal, "got {r}");
    assert!(st.is_ok(), "exact hypot must stay OK, got {st:?}");
}

/// The root of pf-e2ow lives in atan itself (atan2 composes it):
/// atan(2^-545) at target 64 under the inward modes must be
/// pred(2^-545), not the argument (truth = ε − ε³/3 + … < ε).
#[test]
fn atan_tiny_x_directed_modes_round_inward() {
    let x = scaled(1, 64, -545);
    let pred = x.next_down().0;
    for (mode, expected) in [
        (RoundingMode::TowardZero, &pred),
        (RoundingMode::TowardNegative, &pred),
        (NE, &x),
        (RoundingMode::NearestAway, &x),
        (RoundingMode::TowardPositive, &x),
    ] {
        let (r, st) = x.atan_round(64, mode).unwrap();
        assert_eq!(
            r.total_cmp(expected),
            Ordering::Equal,
            "atan(2^-545) {mode:?}: got {r}, want {expected}"
        );
        assert!(st.inexact(), "atan(2^-545) {mode:?}: INEXACT missing");
    }
    // Negative mirror: |atan(x)| < |x|, so the magnitude-shrinking
    // modes are TowardZero and TowardPositive.
    let xn = x.negated();
    let predn = pred.negated();
    let (r_tz, _) = xn.atan_round(64, RoundingMode::TowardZero).unwrap();
    assert_eq!(r_tz.total_cmp(&predn), Ordering::Equal);
    let (r_tp, _) = xn.atan_round(64, RoundingMode::TowardPositive).unwrap();
    assert_eq!(r_tp.total_cmp(&predn), Ordering::Equal);
    let (r_tn, _) = xn.atan_round(64, RoundingMode::TowardNegative).unwrap();
    assert_eq!(r_tn.total_cmp(&xn), Ordering::Equal);
}

/// Depth past every feasible working precision (2·2^40 bits): only
/// the infinitesimal dispatch can resolve this; a precision boost
/// cannot represent it.
#[test]
fn atan_tiny_x_past_every_feasible_precision() {
    let x = scaled(1, 64, -(1_i64 << 40));
    let (r_tz, st) = x.atan_round(64, RoundingMode::TowardZero).unwrap();
    assert_eq!(
        r_tz.total_cmp(&x.next_down().0),
        Ordering::Equal,
        "TZ must be pred, got {r_tz}"
    );
    assert!(st.inexact());
    let (r_ne, _) = x.atan_round(64, NE).unwrap();
    assert_eq!(r_ne.total_cmp(&x), Ordering::Equal);
}

/// pf-e2ow, the named reproducer: atan2(2^-545, 1) at target 64
/// under the inward modes returned the argument itself.
#[test]
fn atan2_tiny_exact_ratio_directed_modes_round_inward() {
    let y = scaled(1, 64, -545);
    let x = BigFloat::try_from_i64_exact(1, 64).unwrap();
    let pred = y.next_down().0;
    for (mode, expected) in [
        (RoundingMode::TowardZero, &pred),
        (RoundingMode::TowardNegative, &pred),
        (NE, &y),
        (RoundingMode::TowardPositive, &y),
    ] {
        let (r, st) = y.atan2_round(&x, 64, mode).unwrap();
        assert_eq!(
            r.total_cmp(expected),
            Ordering::Equal,
            "atan2(2^-545, 1) {mode:?}: got {r}, want {expected}"
        );
        assert!(st.inexact());
    }
    // Negative y mirrors: magnitude shrinks toward zero.
    let yn = y.negated();
    let (r_tp, _) = yn
        .atan2_round(&x, 64, RoundingMode::TowardPositive)
        .unwrap();
    assert_eq!(r_tp.total_cmp(&pred.negated()), Ordering::Equal);
    let (r_tn, _) = yn
        .atan2_round(&x, 64, RoundingMode::TowardNegative)
        .unwrap();
    assert_eq!(r_tn.total_cmp(&yn), Ordering::Equal);
    // And an exact ratio deeper than any feasible precision.
    let yd = scaled(1, 64, -(1_i64 << 40));
    let (r_deep, _) = yd.atan2_round(&x, 64, RoundingMode::TowardZero).unwrap();
    assert_eq!(r_deep.total_cmp(&yd.next_down().0), Ordering::Equal);
}

/// Control: an inexact tiny ratio (2^-600 / 3). The truth's
/// distance to the target grid is carried by the ratio's own
/// expansion (generic position), so the driver's fall-through
/// already rounds it correctly; it must stay correct. mpmath 1.4.1
/// @1000 bits at the exact dyadic inputs.
#[test]
fn atan2_tiny_inexact_ratio_control() {
    let y = scaled(1, 64, -600);
    let x = BigFloat::try_from_i64_exact(3, 64).unwrap();
    for mode in [NE, RoundingMode::TowardZero] {
        let (r, st) = y.atan2_round(&x, 64, mode).unwrap();
        assert_bit_exact(
            "atan2(2^-600, 3)",
            &r,
            "8.0330662170096137258025001157083631214366831816997e-182",
            64,
            mode,
        );
        assert!(st.inexact());
    }
}

// ---------------------------------------------------------------
// pf-jl35 (arc R2, slice R2.2, ADR-0103): the Ziv driver's fixed
// cap (target + 1024) silently fell back on input-encoded
// proximity deeper than the cap — deep directed modes 1 ulp wrong
// with plain INEXACT. The driver now takes a lazily-evaluated
// input-derived certification depth: on legacy-schedule
// exhaustion, one further iteration at target + depth + 64
// certifies under the identical interval test. zeta derives the
// depth from its ADR-0098 conditioning probe; the trig family
// from the reduction residual's exponent. All references mpmath
// 1.4.1 at 2048-9000 bits at the bit-identical parsed inputs.
// The bead's claimed correct value for cos (−1 + 2^-53) was
// itself misadjudicated: the representable just above −1 at p53
// is −(1 − 2^-54) (the binade below 1 has ulp 2^-54).
// ---------------------------------------------------------------

/// pf-jl35, the named zeta reproducer: ζ(1 − 2^-2000 @p2001) → 53
/// under TowardPositive. Truth = −2^2000 + γ + O(2^-2000)
/// (mpmath: ζ + 2^2000 = 0.57721…), strictly inside
/// (−2^2000, −(2^2000 − 2^1947)), so TP must round up to
/// −(2^2000 − 2^1947); the exhausted driver returned −2^2000
/// (1 ulp low, INEXACT, run-verified red at e8b1284). NE stays
/// −2^2000 (γ is ~4946 bits below the half-ulp).
///
/// Release builds only: the conditioning-deep Borwein evaluations
/// cost ~1 min in release but ~15 min in the debug matrix — the
/// row runs in the MPFR full-union release job (whose feature set
/// covers this lane). Any red zeta instance needs its depth past
/// the legacy cap plus the eval's carry slack (shallower depths
/// survive by uncertified fallback: the eval still carries γ), so
/// no cheap debug instance of this kernel's wiring exists; the
/// trig rows below guard the shared driver mechanism in debug.
#[cfg(not(debug_assertions))]
#[test]
fn zeta_near_pole_deep_directed_certifies() {
    let one = BigFloat::try_from_i64_exact(1, 2001).unwrap();
    let (t, _) = scaled(1, 2001, -2000).round_to_precision(2001, NE).unwrap();
    let (s, ss) = one.sub(&t, NE);
    assert!(ss.is_ok(), "1 - 2^-2000 must be exact at p2001");
    let neg_big = scaled(1, 53, 2000).negated();
    let (r_tp, st_tp) = s.zeta_round(53, RoundingMode::TowardPositive).unwrap();
    assert_eq!(
        r_tp.total_cmp(&neg_big.next_up().0),
        Ordering::Equal,
        "zeta TP must be -(2^2000 - 2^1947), got {r_tp}"
    );
    assert!(st_tp.inexact());
    let (r_ne, st_ne) = s.zeta_round(53, NE).unwrap();
    assert_eq!(r_ne.total_cmp(&neg_big), Ordering::Equal, "NE control");
    assert!(st_ne.inexact());
}

/// pf-jl35, the named trig reproducer: cos(RN2048(π)) → 53 under
/// the inward modes. Truth = −1 + 4.227e-1235 (mpmath @9000),
/// strictly inside (−1, −(1 − 2^-54)): `TowardPositive` and
/// `TowardZero` must give `nextUp(−1)` = −(1 − 2^-54); the exhausted
/// driver returned −1. The output needs ~4103 bits, past the old
/// cap; the depth comes from the reduction residual.
#[test]
fn cos_near_pi_deep_directed_certifies() {
    let (x, _) = pfloat::constants::pi(2048, NE);
    let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
    let expected = neg_one.next_up().0;
    for mode in [RoundingMode::TowardPositive, RoundingMode::TowardZero] {
        let (r, st) = x.cos_round(53, mode).unwrap();
        assert_eq!(
            r.total_cmp(&expected),
            Ordering::Equal,
            "cos {mode:?} must be nextUp(-1), got {r}"
        );
        assert!(st.inexact());
    }
    let (r_ne, _) = x.cos_round(53, NE).unwrap();
    assert_eq!(r_ne.total_cmp(&neg_one), Ordering::Equal, "NE control");
}

/// sin's exposed arm is the cos-shape quadrant: x = RN2048(π)/2
/// (exact halving) sits ~2^-2050 from π/2, sin(x) = 1 − 1.06e-1235
/// (mpmath @9000), strictly inside (1 − 2^-54, 1): the inward modes
/// must give pred(1) = 1 − 2^-54, the nearest modes 1.
#[test]
fn sin_near_half_pi_deep_directed_certifies() {
    let (pi, _) = pfloat::constants::pi(2048, NE);
    let (x, sx) = pi.scale_by_pow2(-1);
    assert!(sx.is_ok());
    let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
    let pred = one.next_down().0;
    for mode in [RoundingMode::TowardZero, RoundingMode::TowardNegative] {
        let (r, st) = x.sin_round(53, mode).unwrap();
        assert_eq!(
            r.total_cmp(&pred),
            Ordering::Equal,
            "sin {mode:?} must be pred(1), got {r}"
        );
        assert!(st.inexact());
    }
    let (r_ne, _) = x.sin_round(53, NE).unwrap();
    assert_eq!(r_ne.total_cmp(&one), Ordering::Equal);
}

/// The reciprocal family inherits through the shared helper:
/// sec(RN2048(π)) = −1 − 4.227e-1235 (mpmath @9000), strictly
/// inside (−(1 + 2^-52), −1): `TowardNegative` must give
/// −(1 + 2^-52), the inward modes −1.
#[test]
fn sec_near_pi_deep_directed_certifies() {
    let (x, _) = pfloat::constants::pi(2048, NE);
    let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
    let (r_tn, st_tn) = x.sec_round(53, RoundingMode::TowardNegative).unwrap();
    assert_eq!(
        r_tn.total_cmp(&neg_one.next_down().0),
        Ordering::Equal,
        "sec TN must be -(1 + 2^-52), got {r_tn}"
    );
    assert!(st_tn.inexact());
    let (r_tz, _) = x.sec_round(53, RoundingMode::TowardZero).unwrap();
    assert_eq!(r_tz.total_cmp(&neg_one), Ordering::Equal);
}

/// The reduction range cap is now on the input's exponent alone
/// (`e_x < 4032`); the old working-coupled form (`e_x + working +
/// 64 < 4096`) refused any deep working precision outright — which
/// the deep rung legitimately requests, and which also spuriously
/// refused high TARGETS for mid-range exponents. sin(2^3000) at
/// target 1000 was NaN + INVALID before; it must now compute, and
/// agree with its own higher-precision rounding (the refinement
/// self-consistency oracle).
#[test]
fn sin_high_target_mid_exponent_no_longer_refused() {
    let x = scaled(1, 53, 3000);
    let (r1000, st) = x.sin_round(1000, NE).unwrap();
    assert!(
        !r1000.is_nan() && !st.invalid(),
        "sin(2^3000) at target 1000 must compute, got {r1000} ({st:?})"
    );
    assert!(st.inexact());
    let (r1500, _) = x.sin_round(1500, NE).unwrap();
    let (r1500_down, _) = r1500.round_to_precision(1000, NE).unwrap();
    assert_eq!(
        r1000.total_cmp(&r1500_down),
        Ordering::Equal,
        "refinement self-consistency"
    );
}

// ---------------------------------------------------------------
// pf-fbjn (arc R2, slice R2.2b, ADR-0104): the ADR-0059 tiny-x
// fast-path family lacked the input-precision arm — a
// high-precision input parked next to a rounding-change point puts
// the series correction across the boundary while the
// round_with_infinitesimal residue stays on the input's side, and
// the fast path rounds the wrong side. The arm alone cannot fix it
// (the Ziv fall-through collapses the input onto the same midpoint
// below the input precision); the repair is arm + ADR-0103 depth
// hint, so the driver's deep rung takes the input at full
// precision and certifies the true boundary side. The same hints
// close the atan/atan2/hypot band residues ADR-0102 recorded. All
// directions pinned with mpmath 1.4.1 at 8000-12000 bits.
// ---------------------------------------------------------------

/// Builds `2^-600 · (1 + 2^-53 + tail·2^-1300-ish)` exactly: the
/// p53 NE midpoint above 2^-600, perturbed by a deep tail on the
/// chosen side.
fn parked_at_midpoint(tail_negative: bool, tail_pow: i64, p: u32) -> BigFloat {
    let one = BigFloat::try_from_i64_exact(1, p).unwrap();
    let (a, sa) = one.add(&scaled(1, p, -53), NE);
    assert!(sa.is_ok());
    let t = scaled(1, p, tail_pow);
    let (m, sm) = if tail_negative {
        a.sub(&t, NE)
    } else {
        a.add(&t, NE)
    };
    assert!(sm.is_ok(), "the parked construction must be exact");
    let (x, sx) = m.scale_by_pow2(-600);
    assert!(sx.is_ok());
    x
}

/// The 53-bit grid points bracketing the midpoint at 2^-600.
fn bracket_at_minus_600() -> (BigFloat, BigFloat) {
    let lower = scaled(1, 53, -600);
    (lower.clone(), lower.next_up().0)
}

/// atanh, the family representative (cubic, grows): x parked
/// 2^-2600 BELOW the midpoint at p2001; the +x³/3 correction
/// (≈ 2^-1801.6) carries the truth ABOVE it (mpmath margin
/// +4.67e-543). The broken fast path's residue stayed below and
/// rounded down; the arm + deep rung must round up.
#[test]
fn atanh_parked_high_precision_rounds_the_true_side() {
    let x = parked_at_midpoint(true, -2000, 2001);
    let (lower, upper) = bracket_at_minus_600();
    for mode in [NE, RoundingMode::NearestAway, RoundingMode::TowardPositive] {
        let (r, st) = x.atanh_round(53, mode).unwrap();
        assert_eq!(
            r.total_cmp(&upper),
            Ordering::Equal,
            "atanh {mode:?}: got {r}, want the upper point (truth above the midpoint)"
        );
        assert!(st.inexact());
    }
    let (r_tz, _) = x.atanh_round(53, RoundingMode::TowardZero).unwrap();
    assert_eq!(r_tz.total_cmp(&lower), Ordering::Equal);
}

/// The shrink-direction mirror (asinh, −x³/6): x parked ABOVE the
/// midpoint, truth BELOW it (mpmath margin −4.67e-543 family).
#[test]
fn asinh_parked_high_precision_rounds_the_true_side() {
    let x = parked_at_midpoint(false, -2000, 2001);
    let (lower, _) = bracket_at_minus_600();
    for mode in [NE, RoundingMode::NearestAway, RoundingMode::TowardZero] {
        let (r, st) = x.asinh_round(53, mode).unwrap();
        assert_eq!(
            r.total_cmp(&lower),
            Ordering::Equal,
            "asinh {mode:?}: got {r}, want the lower point (truth below the midpoint)"
        );
        assert!(st.inexact());
    }
}

/// The quadratic pair. expm1 (+x²/2 ≈ 2^-1201): x parked 2^-1900
/// below the midpoint at p1301, truth above (mpmath margin
/// +2.9e-362). log1p (−x²/2): parked above, truth below.
#[test]
fn expm1_log1p_parked_high_precision_round_the_true_side() {
    let xe = parked_at_midpoint(true, -1300, 1301);
    let (lower, upper) = bracket_at_minus_600();
    let (re, ste) = xe.expm1_round(53, NE).unwrap();
    assert_eq!(
        re.total_cmp(&upper),
        Ordering::Equal,
        "expm1 NE: got {re}, want the upper point"
    );
    assert!(ste.inexact());
    let xl = parked_at_midpoint(false, -1300, 1301);
    let (rl, stl) = xl.log1p_round(53, NE).unwrap();
    assert_eq!(
        rl.total_cmp(&lower),
        Ordering::Equal,
        "log1p NE: got {rl}, want the lower point"
    );
    assert!(stl.inexact());
}

/// The arm-and-cap interaction guard (ADR-0104): an arm-failing
/// VERY deep input (p ≥ |e| − 2 with |e| past the old internal
/// boost cap) reaches the driver, whose expm1/log1p closures used
/// to cap their internal cancellation boost at +1024 — the
/// composition then collapsed to exactly 0 and `half_width(0)` = 0
/// certified it with Status OK at the first rung (the
/// review-2026-05-29 certified-zero class, which the precision arm
/// alone would have resurrected). With the arm in place the boost
/// is input-proportional and uncapped. x = (1 + 2^-2997)·2^-3000
/// at p2998: truth = x + x²/2 + … sits strictly between 2^-3000
/// and its successor at p53.
#[test]
fn expm1_log1p_arm_failing_deep_inputs_never_certify_zero() {
    let one = BigFloat::try_from_i64_exact(1, 2998).unwrap();
    let (m, sm) = one.add(&scaled(1, 2998, -2997), NE);
    assert!(sm.is_ok());
    let (x, sx) = m.scale_by_pow2(-3000);
    assert!(sx.is_ok());
    let expect = scaled(1, 53, -3000);
    for mode in [NE, RoundingMode::TowardZero] {
        let (r, st) = x.expm1_round(53, mode).unwrap();
        assert!(!r.is_zero(), "expm1 {mode:?} certified zero");
        assert_eq!(r.total_cmp(&expect), Ordering::Equal, "expm1 {mode:?}");
        assert!(st.inexact());
    }
    let (r_tp, _) = x.expm1_round(53, RoundingMode::TowardPositive).unwrap();
    assert_eq!(
        r_tp.total_cmp(&expect.next_up().0),
        Ordering::Equal,
        "expm1 TP must take the upper neighbour"
    );
    let (rl, stl) = x.log1p_round(53, NE).unwrap();
    assert!(!rl.is_zero(), "log1p certified zero");
    assert_eq!(rl.total_cmp(&expect), Ordering::Equal);
    assert!(stl.inexact());
}

/// The atan band residue ADR-0102 recorded (precision past the
/// arm, depth past the driver cap): x = 2^-600 + 2^-2199 exact at
/// p1600; truth = x − x³/3 sits BELOW 2^-600 (mpmath), so
/// `TowardZero` must give pred(2^-600) where the capped driver
/// returned 2^-600.
#[test]
fn atan_band_high_precision_resolves_through_the_deep_rung() {
    let (x, sx) = scaled(1, 1600, -600).add(&scaled(1, 1600, -2199), NE);
    assert!(sx.is_ok());
    let lower = scaled(1, 53, -600);
    let (r_tz, st) = x.atan_round(53, RoundingMode::TowardZero).unwrap();
    assert_eq!(
        r_tz.total_cmp(&lower.next_down().0),
        Ordering::Equal,
        "atan TZ must be pred(2^-600), got {r_tz}"
    );
    assert!(st.inexact());
    let (r_ne, _) = x.atan_round(53, NE).unwrap();
    assert_eq!(r_ne.total_cmp(&lower), Ordering::Equal, "NE control");
}

/// The hypot hole ADR-0102 recorded (big-operand precision past
/// the dispatch arm, gap past the driver cap): big = 1 + 2^-4990
/// at p5000, small = 2^-700. The truth 1 + ~2^-1401 (mpmath,
/// strictly below the p53 midpoint) was certified falsely EXACT as
/// (1, OK); `TowardPositive` must give `nextUp(1)`, `NearestEven` 1,
/// INEXACT always.
#[test]
fn hypot_hole_high_precision_resolves_through_the_deep_rung() {
    let one = BigFloat::try_from_i64_exact(1, 5000).unwrap();
    let (big, sb) = one.add(&scaled(1, 5000, -4990), NE);
    assert!(sb.is_ok());
    let small = scaled(1, 64, -700);
    let one53 = BigFloat::try_from_i64_exact(1, 53).unwrap();
    let (r_tp, st_tp) = big
        .hypot_round(&small, 53, RoundingMode::TowardPositive)
        .unwrap();
    assert_eq!(
        r_tp.total_cmp(&one53.next_up().0),
        Ordering::Equal,
        "hypot TP must be nextUp(1), got {r_tp}"
    );
    assert!(st_tp.inexact());
    let (r_ne, st_ne) = big.hypot_round(&small, 53, NE).unwrap();
    assert_eq!(r_ne.total_cmp(&one53), Ordering::Equal);
    assert!(
        st_ne.inexact(),
        "INEXACT missing (falsely-exact was the defect)"
    );
}

// ---------------------------------------------------------------
// pf-9761 (arc R2, slice R2.3, ADR-0105): acos near 1 certified
// +0 WITH STATUS OK. acos does not compose asin: its two-branch
// atan form computes 1 − x_w at the Ziv working precision with no
// input-derived boost, so input-encoded proximity to 1 collapses
// x_w to exactly 1, the evaluation returns exact 0,
// half_width(0) = 0 certifies on the first rung, and the
// defensive INEXACT force skipped non-Normal results — the wrong
// zero claimed exactness. Fix: the ADR-0097 asin gap-boost on
// both branches (1 − |x| Sterbenz-exact at the input precision),
// plus the class-wide posture fix: the INEXACT force now covers
// Zero results (every kernel in the family pre-dispatches its
// exact-zero inputs, so a post-driver Zero is always a collapse).
// References mpmath 1.4.1 @4000 bits at the exact dyadic inputs.
// ---------------------------------------------------------------

/// pf-9761, the named reproducer: acos(1 − RN150(π)·2^-202 @p400)
/// → 53 returned +0 with Status OK (truth ≈ 9.887e-31). The same
/// construction as the asin dense-delta row.
#[test]
fn acos_near_one_resolves_input_encoded_proximity() {
    let (pi150, sp) =
        BigFloat::parse_str("1120957716564506572603712206968581818470252692", 150, NE).unwrap();
    assert!(sp.is_ok());
    let (delta, sd) = pi150.scale_by_pow2(-148 - 202);
    assert!(sd.is_ok());
    let one = BigFloat::try_from_i64_exact(1, 400).unwrap();
    let (x, sx) = one.sub(&delta, NE);
    assert!(sx.is_ok(), "1 - RN150(pi)*2^-202 must be exact at p400");
    for mode in [NE, RoundingMode::TowardZero] {
        let (r, st) = x.acos_round(53, mode).unwrap();
        assert!(!r.is_zero(), "acos collapsed to zero ({mode:?})");
        assert_bit_exact(
            "acos(1 - RN150(pi)*2^-202)",
            &r,
            "9.8869052488899701893171944102476337934227554241655e-31",
            53,
            mode,
        );
        assert!(
            st.inexact(),
            "{mode:?}: INEXACT missing (OK was the defect)"
        );
    }
    // The negative branch carries the same collapse shape through
    // 1 + x_w; acos(−(1 − δ)) = π − acos(1 − δ).
    let (rn, stn) = x.negated().acos_round(53, NE).unwrap();
    assert_bit_exact(
        "acos(-(1 - RN150(pi)*2^-202))",
        &rn,
        "3.1415926535897932384626433832785141936722804023562",
        53,
        NE,
    );
    assert!(stn.inexact());
}

/// Shallow control (inside the first working precision's reach):
/// acos(1 − 2^-40 @p53), correct before the fix and bit-pinned.
#[test]
fn acos_near_one_shallow_control() {
    let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
    let (x, sx) = one.sub(&scaled(1, 53, -40), NE);
    assert!(sx.is_ok());
    let (r, st) = x.acos_round(53, NE).unwrap();
    assert_bit_exact(
        "acos(1-2^-40)",
        &r,
        "0.0000013486991523487112367441192231378786709130652500137",
        53,
        NE,
    );
    assert!(st.inexact());
}

// ---------------------------------------------------------------
// pf-06lk (arc R2, slice R2.4, ADR-0106): agm loop exponent
// saturation. The Gauss iteration's mul/sqrt saturate their result
// exponents at the i64 rim per the no-emax contract; the eval
// closure discarded those statuses, the clamp is w-independent, so
// the Ziv interval test certified the corrupted iterates — worst
// case agm(2^(2^62+1000), 3·2^(2^62+1000)) returned a value ~10^31
// wrong with Status OK. AGM is degree-1 homogeneous:
// agm(s·a, s·b) = s·agm(a, b), so the kernel now normalizes the
// operands toward exponent 0 (the midpoint of the operand
// exponents, an exact power-of-two scaling), runs the loop near
// scale 1 where no internal operation can reach the rim, and
// scales the rounded result back exactly. Expected values:
// AGM(2^k, 3·2^k) = AGM(1, 3)·2^k with AGM(1, 3) from mpmath 1.4.1
// @400 bits, parsed at 53 bits then scaled exactly.
// ---------------------------------------------------------------

/// The boundary case (k = 2^62, 47% low pre-fix), the worst case
/// (k = 2^62 + 1000, ~10^31 wrong WITH STATUS OK pre-fix), the
/// negative mirror, and the in-range control (k = 2^62 − 30,
/// correct pre-fix and pinned so the normalization preserves it).
#[test]
fn agm_rim_exponents_certify_the_true_agm() {
    let agm13 = "1.86361678324489654235568903102427059515753285682685372222044";
    let k = 1_i64 << 62;
    for (label, kk, expect_ok_defect) in [
        ("k=2^62", k, false),
        ("k=2^62+1000", k + 1000, true),
        ("k=-(2^62+1000)", -(k + 1000), false),
        ("k=2^62-30", k - 30, false),
    ] {
        let a = scaled(1, 53, kk);
        let b = scaled(3, 53, kk);
        let (r, st) = a.agm(&b, NE);
        let (m, _) = BigFloat::parse_str(agm13, 53, NE).unwrap();
        let (expected, se) = m.scale_by_pow2(kk);
        assert!(se.is_ok());
        assert_eq!(
            r.total_cmp(&expected),
            Ordering::Equal,
            "agm {label}: got {r}, want {expected}"
        );
        assert!(
            st.inexact(),
            "agm {label}: INEXACT missing{}",
            if expect_ok_defect {
                " (Status OK was the worst-case defect)"
            } else {
                ""
            }
        );
    }
}

/// Asserts a Normal result's top-aligned mantissa limb and exponent
/// (the cosh rim-row precedent: Display saturates at astronomical
/// exponents, so structural comparison is the only honest one).
fn assert_parts(label: &str, r: &BigFloat, mantissa: u64, exponent: i64) {
    match r.parts() {
        pfloat::Parts::Normal {
            exponent: e,
            mantissa: m,
            ..
        } => {
            assert_eq!(e, exponent, "{label}: exponent");
            assert_eq!(m, &[mantissa], "{label}: mantissa");
        }
        other => panic!("{label}: expected a Normal, got {other:?}"),
    }
}

/// ADR-0106 verifier refutation 1: the exponent-midpoint's negation
/// overflowed for both operands at the bottom rim, wrapping the
/// normalization shift into an equal-pair Status-OK lie in release
/// (and a debug panic). Expected = RN53(AGM(1, 3))·2^MIN, oracle
/// mantissa from mpmath 1.4.1 @500 bits via exact integer rounding.
#[test]
fn agm_bottom_rim_corner_is_exact_homogeneous() {
    let a = scaled(1, 53, i64::MIN);
    let b = scaled(3, 53, i64::MIN);
    let (r, st) = a.agm(&b, NE);
    assert_parts("agm bottom rim", &r, 17_188_830_925_994_225_664, i64::MIN);
    assert!(st.inexact(), "Status OK was the wrapped-shift defect");
}

/// ADR-0106 verifier refutation 2: rounding the normalized AGM at a
/// target coarser than the operands can carry past max(a, b); at the
/// top rim the carried binade is unrepresentable, the scale-back
/// saturates per the no-emax contract, and the flag must be
/// SURFACED (the first draft debug-asserted it away). Operands
/// (2 − 2^-69)·2^MAX and (2 − 2^-68)·2^MAX at p113: the true AGM
/// sits within 2^-53 of 2^(MAX+1).
#[test]
fn agm_top_rim_carry_surfaces_overflow() {
    let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
    let (a0, sa) = two.sub(&scaled(1, 113, -69), NE);
    assert!(sa.is_ok());
    let (b0, sb) = two.sub(&scaled(1, 113, -68), NE);
    assert!(sb.is_ok());
    let (a, s1) = a0.scale_by_pow2(i64::MAX);
    let (b, s2) = b0.scale_by_pow2(i64::MAX);
    assert!(s1.is_ok() && s2.is_ok());
    let (r_ne, st_ne) = a.agm_round(&b, 53, NE).unwrap();
    assert_parts("agm carry NE", &r_ne, 1 << 63, i64::MAX);
    assert!(
        st_ne.overflow() && st_ne.inexact(),
        "the scale-back saturation flag must be surfaced, got {st_ne:?}"
    );
    // The inward mode stays below the binade: exact scale-back, no
    // overflow; value = the all-ones 53-bit mantissa at MAX.
    let (r_tz, st_tz) = a.agm_round(&b, 53, RoundingMode::TowardZero).unwrap();
    assert_parts("agm carry TZ", &r_tz, u64::MAX << 11, i64::MAX);
    assert!(st_tz.inexact() && !st_tz.overflow());
}

/// ADR-0106 verifier refutation 3: no static normalization keeps a
/// spread ≳ 2^63 loop off the rim (the iterates converge to
/// exponent ≈ half-spread − log₂(half-spread), so near-convergence
/// products approach the full spread): the opposite-rim pair
/// certified a ~2^42-wrong value WITH STATUS OK. Spreads ≥ 2^33 now
/// take the asymptotic branch AGM(a, b) = π·a/(2·ln(4a/b)),
/// correctly roundable there (relative error O((b/a)²), far below
/// every expressible target). Oracle mantissas: π/((4s+4)·ln 2) at
/// exact big-integer s, mpmath 1.4.1 @500 bits, integer-rounded.
#[test]
fn agm_huge_spreads_take_the_asymptotic_branch() {
    // The verifier's measured-boundary case: s = 2^62 + 100.
    let s_exp = (1_i64 << 62) + 100;
    let (r, st) = scaled(1, 53, s_exp).agm(&scaled(1, 53, -s_exp), NE);
    assert_parts(
        "agm s=2^62+100",
        &r,
        10_450_910_948_271_020_032,
        4_611_686_018_427_387_942,
    );
    assert!(st.inexact());

    // The opposite-rim Status-OK lie: agm(2^(MAX−1), 2^(MIN+2)).
    let (r2, st2) = scaled(1, 53, i64::MAX - 1).agm(&scaled(1, 53, i64::MIN + 2), NE);
    assert_parts(
        "agm opposite rims",
        &r2,
        10_450_910_948_271_022_080,
        9_223_372_036_854_775_743,
    );
    assert!(st2.inexact(), "Status OK was the opposite-rim defect");
}

/// The revision verifier's amplification reproducer: SAME-SIGN rim
/// exponents at spread ~2^33. The first asymptotic draft computed
/// ln(big) and ln(small) whole; their near-cancellation amplified
/// the logs' absolute error ~2^31 past the charged half-width and
/// certified the wrong side of a constructed near-tie (1 ulp,
/// INEXACT). The exact exponent split (integer part of L exact,
/// mantissa logs O(1)) kills the amplification; the truth side is
/// pinned by mpmath @2600 bits (tie fraction 0.5000000000009…, so
/// NE rounds UP to …329).
#[test]
fn agm_same_sign_rim_spread_resolves_the_tie_side() {
    let a = scaled(1, 53, 1_i64 << 62);
    let (bm, sp) = BigFloat::parse_str("10384592579223522102099639059340099", 113, NE).unwrap();
    assert!(sp.is_ok(), "the 113-bit integer must parse exactly");
    let (b, sb) = bm.scale_by_pow2(4_611_686_009_837_453_199);
    assert!(sb.is_ok());
    let (r, st) = a.agm_round(&b, 53, NE).unwrap();
    assert_parts(
        "agm same-sign rim tie",
        &r,
        10_450_910_945_837_729_792,
        4_611_686_018_427_387_872,
    );
    assert!(st.inexact());
}

/// The loop/asymptotic seam: one pair just below the 2^33 spread
/// threshold (the loop) and one at it (the branch) — both match the
/// same oracle family (the asymptotic error 4^-s is astronomically
/// below the 53-bit ulp at s ≈ 2^32, so the two paths must agree
/// with the formula and, transitively, with each other).
#[test]
fn agm_spread_seam_is_coherent() {
    for (s_exp, mantissa, exponent) in [
        (
            1_i64 << 32,
            10_450_910_945_837_729_792_u64,
            4_294_967_264_i64,
        ),
        ((1_i64 << 32) - 2, 10_450_910_950_704_314_368, 4_294_967_262),
    ] {
        let (r, st) = scaled(1, 53, s_exp).agm(&scaled(1, 53, -s_exp), NE);
        assert_parts("agm seam", &r, mantissa, exponent);
        assert!(st.inexact());
    }
}

/// ADR-0095's extreme-spread case stays correct under the
/// normalization (the midpoint exponent is ~0 there, so the
/// normalized loop matches the old one; pinned against the
/// refinement-consistency oracle).
#[test]
fn agm_extreme_spread_control() {
    let a = scaled(1, 53, 1000);
    let b = scaled(1, 53, -1000);
    let (r53, st) = a.agm(&b, NE);
    assert!(st.inexact());
    let a_hi = scaled(1, 600, 1000);
    let b_hi = scaled(1, 600, -1000);
    let (r_hi, _) = a_hi.agm_round(&b_hi, 600, NE).unwrap();
    let (r_hi_53, _) = r_hi.round_to_precision(53, NE).unwrap();
    assert_eq!(
        r53.total_cmp(&r_hi_53),
        Ordering::Equal,
        "agm@53 disagrees with round(agm@600 -> 53)"
    );
}

// ---------------------------------------------------------------
// pf-a77o + pf-9wb2 (arc R2, slice R2.5, ADR-0107): the rim and
// ceiling hardening pair. (a) round_with_infinitesimal's residue
// placement saturated near i64::MIN, certifying Status-OK lies
// through every tiny-x dispatch (the ADR-0102 interim rim guards
// existed for exactly this; now removed — the computation lifts
// away from the rim and rounds in the lifted frame). (b) the
// Taylor loops' saturated r² terms parked a few bits below r
// instead of 2|e_r| below, corrupting sums near the rim by ~2^43
// ulps. (c) the driver's half-width exponent arithmetic could
// overflow (debug panic / release garbage). (d) parse_str at the
// documented u32::MAX precision ceiling wrapped its buffer width
// and panicked on valid input. Directions for the tiny-x rows are
// the ADR-0102/0104 series signs (|atan x| < |x|, sinh/atanh grow,
// truth strictly off-grid).
// ---------------------------------------------------------------

/// The rim guards are gone and the lifted infinitesimal rounds the
/// true side at any Normal exponent: atan at the bottom rim under
/// the inward mode takes pred, under the nearest mode the argument
/// — both INEXACT (Status OK / exact-zero results were the
/// verifier-recorded defects at e8b1284..fd2758c baselines).
#[test]
fn tiny_x_dispatches_are_sound_at_the_bottom_rim() {
    for k in [i64::MIN + 1, i64::MIN + 10, i64::MIN + 54] {
        let x = scaled(1, 53, k);
        let (r_tz, st_tz) = x.atan_round(53, RoundingMode::TowardZero).unwrap();
        assert_eq!(
            r_tz.total_cmp(&x.next_down().0),
            Ordering::Equal,
            "atan(2^(MIN+{})) TZ must be pred",
            k - i64::MIN
        );
        assert!(st_tz.inexact());
        let (r_ne, st_ne) = x.atan_round(53, NE).unwrap();
        assert_eq!(r_ne.total_cmp(&x), Ordering::Equal);
        assert!(st_ne.inexact(), "Status OK was the defect");
    }
    // The grow-direction family at the rim (sinh/atanh returned
    // (1 + 2^-9)·x with Status OK).
    let x = scaled(1, 53, i64::MIN + 10);
    for f in [BigFloat::sinh_round, BigFloat::atanh_round] {
        let (r, st) = f(&x, 53, RoundingMode::TowardZero).unwrap();
        assert_eq!(r.total_cmp(&x), Ordering::Equal, "grow-family TZ must be x");
        assert!(st.inexact());
        let (r_tp, _) = f(&x, 53, RoundingMode::TowardPositive).unwrap();
        assert_eq!(r_tp.total_cmp(&x.next_up().0), Ordering::Equal);
    }
    // hypot's dispatch at the rim (was (big, OK)).
    let big = scaled(1, 53, i64::MIN + 30);
    let small = scaled(1, 53, i64::MIN);
    let (rh, sth) = big.hypot_round(&small, 53, NE).unwrap();
    assert_eq!(rh.total_cmp(&big), Ordering::Equal);
    assert!(sth.inexact(), "falsely-exact was the defect");
}

/// The bottom-most borrow: atan(2^MIN)'s truth lies strictly inside
/// (0, `MinPos`), where the no-subnormal grid has nothing below. The
/// inward modes give +0 with UNDERFLOW|INEXACT; the nearest modes
/// give `MinPos` (the truth sits an infinitesimal below it, far above
/// the to-nearest midpoint) with INEXACT and no underflow
/// (after-rounding tininess, the exp-window convention).
#[test]
fn tiny_x_borrow_below_min_pos_is_mode_aware() {
    let x = scaled(1, 53, i64::MIN);
    let (r_tz, st_tz) = x.atan_round(53, RoundingMode::TowardZero).unwrap();
    assert!(r_tz.is_zero() && !r_tz.is_sign_negative());
    assert!(st_tz.underflow() && st_tz.inexact());
    let (r_ne, st_ne) = x.atan_round(53, NE).unwrap();
    assert_eq!(r_ne.total_cmp(&x), Ordering::Equal, "NE must be MinPos");
    assert!(st_ne.inexact() && !st_ne.underflow());
}

/// The saturated-Taylor-term corruption (pf-a77o site (b)):
/// sin(2^(MIN+10)) returned 0.9987·x (~2^43 ulps wrong) because the
/// r³ term's exponent clamped at `i64::MIN`, a few bits below r
/// instead of `2|e_r|`. The deep-tiny evaluations now return the
/// argument (the correct w-bit value); certification past the
/// driver's deep-rung ceiling falls back 1-ulp-honest (TZ pinned at
/// the fall-back side, the documented caveat at beyond-ceiling
/// exponent-encoded depth).
#[test]
fn sin_taylor_rim_term_no_longer_corrupts() {
    let x = scaled(1, 53, i64::MIN + 10);
    let (r_ne, st_ne) = x.sin_round(53, NE).unwrap();
    assert_eq!(
        r_ne.total_cmp(&x),
        Ordering::Equal,
        "sin NE must be the argument, got {r_ne}"
    );
    assert!(st_ne.inexact());
    // The rim-saturated atan2 quotient keeps y's sign now (was −0
    // for positive arguments through the same saturated arithmetic).
    let y = scaled(1, 53, i64::MIN + 10);
    let xq = scaled(1, 53, 1000);
    let (r2, st2) = y.atan2_round(&xq, 53, RoundingMode::TowardZero).unwrap();
    assert!(
        !r2.is_sign_negative() && !r2.is_zero(),
        "atan2 rim quotient must stay positive nonzero, got {r2}"
    );
    assert!(st2.inexact());
}

/// pf-9wb2, the named reproducer: parse_str("0.5", u32::MAX) panicked
/// at the "non-zero quotient" expect — the buffer-width add wrapped
/// at the documented precision ceiling. Release builds only: the
/// honest cost of a u32::MAX-precision request is a ~512 MB
/// computation (~0.6 s release, far slower in the debug matrix).
#[cfg(not(debug_assertions))]
#[test]
fn parse_str_at_the_precision_ceiling_is_exact() {
    let (v, st) = BigFloat::parse_str("0.5", u32::MAX, NE).unwrap();
    assert!(st.is_ok(), "0.5 is dyadic: the parse must be exact");
    let half = scaled(1, 53, -1);
    assert_eq!(
        v.partial_cmp(&half).0,
        Some(Ordering::Equal),
        "got a wrong value at the ceiling"
    );
}

/// ADR-0101 verifier round 2: exp2 at the exact integers past the
/// upstream dispatch's i64 magnitude cap. 2^(-2^63) = `MinPos` is
/// exactly representable: `Status::OK`, every mode. One deeper,
/// 2^(-2^63 - 1) sits EXACTLY at the `MinPos/2` tie: NE resolves to
/// +0 (the zero-significand-even convention), away/upward to
/// `MinPos`, with `UNDERFLOW|INEXACT`.
#[test]
fn exp2_exact_integers_past_the_i64_cap() {
    let x_min = scaled(1, 64, 63).negated();
    for mode in ALL_MODES {
        let (r, st) = x_min.exp2_round(53, mode).unwrap();
        assert_eq!(
            r.total_cmp(&min_pos(53)),
            Ordering::Equal,
            "2^(-2^63) must be exactly MinPos under {mode:?}, got {r}"
        );
        assert!(st.is_ok(), "exact power must be OK, got {st:?}");
    }
    let one = BigFloat::try_from_i64_exact(1, 66).unwrap();
    let (x_tie, sx) = x_min.sub(&one, NE);
    assert!(sx.is_ok());
    for mode in ALL_MODES {
        let (r, st) = x_tie.exp2_round(53, mode).unwrap();
        let expect_minpos = matches!(
            mode,
            RoundingMode::NearestAway | RoundingMode::TowardPositive
        );
        if expect_minpos {
            assert_eq!(r.total_cmp(&min_pos(53)), Ordering::Equal, "{mode:?}");
        } else {
            assert!(r.is_zero(), "{mode:?} must give +0, got {r}");
        }
        assert!(st.underflow() && st.inexact(), "{mode:?}: {st:?}");
    }
}

// ===============================================================
// Arc R3, slice R3.1 (pf-hkoj, ADR-0109): zeta_fe near the trivial
// zeros certified wrong VALUES. The trivial zeros zeta(-2n) = +0 are
// dispatched exactly, but the NEIGHBOURHOOD s = -2n - eps feeds the
// functional-equation branch, where sin(pi s/2) carries the whole
// result (zeta(-2n - eps) ~ -eps*zeta'(-2n), magnitude |eps|). The
// flat +96 working boost rounded s to the working width before sin
// was formed, collapsing the eps proximity; the half-width model was
// then violated and the first Ziv rung certified a value |e(eps)|
// orders too small. Fix: pre-boost by pole_proximity_depth(s) (the
// ADR-0098 lgamma/digamma reflection precedent), bounded by the
// input precision. Truths: mpmath 1.4.1, two-precision cross-checked,
// at the bit-identical exact input -2 - 2^-k.
// ===============================================================

/// `s = -2 - 2^-200` (exact at p264). zeta(-2-eps) is
/// +eps*zeta(3)/(4pi^2), positive; the broken kernel returned
/// 1.8949e-62 against the truth 1.89481e-62 (~5.5e-5 relative, ~39
/// wrong bits at p53), certified at the first rung. The true value's
/// p53 scaled mantissa has fraction above one half, so NE rounds up
/// and `TowardZero` truncates down: the
/// two modes differ by 1 ULP, and `assert_bit_exact` re-rounds the
/// reference per mode, making this a directed-rounding check on both.
/// Debug-cheap: the boost lifts the internal working precision to
/// ~target+96+200.
#[test]
fn zeta_fe_near_trivial_zero_shallow_certifies() {
    let neg2 = BigFloat::try_from_i64_exact(-2, 264).unwrap();
    let (eps, se) = BigFloat::try_from_i64_exact(1, 264)
        .unwrap()
        .scale_by_pow2(-200);
    assert!(se.is_ok());
    let (s, ss) = neg2.sub(&eps, NE);
    assert!(ss.is_ok(), "-2 - 2^-200 must be exact at p264");
    for mode in [NE, RoundingMode::TowardZero] {
        let (r, st) = s.zeta_round(53, mode).unwrap();
        assert!(
            !r.is_sign_negative() && !r.is_zero(),
            "zeta(-2-2^-200) must be a positive normal, got {r}"
        );
        assert_bit_exact(
            "zeta(-2-2^-200)",
            &r,
            "1.89481213461680241430670492056828379441453045590e-62",
            53,
            mode,
        );
        assert!(st.inexact());
    }
}

/// The deep case: `s = -2 - 2^-1200` (exact at p1264). The broken
/// kernel returned the same 9.944e-67 under both NE and `TowardZero`
/// (the collapse erased the mode distinction) where the truth is
/// 1.76836e-363 — ~296 orders of magnitude wrong, certified at the
/// first rung (the named pf-hkoj reproducer). Once fixed the two
/// modes differ by 1 ULP, each checked against its per-mode reference
/// rounding. The proximity
/// (1200 bits) exceeds every legacy Ziv working precision, so only
/// the pole_proximity_depth pre-boost recovers it. The result's
/// astronomical exponent is pinned structurally first, then the
/// mantissa bit-exactly.
///
/// Release builds only (the pf-jl35 zeta-pole precedent): the
/// ~1400-bit conditioning-boosted Borwein evaluations cost ~120 s in
/// the debug matrix but ~10 s in the MPFR full-union release job
/// (whose feature set covers this lane). The shallow row above
/// exercises the same boost mechanism debug-cheaply (~4 s).
#[cfg(not(debug_assertions))]
#[test]
fn zeta_fe_near_trivial_zero_deep_certifies() {
    let neg2 = BigFloat::try_from_i64_exact(-2, 1264).unwrap();
    let (eps, se) = BigFloat::try_from_i64_exact(1, 1264)
        .unwrap()
        .scale_by_pow2(-1200);
    assert!(se.is_ok());
    let (s, ss) = neg2.sub(&eps, NE);
    assert!(ss.is_ok(), "-2 - 2^-1200 must be exact at p1264");
    for mode in [NE, RoundingMode::TowardZero] {
        let (r, st) = s.zeta_round(53, mode).unwrap();
        assert!(
            !r.is_sign_negative() && !r.is_zero() && !r.is_infinite(),
            "zeta(-2-2^-1200) must be a tiny positive normal, got {r}"
        );
        assert_bit_exact(
            "zeta(-2-2^-1200)",
            &r,
            "1.76835922913628530304569635040298442138500151352e-363",
            53,
            mode,
        );
        assert!(st.inexact());
    }
}

/// Odd-integer control: at `s = -3 - 2^-200` the same proximity boost
/// fires (depth 200), but sin(pi s/2) PEAKS there rather than
/// vanishing (zeta(-3) = +1/120 is a nonzero rational), so there is
/// no cancellation and the value must stay the smooth zeta(-3)
/// rounding — confirming the boost only over-provisions harmlessly
/// near the odd integers. The eps perturbation sits ~140 bits below
/// the p53 ulp, so the result is RN53(1/120), INEXACT.
#[test]
fn zeta_fe_near_odd_integer_is_unperturbed() {
    let neg3 = BigFloat::try_from_i64_exact(-3, 264).unwrap();
    let (eps, se) = BigFloat::try_from_i64_exact(1, 264)
        .unwrap()
        .scale_by_pow2(-200);
    assert!(se.is_ok());
    let (s, ss) = neg3.sub(&eps, NE);
    assert!(ss.is_ok(), "-3 - 2^-200 must be exact at p264");
    let (r, st) = s.zeta_round(53, NE).unwrap();
    assert_bit_exact(
        "zeta(-3-2^-200)",
        &r,
        "0.00833333333333333333333333333333333333333333333333",
        53,
        NE,
    );
    assert!(st.inexact());
}

// ===============================================================
// Arc R3, slice R3.2 (pf-0r1l, ADR-0110): li / Ei / digamma at their
// deep INTERIOR zeros certified wrong values (Ei and digamma with the
// WRONG SIGN). The root cause is the input-encoded proximity to an
// irrational zero being rounded away before the kernel evaluates:
//   - Ei: ei_series rounds x to a flat working width and had no
//     near-zero cancellation handling, so it SATURATED (every deep
//     input returned the same ~-1.76e-74) and, when the input sat the
//     far side of the zero, returned the wrong sign.
//   - li / digamma: cancellation_boosted's 12-iteration LINEAR crawl
//     saturated at ~1492 bits, so inputs deeper than that never
//     reached the realised cancellation (geometric growth fixes it,
//     ADR-0110).
//   - digamma additionally hit the z_min = 2^28 shift bomb (~28 min)
//     once the Ziv working precision climbed past the Stirling table
//     reach; the Spouge-derivative dispatch removes it.
// Truths: mpmath 1.4.1, two-precision cross-checked, at the
// bit-identical p-bit roundings of each zero (round-trip decimals).
// Inputs are the closest p-bit dyadic to each zero (proximity ~2^-p),
// so parse_str(.., p, NE) recovers the identical dyadic mpmath used.
// ===============================================================

/// RN256 of Ei's zero x0 ~ 0.37250741... (83-digit round-trip).
const EI_ZERO_RN256: &str =
    "0.37250741078136663446199186658011913353568949777165405155565743524220012063620033168";
/// RN1500 of Ei's zero (457-digit round-trip; the wrong-sign case).
const EI_ZERO_RN1500: &str = "0.3725074107813666344619918665801191335356894977716540515556574352422001206362018543849260499515489423924647410089784888971884859964513190909730851441030323246757175996464553431492013427280264636400043516796895802963952541696002687969488523429259713357476646792904951276419139147433560923961567399937198579084410146164504254573975974766449866360285719219941913988999613476230720360214336149937872004295614699239838216968183789311782018826806090899072536810635";
/// RN256 of li's zero, the Ramanujan-Soldner constant mu ~ 1.4513692...
const LI_ZERO_RN256: &str =
    "1.4513692348833810502839684858920274494930322836480158630930045576624255957545243600";
/// RN2000 of li's zero (608-digit round-trip; the ~139-orders case).
/// Used only by the release-gated deep row below.
#[cfg(not(debug_assertions))]
const LI_ZERO_RN2000: &str = "1.4513692348833810502839684858920274494930322836480158630930045576624255957545178356595313577110868288470407515709706492030714335702042347848831900391108409842208865034838259375973333131396904110893875200271409693094298719306951483277641087216382996708422369126578452842707230772241061811131318318124512352191633404002577360159907445904017251368918339612155074307424793709893064395293020516241143642020809737507918258956251798522220055305868862328101012768516900105875608393150662138717113369733580730682011316746368408580505641371784257812013255078429486697408350823301767055762948648619812836567241287407428";
/// RN1500 of digamma's positive root r ~ 1.4616321... (457-digit).
/// Used only by the release-gated deep row below.
#[cfg(not(debug_assertions))]
const DG_ROOT_RN1500: &str = "1.461632144968362341262659542325721328468196204006446351295988408598786440353801810243074992733725592750556793365533053341617365778466985829177168381645024652542618792044384381978333559773961976074719431934937175414059451930109963724166527772172791673250880463960076932978144901475185803414306536810631010706016949785457933765577116113646852653864407737258989068226295819675052911994431197220725866405648207495227280806664927802646725469139476123636535751660";

/// Ei at the p256 rounding of its zero. The broken kernel SATURATED
/// (`ei_series` rounds x to ~target+128 and had no near-zero
/// cancellation boost), returning the same ~-1.76e-74 for every input
/// deeper than ~245 bits — here ~4 orders too large. Fixed by the
/// zero-window `cancellation_boosted` wrap. Debug-cheap (~1 ms).
#[test]
fn ei_at_zero_shallow_resolves_proximity() {
    let (x, _) = BigFloat::parse_str(EI_ZERO_RN256, 256, NE).unwrap();
    let (r, st) = x.ei_round(53, NE).unwrap();
    assert_bit_exact(
        "Ei(RN256 zero)",
        &r,
        "-5.93278017227464494659448410905778474751721681552e-78",
        53,
        NE,
    );
    assert!(st.inexact());
}

/// The Ei wrong-SIGN case: the p1500 rounding of the zero sits the
/// positive side, truth +1.74495e-452, but the saturated kernel
/// returned the negative ~-1.76e-74 (wrong sign AND 378 orders). The
/// fix preserves the input-encoded proximity, so the sign and
/// magnitude are both right. ~130 ms release / ~1.5 s debug.
#[test]
fn ei_at_zero_deep_keeps_the_sign() {
    let (x, _) = BigFloat::parse_str(EI_ZERO_RN1500, 1500, NE).unwrap();
    let (r, st) = x.ei_round(53, NE).unwrap();
    assert!(
        !r.is_sign_negative() && !r.is_zero(),
        "Ei(RN1500 zero) must be a positive normal (was wrong-sign), got {r}"
    );
    assert_bit_exact(
        "Ei(RN1500 zero)",
        &r,
        "1.74495224423398994478629394787165905009655137054e-452",
        53,
        NE,
    );
    assert!(st.inexact());
}

/// li shallow control: the p256 rounding of the Soldner constant was
/// already correct pre-fix (`cancellation_boosted` reached the ~258-bit
/// depth within its linear crawl). It must stay correct under the
/// geometric-growth change — a guard that the change preserves the
/// shallow path. Debug-cheap.
#[test]
fn li_at_soldner_shallow_control() {
    let (x, _) = BigFloat::parse_str(LI_ZERO_RN256, 256, NE).unwrap();
    let (r, st) = x.li_round(53, NE).unwrap();
    assert_bit_exact(
        "li(RN256 zero)",
        &r,
        "1.75145860178909040777584303367206605065164608427e-77",
        53,
        NE,
    );
    assert!(st.inexact());
}

/// li at the p2000 rounding of the Soldner constant: ~2000-bit
/// proximity, past the legacy linear crawl's ~1492 saturation, so the
/// broken kernel was ~139 orders wrong (NE and TZ identical — the tell
/// of a collapse). Geometric `cancellation_boosted` reaches the depth.
/// Release-gated (~2.8 s release / ~30 s+ debug).
#[cfg(not(debug_assertions))]
#[test]
fn li_at_soldner_deep_certifies() {
    let (x, _) = BigFloat::parse_str(LI_ZERO_RN2000, 2000, NE).unwrap();
    let (r, st) = x.li_round(53, NE).unwrap();
    assert!(
        !r.is_sign_negative() && !r.is_zero() && !r.is_infinite(),
        "li(RN2000 zero) must be a tiny positive normal, got {r}"
    );
    assert_bit_exact(
        "li(RN2000 zero)",
        &r,
        "3.01155405456412427246887850150149437719766046662e-603",
        53,
        NE,
    );
    assert!(st.inexact());
}

/// digamma at a high TARGET away from its root: digamma(2.5) at p1024
/// drives the Ziv working precision past the Stirling table reach, so
/// the broken kernel shifted up to `z_min` = 2^28 — a ~268M-term sum
/// costing ~28 minutes. The Spouge-derivative dispatch makes it ~0.3 s
/// (release). This guards the `spouge_digamma` path's correctness AND
/// its cost in the debug matrix (it is NOT release-gated precisely so
/// a spouge regression is caught fast); ~3-4 s debug.
#[test]
fn digamma_high_precision_uses_spouge() {
    let x = BigFloat::parse_str("2.5", 1024, NE).unwrap().0;
    let (r, st) = x.digamma_round(1024, NE).unwrap();
    // 330-digit reference (mpmath @1200 bits): assert_bit_exact at
    // p1024 needs ~310 digits to round-trip. This pins the value to
    // 1024 bits, catching the ~0.1·w S-sum cancellation loss the
    // internal absorption fixes (a too-short reference would mask it).
    assert_bit_exact(
        "digamma(2.5)@1024",
        &r,
        "0.703156640645243187225690333667911099473507062006232559619539412795011695949612564517992949382082542068032257375718212718335122180663862716377151165751591206551711191266986458548378646626288061477325050427202466808022754863088151740191007734417696879630107020674946735984012026572311494321432437157666747120129795446926129310661957",
        1024,
        NE,
    );
    assert!(st.inexact());
}

/// The Spouge sum cancellation GROWS with the argument (~0.1·w at
/// z≈2.5, ~0.4·w at z≈1e6, found by this slice's adversarial verifier),
/// so a fixed internal-precision margin under-covers large z and the
/// first draft of `spouge_digamma` certified `digamma(1000000)@1024`
/// ~190 bits short — a silent wrong value. The fix charges the measured
/// depth into the scale and re-drives through `cancellation_boosted` (the
/// whole Spouge regime, not just the root window). This large-z row
/// guards it; the small-z `digamma(2.5)@1024` row above did not.
/// Release-gated (~15 s release, dominated by the one-time
/// coefficient computation at the boosted precision).
#[cfg(not(debug_assertions))]
#[test]
fn digamma_large_argument_spouge_charges_sum_cancellation() {
    let x = BigFloat::parse_str("1000000.0", 1200, NE).unwrap().0;
    let (r, st) = x.digamma_round(1024, NE).unwrap();
    assert_bit_exact(
        "digamma(1000000)@1024",
        &r,
        "13.8155100579641907707746154031061852456026406778043880546126658109280907807403303326313511126554723511937285048862727304193036600881003872584088322029046881303856424221003292038345233650366929368941164152115961264901240562762800629222179573243109739999916937581752554743770949455676902220820754963461901359068818959949729948464555",
        1024,
        NE,
    );
    assert!(st.inexact());
}

/// The digamma wrong-SIGN + cost case: the p1500 rounding of the
/// positive root. The legacy linear crawl saturated below the ~1500
/// depth, returning the wrong sign, AND each deep evaluation hit the
/// `z_min` = 2^28 shift bomb (~28 min). Geometric `cancellation_boosted`
/// reaches the depth and the Spouge dispatch removes the bomb: truth
/// +7.84636e-453, ~2 s release. Release-gated (~25 s debug).
#[cfg(not(debug_assertions))]
#[test]
fn digamma_at_root_deep_keeps_the_sign() {
    let (x, _) = BigFloat::parse_str(DG_ROOT_RN1500, 1500, NE).unwrap();
    let (r, st) = x.digamma_round(53, NE).unwrap();
    assert!(
        !r.is_sign_negative() && !r.is_zero(),
        "digamma(RN1500 root) must be a positive normal (was wrong-sign), got {r}"
    );
    assert_bit_exact(
        "digamma(RN1500 root)",
        &r,
        "7.84635712238823226330795968835636586463977510063e-453",
        53,
        NE,
    );
    assert!(st.inexact());
}
