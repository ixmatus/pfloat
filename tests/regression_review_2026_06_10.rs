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
