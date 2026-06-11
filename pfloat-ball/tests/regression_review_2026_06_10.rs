//! Regression guards for the 2026-06-10 workspace deep review,
//! pfloat-ball slice (epic pf-8iji: pf-m37w bracket-saturation
//! unsoundness, pf-1bqy parse_decimal collapse). Law 1 (the ball
//! contains the truth) is the property under test; each test began
//! red and lands with its fix (ADR-0099).

#![cfg(feature = "big")]

use pfloat::BigFloat;
use pfloat_ball::Ball;

fn pow2_point(k: i64) -> Ball<BigFloat> {
    let (m, s) = BigFloat::try_from_i64_exact(1, 53)
        .unwrap()
        .scale_by_pow2(k);
    assert!(s.is_ok());
    Ball::point(m).unwrap()
}

/// pf-m37w: squaring a point ball at 2^k with 2k past i64::MAX
/// saturates all three directed products to the same clamped finite
/// value; the bracket spread vanishes and the "exact" zero-radius
/// ball EXCLUDES the truth (~2^(2k) > MaxFinite). A sound enclosure
/// must be unbounded above. NOTE: the review's reproducer used
/// k = 4611686018427387900, whose square (exponent 2k = i64::MAX−7)
/// is REPRESENTABLE and exactly computed — another misadjudicated
/// reproducer; k = 2^62 genuinely crosses the rim.
#[test]
fn ball_mul_overflow_saturation_widens() {
    let b = pow2_point(4_611_686_018_427_387_904);
    let (r, st) = b.mul(&b);
    assert!(
        !r.is_exact(),
        "saturated product claimed exactness: mid={} rad={:?}",
        r.midpoint(),
        r.radius()
    );
    assert!(
        r.is_entire() || r.upper().is_infinite(),
        "truth exceeds every representable; the upper end must be unbounded (got upper={})",
        r.upper()
    );
    assert!(st.overflow(), "OVERFLOW must be surfaced, got {st:?}");
}

/// The underflow mirror: 2^k with 2k below i64::MIN clamps the
/// products UP to the exponent floor, overstating the magnitude;
/// the ball must stretch down past the truth (which lies between 0
/// and the clamped value).
#[test]
fn ball_mul_underflow_saturation_widens() {
    let b = pow2_point(-4_611_686_018_427_387_905);
    let (r, st) = b.mul(&b);
    assert!(!r.is_exact(), "saturated product claimed exactness");
    // The truth is a positive value below the clamped midpoint;
    // soundness needs lower() <= truth. With no representable
    // positive below MinPos, reaching zero (or below) suffices.
    let lower = r.lower();
    assert!(
        lower.is_zero() || lower.is_sign_negative(),
        "lower end must reach zero, got {lower}"
    );
    assert!(st.underflow(), "UNDERFLOW must be surfaced, got {st:?}");
}

/// Division through the same rim: huge / tiny overflows the result
/// exponent.
#[test]
fn ball_div_overflow_saturation_widens() {
    let num = pow2_point(4_611_686_018_427_387_904);
    let den = pow2_point(-4_611_686_018_427_387_905);
    let (r, st) = num.div(&den);
    assert!(!r.is_exact(), "saturated quotient claimed exactness");
    assert!(
        r.is_entire() || r.upper().is_infinite(),
        "upper end must be unbounded, got {}",
        r.upper()
    );
    assert!(st.overflow(), "OVERFLOW must be surfaced, got {st:?}");
}

/// pf-1bqy: pfloat's parse saturates past-budget tiny magnitudes to
/// +0 (mode-blind: pf-mw6u, separate arc) but flags UNDERFLOW; the
/// directed-pair interval then collapsed to the exact [0, 0] ball,
/// excluding the positive truth. The ball must keep 0 ≤ truth ≤
/// upper with a positive upper end at least the value itself
/// (witness: 2^-3321932 ≤ 0.1e-1000000).
#[test]
fn ball_parse_decimal_past_cap_tiny_is_contained() {
    let b = Ball::<BigFloat>::parse_decimal("0.1e-1000000", 53).unwrap();
    assert!(!b.is_exact(), "collapsed [0,0] claimed exactness");
    let (witness, s) = BigFloat::try_from_i64_exact(1, 53)
        .unwrap()
        .scale_by_pow2(-3_321_932);
    assert!(s.is_ok());
    let upper = b.upper();
    assert!(
        matches!(
            upper.partial_cmp(&witness).0,
            Some(core::cmp::Ordering::Greater | core::cmp::Ordering::Equal)
        ),
        "upper end {upper} excludes the truth (>= {witness})"
    );
    let lower = b.lower();
    assert!(
        lower.is_zero() || lower.is_sign_negative(),
        "lower end must not exceed the tiny truth, got {lower}"
    );

    // The negative literal mirrors.
    let bn = Ball::<BigFloat>::parse_decimal("-0.1e-1000000", 53).unwrap();
    assert!(!bn.is_exact());
    let lower_n = bn.lower();
    assert!(
        matches!(
            lower_n.partial_cmp(&witness.negated()).0,
            Some(core::cmp::Ordering::Less | core::cmp::Ordering::Equal)
        ),
        "negative lower end {lower_n} excludes the truth"
    );

    // Control: an in-budget parse stays tight and inexact-only.
    let c = Ball::<BigFloat>::parse_decimal("0.5", 53).unwrap();
    assert!(c.is_exact(), "0.5 is dyadic; the interval is a point");
}

/// Control (also the review's BL3 input, which never crossed the
/// rim): a representable power-of-two square stays an exact ball.
#[test]
fn ball_mul_review_input_is_genuinely_exact() {
    let b = pow2_point(4_611_686_018_427_387_900);
    let (r, st) = b.mul(&b);
    assert!(r.is_exact(), "2^(2k) at 2k = i64::MAX-7 is exact");
    assert!(st.is_ok(), "no saturation, got {st:?}");
}

/// Control: ordinary products keep their tight radius and flags.
#[test]
fn ball_mul_generic_unchanged() {
    let a = Ball::point(BigFloat::try_from_i64_exact(10, 53).unwrap()).unwrap();
    let (p, st) = a.mul(&a);
    assert!(p.is_exact(), "10*10 = 100 exactly");
    assert!(!st.overflow() && !st.underflow());
}

/// The adversarial verifier's refutation of this slice's first
/// draft: the RADIUS pipeline had the same discarded-status holes
/// one level down. Here the three directed midpoint products are
/// exact (no widening trigger), but the propagated-radius term
/// |a|·rb (true exponent i64::MAX + 10) overflow-clamps DOWN with
/// its OVERFLOW discarded, and the under-sized radius excluded a
/// member of the product set with Status OK.
#[test]
fn ball_mul_radius_term_overflow_widens() {
    use pfloat_ball::Mag;
    let a = pow2_point(i64::MAX - 10);
    let b = Ball::new(
        pfloat::BigFloat::try_from_i64_exact(1, 53).unwrap(),
        Mag::from_pow2(20),
    )
    .unwrap();
    let (r, _) = a.mul(&b);
    // The product set contains 2^(MAX-10)·(1 + 2^20) ≈ 2^(MAX+10),
    // beyond every representable: the enclosure must be unbounded.
    assert!(
        r.is_entire() || r.upper().is_infinite(),
        "radius-term overflow must widen to an unbounded enclosure"
    );
}

/// The div mirror: the denominator LOWER bound blo·|b| (true
/// exponent below i64::MIN) underflow-clamps UP, under-sizing prop;
/// the representable truth 1/blo = 2^(−k0+10) escaped the ball with
/// Status OK.
#[test]
fn ball_div_denominator_underflow_widens() {
    use pfloat_ball::Mag;
    let k0: i64 = -4_611_686_018_427_387_900;
    let (mid, s) = pfloat::BigFloat::try_from_i64_exact(1, 53)
        .unwrap()
        .scale_by_pow2(k0);
    assert!(s.is_ok());
    // radius just under the midpoint: blo = 2^(k0-10) stays positive.
    let (rad_bf, s2) = pfloat::BigFloat::try_from_i64_exact(1023, 53)
        .unwrap()
        .scale_by_pow2(k0 - 10);
    assert!(s2.is_ok());
    let den = Ball::new(mid, Mag::from_bigfloat_ceil(&rad_bf)).unwrap();
    let num = Ball::point(pfloat::BigFloat::try_from_i64_exact(1, 53).unwrap()).unwrap();
    let (r, _) = num.div(&den);
    // 1/blo = 2^(-(k0-10)) is representable and in the quotient set.
    let (truth_member, s3) = pfloat::BigFloat::try_from_i64_exact(1, 53)
        .unwrap()
        .scale_by_pow2(-(k0 - 10));
    assert!(s3.is_ok());
    let upper = r.upper();
    assert!(
        upper.is_infinite()
            || matches!(
                upper.partial_cmp(&truth_member).0,
                Some(core::cmp::Ordering::Greater | core::cmp::Ordering::Equal)
            ),
        "quotient-set member {truth_member} escaped the ball (upper {upper})"
    );
}
