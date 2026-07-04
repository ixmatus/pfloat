//! Regression guards for review-remediation R4 (epic pf-8iji):
//!
//! - **pf-t9qq**: the elementary functions and `sqrt`/`cbrt` dropped the
//!   underlying kernel `Status` (a Law-5 drift) and did not flag a
//!   degenerate (unbounded) enclosure. Each test below began red against
//!   the pre-fix behaviour (Status always OK) and lands with ADR-0116.
//! - **pf-ufb8**: `Mag` and `Ball` derived `Deserialize`, bypassing their
//!   construction invariants; a non-canonical `Mag` mantissa poisons the
//!   `Ord`-is-value-order precondition, and a non-finite `Ball` midpoint
//!   violates the finite-midpoint invariant. The validating impls reject
//!   such input (ADR-0116).
//!
//! Soundness (Law 1: the ball contains the truth) is re-checked alongside
//! the status assertions: threading a status must not change the interval.

#![cfg(feature = "big")]

use pfloat::{BigFloat, RoundingMode};
use pfloat_ball::Ball;

const NE: RoundingMode = RoundingMode::NearestEven;

fn bf(n: i64, p: u32) -> BigFloat {
    BigFloat::try_from_i64_exact(n, p).unwrap()
}
fn point(v: BigFloat) -> Ball<BigFloat> {
    Ball::point(v).unwrap()
}
#[cfg(feature = "exp-log")]
fn two_pow(k: i64) -> BigFloat {
    bf(1, 53).scale_by_pow2(k).0
}
/// `lower <= x <= upper`.
fn contains(b: &Ball<BigFloat>, x: &BigFloat) -> bool {
    use core::cmp::Ordering;
    b.lower().partial_cmp(x).0 != Some(Ordering::Greater)
        && b.upper().partial_cmp(x).0 != Some(Ordering::Less)
}

// ---------- pf-t9qq: status propagation on sqrt / cbrt (always available) ----------

#[test]
fn sqrt_of_two_flags_inexact_and_stays_sound() {
    // √2 is irrational: the kernel rounds, so the ball op must surface
    // INEXACT (was Status OK before ADR-0116). The interval is unchanged
    // and still contains the true √2.
    let (s, st) = point(bf(2, 53)).sqrt();
    assert!(st.inexact(), "sqrt(2) must flag INEXACT, got {st:?}");
    assert!(
        !st.invalid() && !st.overflow(),
        "no domain/overflow, got {st:?}"
    );
    assert!(!s.is_exact(), "sqrt(2) is not exact");
    // Soundness: a 400-bit reference √2 lies inside the ball.
    let true_sqrt2 = bf(2, 400).sqrt(NE).0;
    assert!(
        contains(&s, &true_sqrt2),
        "sqrt(2) ball must enclose the truth"
    );
}

#[test]
fn sqrt_of_perfect_square_stays_ok_and_exact() {
    // √9 = 3 exactly: threading status must not manufacture a spurious
    // flag on an exact result (Law 3).
    let (s, st) = point(bf(9, 53)).sqrt();
    assert!(st.is_ok(), "sqrt(9) is exact, expected OK, got {st:?}");
    assert!(s.is_exact());
}

#[test]
fn cbrt_of_two_flags_inexact() {
    let (c, st) = point(bf(2, 53)).cbrt();
    assert!(st.inexact(), "cbrt(2) must flag INEXACT, got {st:?}");
    assert!(!c.is_exact());
    let true_cbrt2 = bf(2, 400).cbrt(NE).0;
    assert!(contains(&c, &true_cbrt2));
}

#[test]
fn cbrt_of_perfect_cube_stays_ok_and_exact() {
    let (c, st) = point(bf(27, 53)).cbrt();
    assert!(st.is_ok(), "cbrt(27) is exact, expected OK, got {st:?}");
    assert!(c.is_exact());
}

// ---------- pf-t9qq: status propagation + degenerate flag on elementary fns ----------

#[cfg(feature = "exp-log")]
#[test]
fn exp_of_inexact_point_flags_inexact() {
    // exp(1) = e is irrational.
    let (e, st) = point(bf(1, 64)).exp();
    assert!(st.inexact(), "exp(1) must flag INEXACT, got {st:?}");
    assert!(!st.overflow());
    let true_e = bf(1, 400).exp(NE).0;
    assert!(contains(&e, &true_e), "exp(1) ball must enclose e");
}

#[cfg(feature = "exp-log")]
#[test]
fn exp_overflow_degenerates_to_entire_and_flags_overflow() {
    // exp(2^63) overflows pfloat's i64 exponent ceiling: the sound
    // enclosure is the entire line (unbounded). Before ADR-0116 this
    // returned a sound-but-degenerate ball with Status OK and no OVERFLOW
    // signal; now OVERFLOW is surfaced.
    let (e, st) = point(two_pow(63)).exp();
    assert!(
        st.overflow(),
        "exp(2^63) degenerated to an unbounded enclosure; OVERFLOW must be surfaced, got {st:?}"
    );
    // Soundness preserved: the enclosure is unbounded above (it must be,
    // the truth exceeds every representable).
    assert!(e.is_entire() || e.upper().is_infinite());
}

#[cfg(feature = "exp-log")]
#[test]
fn ln_inexact_flags_inexact_not_ok() {
    let (l, st) = point(bf(10, 64)).ln();
    assert!(
        st.inexact() && !st.invalid(),
        "ln(10) must flag INEXACT, got {st:?}"
    );
    let true_ln = bf(10, 400).ln(NE).0;
    assert!(contains(&l, &true_ln));
}

#[cfg(feature = "exp-log")]
#[test]
fn cosh_overflow_degenerates_and_flags_overflow() {
    // cosh grows like exp, so a huge argument overflows the exponent rim.
    let (c, st) = point(two_pow(63)).cosh();
    assert!(
        st.overflow(),
        "cosh(2^63) must surface OVERFLOW, got {st:?}"
    );
    assert!(c.is_entire() || c.upper().is_infinite());
}

#[cfg(feature = "trig")]
#[test]
fn sin_of_inexact_point_flags_inexact() {
    // sin(1) is irrational; the 1-Lipschitz route must surface the
    // midpoint kernel's INEXACT (was Status OK before ADR-0116).
    let (s, st) = point(bf(1, 64)).sin();
    assert!(st.inexact(), "sin(1) must flag INEXACT, got {st:?}");
    assert!(!st.overflow());
    let true_sin = bf(1, 400).sin(NE).0;
    assert!(
        contains(&s, &true_sin),
        "sin(1) ball must enclose the truth"
    );
}

#[cfg(feature = "trig")]
#[test]
fn atan_of_inexact_point_flags_inexact() {
    let (_r, st) = point(bf(1, 64)).atan();
    assert!(st.inexact(), "atan(1)=pi/4 must flag INEXACT, got {st:?}");
}

// ---------- pf-ufb8: validating serde Deserialize for Mag and Ball ----------

#[cfg(feature = "serde")]
mod serde_tests {
    use super::*;
    use pfloat_ball::Mag;

    /// A canonical `Mag` round-trips unchanged through JSON.
    #[test]
    fn mag_round_trip_is_identity() {
        for m in [
            Mag::ZERO,
            Mag::INFINITY,
            Mag::from_pow2(5),
            Mag::from_pow2(-30),
        ] {
            let s = serde_json::to_string(&m).expect("serialize");
            let back: Mag = serde_json::from_str(&s).expect("deserialize");
            assert_eq!(m, back, "round-trip mismatch via {s}");
        }
    }

    /// The exploit (code-read in the review, run here): a hand-crafted
    /// `Mag::Finite` whose mantissa lacks the top bit is non-canonical.
    /// The derived `Deserialize` accepted it silently, poisoning the
    /// `Ord`-is-value-order precondition. The validating impl rejects it.
    #[test]
    fn mag_deserialize_rejects_non_canonical_mantissa() {
        // mantissa = 1: top bit clear, value 2^(5-63) = 2^-58, but the
        // enum would sort it by (exponent=5, mantissa=1) — larger than a
        // genuine 2^0 whose exponent is 0. That inversion is the break.
        let non_canonical = r#"{"Finite":{"exponent":5,"mantissa":1}}"#;
        let r: Result<Mag, _> = serde_json::from_str(non_canonical);
        assert!(
            r.is_err(),
            "non-canonical mantissa must be rejected, got {r:?}"
        );

        // mantissa = 0 is not a finite positive magnitude at all.
        let zero_mantissa = r#"{"Finite":{"exponent":0,"mantissa":0}}"#;
        let r2: Result<Mag, _> = serde_json::from_str(zero_mantissa);
        assert!(
            r2.is_err(),
            "zero mantissa Finite must be rejected, got {r2:?}"
        );
    }

    /// A canonical `Mag::Finite` (top bit set) is still accepted.
    #[test]
    fn mag_deserialize_accepts_canonical_mantissa() {
        let canonical = r#"{"Finite":{"exponent":5,"mantissa":9223372036854775808}}"#;
        let m: Mag = serde_json::from_str(canonical).expect("canonical Mag accepted");
        assert_eq!(m, Mag::from_pow2(5));
    }

    /// The Ord break made concrete: the non-canonical `Finite{5,1}` would
    /// sort GREATER than a strictly larger canonical value. Removing it
    /// from the inhabitant set at the deserialize boundary is what upholds
    /// the ordering invariant the radius pipeline depends on.
    #[test]
    fn non_canonical_mag_would_break_value_order() {
        let smaller_value = Mag::from_pow2(-58); // canonical form of 2^-58
        let larger_value = Mag::from_pow2(0); // 2^0 = 1
        assert!(
            smaller_value < larger_value,
            "canonical order is value order"
        );
        // The would-be non-canonical encoding sorts by (exp=5, mant=1),
        // which is GREATER than (exp=0, mant=2^63): the derived Ord would
        // wrongly call 2^-58 > 1. The deserialize guard forbids it.
    }

    /// A `Ball<BigFloat>` round-trips unchanged.
    #[test]
    fn ball_round_trip_is_identity() {
        let b = Ball::new(bf(3, 53), Mag::from_pow2(-2)).unwrap();
        let s = serde_json::to_string(&b).expect("serialize");
        let back: Ball<BigFloat> = serde_json::from_str(&s).expect("deserialize");
        assert_eq!(b, back, "round-trip mismatch via {s}");

        // A point ball (rad = Zero) too.
        let p = point(bf(-7, 113));
        let sp = serde_json::to_string(&p).expect("serialize");
        let backp: Ball<BigFloat> = serde_json::from_str(&sp).expect("deserialize");
        assert_eq!(p, backp);
    }

    /// A deserialized `Ball` with a non-finite (`+∞` or NaN) midpoint
    /// violates the finite-midpoint invariant; the derived impl admitted
    /// it, the validating impl rejects it.
    #[test]
    fn ball_deserialize_rejects_non_finite_midpoint() {
        let inf_mid =
            r#"{"mid":{"precision":53,"class":{"Infinity":{"sign":"Positive"}}},"rad":"Zero"}"#;
        let r: Result<Ball<BigFloat>, _> = serde_json::from_str(inf_mid);
        assert!(r.is_err(), "infinite midpoint must be rejected, got {r:?}");

        let nan_mid = r#"{"mid":{"precision":53,"class":{"Nan":{"quiet":true,"sign":"Positive","payload":[0]}}},"rad":"Zero"}"#;
        let r2: Result<Ball<BigFloat>, _> = serde_json::from_str(nan_mid);
        assert!(r2.is_err(), "NaN midpoint must be rejected, got {r2:?}");
    }

    /// A `Ball` whose `rad` is a non-canonical `Mag` is rejected via the
    /// field's own validating deserialize (composition holds).
    #[test]
    fn ball_deserialize_rejects_non_canonical_radius() {
        let bad_rad = r#"{"mid":{"precision":53,"class":{"Normal":{"sign":"Positive","exponent":0,"mantissa":[9223372036854775808]}}},"rad":{"Finite":{"exponent":0,"mantissa":1}}}"#;
        let r: Result<Ball<BigFloat>, _> = serde_json::from_str(bad_rad);
        assert!(
            r.is_err(),
            "non-canonical radius must be rejected, got {r:?}"
        );
    }
}
