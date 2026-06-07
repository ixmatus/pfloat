//! Tier-3 blocking lane: the Fundamental Theorem of Interval Arithmetic
//! (FTIA) as a pure-Rust self-consistency property, plus the ball edge
//! cases. Runs on every push under the `trig` feature.
//!
//! For a ball operation implementing real `f`, pick witnesses inside the
//! input ball(s), evaluate `f(witness)` with the trusted pfloat scalar
//! kernel at a much higher precision, and assert the witness result lies
//! in the output ball. The claim is deliberately bounded: because the
//! midpoint is produced by the same kernel that supplies the witness
//! oracle, this verifies that **the radius covers the kernel's own
//! residual**, not that the kernel is correct (that is Phase 1's job).
//! It is the blocking self-consistency lane, not the independent
//! soundness backstop (that is the per-release Arb containment lane,
//! ADR-0078).
//!
//! Required features: `trig` (the full elementary surface).

#![cfg(feature = "trig")]

use core::cmp::Ordering;
use pfloat::{BigFloat, RoundingMode, Status};
use pfloat_ball::{Ball, Mag};

mod common;
use common::{bf, contains, random_ball, witnesses, Rng};

const TN: RoundingMode = RoundingMode::TowardNegative;
const TP: RoundingMode = RoundingMode::TowardPositive;

/// Assert FTIA for a unary op: `f(w) ∈ result` for every witness `w`.
fn check_unary(
    a: &Ball<BigFloat>,
    result: &Ball<BigFloat>,
    f: impl Fn(&BigFloat, RoundingMode) -> (BigFloat, Status),
    label: &str,
) {
    if result.is_entire() {
        return; // entire encloses everything
    }
    for w in witnesses(a, 400) {
        // Bracket the true f(w): TN ≤ f(w) ≤ TP at 400 bits.
        let lo = f(&w, TN).0;
        let hi = f(&w, TP).0;
        if !lo.is_finite() || !hi.is_finite() {
            continue; // domain/overflow witness; the op handles it via flags
        }
        assert!(
            result.upper().partial_cmp(&lo).0 != Some(Ordering::Less),
            "{label}: result.upper < f(w) (lower bracket) — UNSOUND"
        );
        assert!(
            result.lower().partial_cmp(&hi).0 != Some(Ordering::Greater),
            "{label}: result.lower > f(w) (upper bracket) — UNSOUND"
        );
    }
}

#[test]
fn ftia_arithmetic_self_consistency() {
    let mut rng = Rng(0x1234_5678_9abc_def0);
    for _ in 0..1000 {
        let p = [24u32, 53, 113][(rng.next() % 3) as usize];
        let a = random_ball(&mut rng, p);
        let b = random_ball(&mut rng, p);
        // For each binary op, the true f(x,y) for the witness pair must be
        // enclosed. We fix y = b's witnesses and x = a's witnesses.
        let (sum, _) = a.add(&b);
        let (diff, _) = a.sub(&b);
        let (prod, _) = a.mul(&b);
        let (quot, _) = a.div(&b);
        for x in witnesses(&a, 400) {
            for y in witnesses(&b, 400) {
                if !sum.is_entire() {
                    let (lo, _) = x.add(&y, TN);
                    let (hi, _) = x.add(&y, TP);
                    assert!(contains_bracket(&sum, &lo, &hi), "add unsound");
                }
                if !diff.is_entire() {
                    let (lo, _) = x.sub(&y, TN);
                    let (hi, _) = x.sub(&y, TP);
                    assert!(contains_bracket(&diff, &lo, &hi), "sub unsound");
                }
                if !prod.is_entire() {
                    let (lo, _) = x.mul(&y, TN);
                    let (hi, _) = x.mul(&y, TP);
                    assert!(contains_bracket(&prod, &lo, &hi), "mul unsound");
                }
                if !quot.is_entire() && !y.is_zero() {
                    let (lo, _) = x.div(&y, TN);
                    let (hi, _) = x.div(&y, TP);
                    if lo.is_finite() && hi.is_finite() {
                        assert!(contains_bracket(&quot, &lo, &hi), "div unsound");
                    }
                }
            }
        }
    }
}

/// `result` contains the bracket `[lo, hi]` of a true value: `result`'s
/// own endpoints must straddle it.
fn contains_bracket(result: &Ball<BigFloat>, lo: &BigFloat, hi: &BigFloat) -> bool {
    result.lower().partial_cmp(hi).0 != Some(Ordering::Greater)
        && result.upper().partial_cmp(lo).0 != Some(Ordering::Less)
}

#[test]
fn ftia_unary_self_consistency() {
    let mut rng = Rng(0xfeed_face_cafe_babe);
    for _ in 0..700 {
        let p = [24u32, 53, 113][(rng.next() % 3) as usize];
        let a = random_ball(&mut rng, p);
        let (s, _) = a.sqrt();
        check_unary(&a, &s, BigFloat::sqrt, "sqrt");
        let (c, _) = a.cbrt();
        check_unary(&a, &c, BigFloat::cbrt, "cbrt");
        let (e, _) = a.exp();
        check_unary(&a, &e, BigFloat::exp, "exp");
        let (si, _) = a.sin();
        check_unary(&a, &si, BigFloat::sin, "sin");
        let (co, _) = a.cos();
        check_unary(&a, &co, BigFloat::cos, "cos");
        let (at, _) = a.atan();
        check_unary(&a, &at, BigFloat::atan, "atan");
        // ln only for positive balls.
        if a.lower().is_sign_positive() && !a.lower().is_zero() {
            let (l, _) = a.ln();
            check_unary(&a, &l, BigFloat::ln, "ln");
        }
    }
}

/// Dispatch the ball unary ops the seeded self-consistency lane covers.
fn ball_unary_ftia(a: &Ball<BigFloat>, fn_id: &str) -> Ball<BigFloat> {
    match fn_id {
        "exp" => a.exp().0,
        "sin" => a.sin().0,
        "cos" => a.cos().0,
        "atan" => a.atan().0,
        "ln" => a.ln().0,
        other => panic!("seeded ftia: unhandled {other}"),
    }
}

#[test]
fn ftia_unary_self_consistency_hard_to_round() {
    // Seed the self-consistency lane with the Lefevre-Muller / CORE-MATH
    // hard-to-round corpus: a midpoint where f(mid) sits near a rounding
    // boundary stresses the radius far more than a random integer midpoint.
    // Same bounded claim as ftia_unary_self_consistency (the radius covers the
    // kernel's residual), with the hardest available inputs. pf-vcqh.
    let mut rng = Rng(0x4c4d_5345_4544_0001);
    type ScalarKernel = fn(&BigFloat, RoundingMode) -> (BigFloat, Status);
    let unary: [(&str, ScalarKernel); 5] = [
        ("exp", BigFloat::exp),
        ("sin", BigFloat::sin),
        ("cos", BigFloat::cos),
        ("atan", BigFloat::atan),
        ("ln", BigFloat::ln),
    ];
    for (fn_id, kernel) in unary {
        let cases = common::lm_cases_for(fn_id).expect("seeded fn has a corpus");
        let mut used = 0u32;
        for &(xbits, _) in cases {
            if !common::is_finite_nonzero_f64(xbits) {
                continue;
            }
            let a = common::seeded_ball(&mut rng, xbits);
            let r = ball_unary_ftia(&a, fn_id);
            check_unary(&a, &r, kernel, fn_id);
            used += 1;
        }
        assert!(used >= 5, "{fn_id}: only {used} seeded hard-to-round balls");
    }
}

#[test]
fn seeded_midpoint_is_bit_exact() {
    // The bit-exact builder must reproduce the binary64 value exactly, so a
    // hard-to-round seed is not softened into a near-miss that would defeat
    // the point of seeding. Round-trip a spread of corpus inputs through
    // bf_of_f64_bits(bits, 53) and back to f64: the bits must match.
    for fn_id in ["exp", "sin", "ln", "tanh", "log2"] {
        let cases = common::lm_cases_for(fn_id).expect("seeded fn has a corpus");
        let mut checked = 0u32;
        for &(xbits, _) in cases.iter().take(8) {
            if !common::is_finite_nonzero_f64(xbits) {
                continue;
            }
            let mid = common::bf_of_f64_bits(xbits, 53);
            let (back, _) = mid.to_f64_round(RoundingMode::NearestEven);
            assert_eq!(
                back.to_bits(),
                xbits,
                "{fn_id}: seed 0x{xbits:016x} not reproduced bit-exact (got 0x{:016x})",
                back.to_bits()
            );
            checked += 1;
        }
        assert!(checked >= 5, "{fn_id}: only {checked} seeds round-tripped");
    }
}

// ---------- edge cases (slice 10 (d)) ----------

#[test]
fn degenerate_ball_reduces_to_scalar_kernel() {
    // An exact (radius-0) ball's op equals the scalar kernel's exact
    // result, with a zero radius when the kernel is exact.
    let a = Ball::point(bf(6, 53)).unwrap();
    let b = Ball::point(bf(7, 53)).unwrap();
    let (p, _) = a.mul(&b);
    assert!(p.is_exact());
    assert_eq!(
        p.midpoint().partial_cmp(&bf(42, 53)).0,
        Some(Ordering::Equal)
    );
    let (s, _) = Ball::point(bf(9, 53)).unwrap().sqrt();
    assert!(s.is_exact() && s.midpoint().partial_cmp(&bf(3, 53)).0 == Some(Ordering::Equal));
}

#[test]
fn zero_straddling_ball_mul_contains_zero() {
    // [0 ± 1] · [0 ± 1] must contain 0 and the extreme products ±1.
    let z = Ball::new(bf(0, 53), Mag::from_pow2(0)).unwrap();
    let (p, _) = z.mul(&z);
    assert!(contains(&p, &bf(0, 53)));
    assert!(contains(&p, &bf(1, 53)));
    assert!(contains(&p, &bf(-1, 53)));
}

#[test]
fn conversion_preserves_containment_both_directions() {
    // ball -> [lower, upper] -> from_interval must contain the original.
    let mut rng = Rng(0x0bad_c0de_0bad_c0de);
    for _ in 0..500 {
        let a = random_ball(&mut rng, 53);
        if a.is_entire() {
            continue;
        }
        let lo = a.lower();
        let hi = a.upper();
        let reboxed = Ball::from_interval(&lo, &hi).unwrap();
        // The re-boxed ball must contain the original's midpoint and both
        // of its endpoints (it re-encloses [lo, hi] ⊇ the original).
        assert!(contains(&reboxed, a.midpoint()));
        assert!(contains(&reboxed, &lo));
        assert!(contains(&reboxed, &hi));
    }
}

#[test]
fn entire_inputs_never_panic() {
    let entire = Ball::new(bf(0, 53), Mag::INFINITY).unwrap();
    let pt = Ball::point(bf(2, 53)).unwrap();
    // A sweep of ops with an entire operand must not panic.
    let _ = entire.add(&pt);
    let _ = entire.mul(&pt);
    let _ = pt.div(&entire);
    let _ = entire.sqrt();
    let _ = entire.exp();
    let _ = entire.ln();
    let _ = entire.sin();
    let _ = entire.atan();
    let _ = entire.cbrt();
}
