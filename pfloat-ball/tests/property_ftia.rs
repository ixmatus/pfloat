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

const NE: RoundingMode = RoundingMode::NearestEven;
const TN: RoundingMode = RoundingMode::TowardNegative;
const TP: RoundingMode = RoundingMode::TowardPositive;

/// Deterministic xorshift64 PRNG (no `rand` dependency; fixed seeds keep
/// the lane reproducible).
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    /// A signed integer in `[-range, range]`.
    fn int(&mut self, range: i64) -> i64 {
        (self.next() % (2 * range as u64 + 1)) as i64 - range
    }
}

fn bf(n: i64, p: u32) -> BigFloat {
    BigFloat::try_from_i64_exact(n, p).unwrap()
}

/// `lower <= x <= upper`.
fn contains(b: &Ball<BigFloat>, x: &BigFloat) -> bool {
    b.lower().partial_cmp(x).0 != Some(Ordering::Greater)
        && b.upper().partial_cmp(x).0 != Some(Ordering::Less)
}

/// A random ball at precision `p`: integer midpoint in `[-range, range]`
/// scaled by `2^scale`, radius `2^radexp` (or exact / entire).
fn random_ball(rng: &mut Rng, p: u32) -> Ball<BigFloat> {
    let m = rng.int(1 << 20);
    let scale = rng.int(40);
    let (mid, _) = bf(m, p).scale_by_pow2(scale);
    let rad = match rng.next() % 8 {
        0 => Mag::ZERO,
        _ => Mag::from_pow2(scale + rng.int(8) - 30),
    };
    Ball::new(mid, rad).unwrap()
}

/// Witnesses inside `[mid - rad, mid + rad]`, reconstructed EXACTLY (not
/// the outward-rounded `lower()`/`upper()`): `mid` and `mid ± rad·t` for
/// dyadic `t ∈ {0, ±1/2, ±1}` at high precision.
fn witnesses(b: &Ball<BigFloat>, work: u32) -> Vec<BigFloat> {
    let mid = b.midpoint().round_to_precision(work, NE).unwrap().0;
    let mut out = vec![mid.clone()];
    if let Mag::Finite { .. } = b.radius() {
        let rad = b
            .radius()
            .to_bigfloat()
            .round_to_precision(work, NE)
            .unwrap()
            .0;
        for &(num, den_pow) in &[(1i64, 0u32), (1, 1)] {
            let (scaled, _) = rad.scale_by_pow2(-(den_pow as i64));
            let (scaled, _) = scaled.mul(&bf(num, work), NE);
            out.push(mid.add(&scaled, NE).0);
            out.push(mid.sub(&scaled, NE).0);
        }
    }
    out
}

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
