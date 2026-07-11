//! Algebraic-identity cross-checks for the complex elementary core (ADR-0092,
//! the C5 verification pass).
//!
//! These are CONSISTENCY checks, cheap and oracle-free: each holds to a few
//! ULPs because every operation rounds, so they catch gross errors (a wrong
//! branch, a wrong magnitude, a dropped term) without an external reference.
//! The bit-exact componentwise correct-rounding claim is pinned separately, by
//! the independent acb differential (`differential_acb.rs`); here the round
//! trips and the homomorphism law cross-tie `csqrt`, `cexp`, and `clog` to one
//! another so a defect in any one breaks an identity the others witness.
//!
//! The grid deliberately includes the cancellation regimes the kernels name as
//! their failure mode: `clog` near `|z| = 1`, `csqrt` near the negative-real
//! cut, and `cexp` near `y = kπ/2`, plus Lefevre-Muller hard-to-round seeds for
//! the scalar sub-kernels (`exp`, `ln`, `sin`, `cos`) the components compose
//! through.

#![cfg(feature = "trig")]

mod common;

use common::{bf, bf_of_f64_bits, close, is_finite_nonzero_f64, lm_inputs_for, NE};
use pfloat::BigFloat;
use pfloat_complex::Complex;

const P: u32 = 200;
/// Absolute tolerance `2^-(P-40) = 2^-160`. A round trip through two or three
/// roundings at `P = 200` bits accumulates error far below this, while a
/// wrong-branch or wrong-magnitude error is `O(1)`, so the band catches the
/// latter without false-failing on the former.
const SLACK: i64 = 40;

fn c(re: i64, im: i64) -> Complex<BigFloat> {
    Complex::new(bf(re, P), bf(im, P))
}

/// `re + im` from dyadic rationals `rn/rd + (in_/id) i` at precision `P`.
fn cq(rn: i64, rd: i64, in_: i64, id: i64) -> Complex<BigFloat> {
    let re = bf(rn, P).div(&bf(rd, P), NE).0;
    let im = bf(in_, P).div(&bf(id, P), NE).0;
    Complex::new(re, im)
}

/// Assert `got` is within tolerance of `want`, componentwise.
#[track_caller]
fn assert_close(got: &Complex<BigFloat>, want: &Complex<BigFloat>, who: &str) {
    assert!(
        close(&got.re, &want.re, P, SLACK),
        "{who}: re {} not within tol of {}",
        got.re,
        want.re
    );
    assert!(
        close(&got.im, &want.im, P, SLACK),
        "{who}: im {} not within tol of {}",
        got.im,
        want.im
    );
}

/// A deterministic grid of finite nonzero complex values, covering all four
/// quadrants, the axes, and a near-unit-circle / near-cut spread.
fn grid() -> Vec<Complex<BigFloat>> {
    let mut g = vec![
        c(3, 4),
        c(2, 3),
        c(-2, 3),
        c(-2, -3),
        c(2, -3),
        c(5, 0),
        c(0, 5),
        c(-5, 0),
        c(0, -5),
        c(1, 1),
        c(7, -2),
        cq(1, 2, 1, 3),
        cq(-1, 4, 5, 7),
        cq(11, 8, -3, 16),
    ];
    // Near the unit circle: |z| within 2^-50 of 1 (the clog cancellation regime).
    let tiny = bf(1, P).scale_by_pow2(-50).0;
    let one = bf(1, P);
    g.push(Complex::new(one.clone().add(&tiny, NE).0, tiny.clone()));
    g.push(Complex::new(one.clone(), tiny.clone()));
    g.push(Complex::new(one.sub(&tiny, NE).0, tiny.clone()));
    // Near the negative-real cut: -1 + tiny i (the csqrt cancellation regime).
    g.push(Complex::new(bf(-1, P), tiny.clone()));
    g.push(Complex::new(bf(-3, P), tiny));
    g
}

#[test]
fn csqrt_squared_round_trips() {
    // csqrt(z)^2 = z to rounding, every branch (including near the cut).
    for z in grid() {
        let w = z.sqrt(NE).0;
        let w2 = w.mul(&w, NE).0;
        assert_close(&w2, &z, "csqrt(z)^2");
    }
}

#[test]
fn cexp_clog_inverse() {
    // cexp(clog z) = z for z != 0 (clog = ln|z| + i*arg, cexp inverts it).
    for z in grid() {
        let lz = z.log(NE).0;
        let back = lz.exp(NE).0;
        assert_close(&back, &z, "cexp(clog z)");
    }
}

#[test]
fn clog_cexp_inverse_on_principal_strip() {
    // clog(cexp z) = z for z on the principal strip |Im z| < pi. The grid keeps
    // |im| <= 3 < pi, including im near pi/2 (the cexp cancellation regime).
    let half_pi = bf(1, P).atan2(&bf(0, P), NE).0; // atan2(1, 0) = pi/2
    let near_half_pi = half_pi.sub(&bf(1, P).scale_by_pow2(-40).0, NE).0;
    let strip = vec![
        c(0, 0).add(&c(2, 1), NE).0,
        c(1, -1),
        c(-2, 2),
        c(3, -3),
        Complex::new(bf(0, P), near_half_pi.clone()),
        Complex::new(bf(1, P), near_half_pi.negated()),
        cq(-3, 2, 5, 4),
    ];
    for z in strip {
        let ez = z.exp(NE).0;
        let back = ez.log(NE).0;
        assert_close(&back, &z, "clog(cexp z)");
    }
}

#[test]
fn cexp_additive_homomorphism() {
    // cexp(z + w) = cexp(z) * cexp(w).
    let g = grid();
    for (i, z) in g.iter().enumerate() {
        // Pair each z with a small w so z + w stays in a sane range.
        let w = &g[(i * 7 + 3) % g.len()];
        let zw = z.add(w, NE).0;
        let lhs = zw.exp(NE).0;
        let rhs = z.exp(NE).0.mul(&w.exp(NE).0, NE).0;
        assert_close(&lhs, &rhs, "cexp(z+w) vs cexp(z)cexp(w)");
    }
}

#[test]
fn modulus_of_cexp_is_exp_real_part() {
    // |cexp(x + iy)| = e^x exactly (to rounding), a check independent of the
    // round trips: it pins the real-part magnitude without inverting.
    for z in grid() {
        let w = z.exp(NE).0;
        let modulus = w.re.hypot(&w.im, NE).0;
        let ex = z.re.exp_round(P, NE).unwrap().0;
        assert!(
            close(&modulus, &ex, P, SLACK),
            "|cexp({}, {})| = {} not within tol of e^x = {}",
            z.re,
            z.im,
            modulus,
            ex
        );
    }
}

#[test]
fn lefevre_muller_seeded_round_trips() {
    // Seed the components with hard-to-round binary64 inputs for the scalar
    // sub-kernels (exp -> cexp's x, ln -> clog, sin/cos -> cexp's y), so each
    // composed component sits boundary-close. The identity still holds; a
    // wrong-branch or lost-bit defect at a hard case breaks it.
    let exp_in = lm_inputs_for("exp").expect("exp corpus");
    let sin_in = lm_inputs_for("sin").expect("sin corpus");
    // Bound the seeds so e^x does not overflow the round trip: take the small
    // ones (|x| <= 8) and pair with a hard sin-input imaginary part scaled into
    // the principal strip.
    let mut used = 0u32;
    for (k, &(xbits, _)) in exp_in.iter().enumerate() {
        if !is_finite_nonzero_f64(xbits) {
            continue;
        }
        let x = bf_of_f64_bits(xbits, P);
        // Skip large-magnitude exp inputs (the round trip would overflow).
        if matches!(
            x.abs().partial_cmp(&bf(8, P)).0,
            Some(core::cmp::Ordering::Greater)
        ) {
            continue;
        }
        let (ybits, _) = sin_in[k % sin_in.len()];
        if !is_finite_nonzero_f64(ybits) {
            continue;
        }
        // Scale the imaginary seed into (-pi/2, pi/2) by halving until small.
        let mut y = bf_of_f64_bits(ybits, P);
        while matches!(
            y.abs().partial_cmp(&bf(1, P)).0,
            Some(core::cmp::Ordering::Greater)
        ) {
            y = y.scale_by_pow2(-1).0;
        }
        let z = Complex::new(x, y);
        // cexp(clog z) = z (z != 0 by construction: x is finite nonzero).
        let back = z.log(NE).0.exp(NE).0;
        assert_close(&back, &z, "L-M seeded cexp(clog z)");
        used += 1;
    }
    assert!(used >= 10, "L-M seeding covered only {used} cases");
}
