//! MPFR differential: `pfloat::constants::*` match the
//! correctly-rounded MPFR value of each constant bit-for-bit at every
//! tested precision and IEEE rounding mode.
//!
//! Each constant is driven through the Ziv interval test
//! (`crate::math::ziv`) inside pfloat, so the correctly-rounded claim
//! holds under all five modes. This lane is the rigorous cross-check
//! against an independent oracle (slice C3a, ADR-0067).
//!
//! `π`, `ln 2`, and `γ` use MPFR's own correctly-rounded constant
//! generators (`Constant::Pi`, `Constant::Log2`, `Constant::Euler`).
//! `π/2` and `π/4` are exact power-of-two scalings of `π`. The
//! composed constants (`2/π`, `ln 10`, `ln 2π`, `2/√π`) have no direct
//! MPFR generator, so the oracle evaluates them at a high working
//! precision and rounds down to the target, the standard
//! oracle-bracket approach; `mpfr_oracle_for_mode` then synthesizes
//! `NearestAway`, which MPFR lacks.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_to_rug, mpfr_oracle_for_mode, BIT_EXACT_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};
use pfloat::{constants, BigFloat, RoundingMode};
use rug::float::{Constant, Round};
use rug::Float;

/// Extra working bits the composed-constant oracle carries before
/// rounding to the target. 192 bits resolves every realistic target
/// unambiguously (the constants are irrational, so a 192-bit run at a
/// rounding boundary does not occur for the tested precisions).
const ORACLE_GUARD: u32 = 192;

/// Cross-check one pfloat constant against its MPFR oracle across the
/// transcendental precisions and all five IEEE rounding modes.
fn check<P, O>(name: &str, pf: P, op: O)
where
    P: Fn(u32, RoundingMode) -> BigFloat,
    O: Fn(u32, Round) -> Float,
{
    for &p in TRANSCENDENTAL_PRECISIONS {
        for &mode in BIT_EXACT_ROUNDING_MODES {
            let pf_rug = bigfloat_to_rug(&pf(p, mode));
            let oracle = mpfr_oracle_for_mode(&op, mode, p);
            assert_eq!(pf_rug, oracle, "{name} at p={p}, mode={mode:?}");
        }
    }
}

/// Round a high-precision MPFR value to `p` bits under `round`.
fn round_to(hi: &Float, p: u32, round: Round) -> Float {
    Float::with_val_round(p, hi, round).0
}

#[test]
fn pi_matches_mpfr() {
    check(
        "pi",
        |p, mode| constants::pi(p, mode).0,
        |p, round| Float::with_val_round(p, Constant::Pi, round).0,
    );
}

#[test]
fn pi_over_2_matches_mpfr() {
    check(
        "pi_over_2",
        |p, mode| constants::pi_over_2(p, mode).0,
        |p, round| {
            // π/2 shares π's mantissa; halving is an exact exponent
            // shift, so the directed rounding of π carries through.
            let (mut f, _) = Float::with_val_round(p, Constant::Pi, round);
            f >>= 1;
            f
        },
    );
}

#[test]
fn pi_over_4_matches_mpfr() {
    check(
        "pi_over_4",
        |p, mode| constants::pi_over_4(p, mode).0,
        |p, round| {
            let (mut f, _) = Float::with_val_round(p, Constant::Pi, round);
            f >>= 2;
            f
        },
    );
}

#[test]
fn two_over_pi_matches_mpfr() {
    check(
        "two_over_pi",
        |p, mode| constants::two_over_pi(p, mode).0,
        |p, round| {
            let w = p + ORACLE_GUARD;
            let pi_hi = Float::with_val(w, Constant::Pi);
            let two = Float::with_val(w, 2u32);
            let hi = Float::with_val(w, &two / &pi_hi);
            round_to(&hi, p, round)
        },
    );
}

#[test]
fn ln_2_matches_mpfr() {
    check(
        "ln_2",
        |p, mode| constants::ln_2(p, mode).0,
        |p, round| Float::with_val_round(p, Constant::Log2, round).0,
    );
}

#[test]
fn ln_10_matches_mpfr() {
    check(
        "ln_10",
        |p, mode| constants::ln_10(p, mode).0,
        |p, round| {
            let w = p + ORACLE_GUARD;
            let ten = Float::with_val(w, 10u32);
            let hi = Float::with_val(w, ten.ln_ref());
            round_to(&hi, p, round)
        },
    );
}

#[test]
fn euler_gamma_matches_mpfr() {
    check(
        "euler_gamma",
        |p, mode| constants::euler_gamma(p, mode).0,
        |p, round| Float::with_val_round(p, Constant::Euler, round).0,
    );
}

#[test]
fn two_over_sqrt_pi_matches_mpfr() {
    check(
        "two_over_sqrt_pi",
        |p, mode| constants::two_over_sqrt_pi(p, mode).0,
        |p, round| {
            let w = p + ORACLE_GUARD;
            let pi_hi = Float::with_val(w, Constant::Pi);
            let sqrt_pi = Float::with_val(w, pi_hi.sqrt_ref());
            let two = Float::with_val(w, 2u32);
            let hi = Float::with_val(w, &two / &sqrt_pi);
            round_to(&hi, p, round)
        },
    );
}

#[test]
fn ln_2pi_matches_mpfr() {
    check(
        "ln_2pi",
        |p, mode| constants::ln_2pi(p, mode).0,
        |p, round| {
            let w = p + ORACLE_GUARD;
            let pi_hi = Float::with_val(w, Constant::Pi);
            let two = Float::with_val(w, 2u32);
            let two_pi = Float::with_val(w, &two * &pi_hi);
            let hi = Float::with_val(w, two_pi.ln_ref());
            round_to(&hi, p, round)
        },
    );
}
