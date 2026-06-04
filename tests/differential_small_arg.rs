//! MPFR differential: the small-argument fast-paths (ADR-0059) for
//! `atanh`, `asinh`, `sinh`, and `tanh` stay bit-exact correctly rounded
//! under every IEEE rounding mode across the tiny-x activation band and
//! its boundary.
//!
//! Each fast-path activates at `x.exponent <= -(target + 2)` and returns
//! `round_with_infinitesimal(x, …)` in place of the full Ziv
//! composition. The integer-input lanes (`differential_asinh` and
//! friends) sweep only `|x| >= 1`, so they never reach this band. This
//! lane sweeps two families of exact dyadic inputs against the MPFR
//! oracle, both signs, every mode:
//!
//! 1. `±2^exp` straddling each precision's activation edge — a few
//!    exponents above (still the Ziv path, a boundary-handoff check) and
//!    many at or below (the fast-path). Powers of two are the hardest
//!    case for `round_with_infinitesimal`: the residue add triggers a
//!    borrow renormalisation.
//! 2. `±M·2^exp` deep in the band for a spread of mantissa patterns `M`,
//!    so the full-mantissa (no-borrow) case is covered too.
//!
//! A mismatch on any input is the strict-revert trigger for that
//! function.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_pow2, bigfloat_to_rug, mpfr_oracle_for_mode, rug_pow2, BIT_EXACT_ROUNDING_MODES,
};
use pfloat::{BigFloat, RoundingMode};

/// Precisions whose `-(p + 2)` activation edges the sweep brackets.
const BAND_PRECISIONS: &[u32] = &[24, 53, 113, 200];

/// Offsets added to the activation edge `-(p + 2)`. Positive offsets sit
/// just above the edge (the Ziv path); zero and negative sit at or below
/// it (the fast-path), out to a deep-subnormal depth.
const EDGE_OFFSETS: &[i64] = &[3, 2, 1, 0, -1, -2, -5, -20, -80, -200];

/// Mantissa patterns for the deep-in-band sweep. All fit in 24 bits, so
/// `M·2^exp` is exact at every `BAND_PRECISIONS` entry: one bit, low
/// odds, alternating runs, and the full 24-bit width.
const MANTISSAS: &[u64] = &[1, 3, 5, 0x15, 0x7F, 0x5555, 0xAAAA, 0x80_0001, 0xFF_FFFF];

fn signed_bf(v: BigFloat, neg: bool) -> BigFloat {
    if neg {
        v.negated()
    } else {
        v
    }
}

fn signed_rug(v: rug::Float, neg: bool) -> rug::Float {
    if neg {
        -v
    } else {
        v
    }
}

/// `M·2^exp` at precision `p`, exact when `M` fits in `p` bits (the
/// multiply only shifts `M`'s exponent).
fn bf_scaled(m: u64, exp: i64, p: u32) -> BigFloat {
    let m_bf = BigFloat::try_from_i64_exact(m as i64, p).expect("precision >= 1");
    m_bf.mul(&bigfloat_pow2(exp, p), RoundingMode::NearestEven)
        .0
}

fn rug_scaled(m: u64, exp: i64, p: u32) -> rug::Float {
    rug_pow2(exp, p) * rug::Float::with_val(p, m)
}

/// Generate one band-sweep test per function: `$bf` is the `BigFloat`
/// method, `$rug` the matching `rug::Float` reference method.
macro_rules! band_lane {
    ($test:ident, $name:literal, $bf:ident, $rug:ident) => {
        #[test]
        fn $test() {
            // Compare kernel(x) against the MPFR oracle for every mode.
            // `mk_rug` rebuilds the exact dyadic input at the oracle's
            // working precision (which differs from `p` for NearestAway).
            let compare =
                |x: &BigFloat, mk_rug: &dyn Fn(u32) -> rug::Float, p: u32, label: &str| {
                    for &mode in BIT_EXACT_ROUNDING_MODES {
                        let (r, _status) = x.$bf(mode);
                        let bf_r = bigfloat_to_rug(&r);
                        let rug_r = mpfr_oracle_for_mode(
                            |prec, round| {
                                rug::Float::with_val_round(prec, mk_rug(prec).$rug(), round).0
                            },
                            mode,
                            p,
                        );
                        assert_eq!(bf_r, rug_r, "{} {label} mode={mode:?}", $name);
                    }
                };

            for &p in BAND_PRECISIONS {
                let edge = -(i64::from(p) + 2);

                // Family 1: ±2^exp straddling the activation edge.
                for &off in EDGE_OFFSETS {
                    let exp = edge + off;
                    for neg in [false, true] {
                        let x = signed_bf(bigfloat_pow2(exp, p), neg);
                        let s = if neg { "-" } else { "" };
                        let label = format!("{s}2^{exp}@p{p}");
                        compare(&x, &|prec| signed_rug(rug_pow2(exp, prec), neg), p, &label);
                    }
                }

                // Family 2: ±M·2^base deep in the band, full mantissas.
                // base = edge − 24 keeps exponent(M·2^base) ≤ edge for
                // every 24-bit M, so all inputs are in-band.
                let base = edge - 24;
                for &m_pat in MANTISSAS {
                    for neg in [false, true] {
                        let x = signed_bf(bf_scaled(m_pat, base, p), neg);
                        let s = if neg { "-" } else { "" };
                        let label = format!("{s}{m_pat:#x}*2^{base}@p{p}");
                        compare(
                            &x,
                            &|prec| signed_rug(rug_scaled(m_pat, base, prec), neg),
                            p,
                            &label,
                        );
                    }
                }
            }
        }
    };
}

band_lane!(atanh_band_matches_mpfr, "atanh", atanh, atanh_ref);
band_lane!(asinh_band_matches_mpfr, "asinh", asinh, asinh_ref);
band_lane!(sinh_band_matches_mpfr, "sinh", sinh, sinh_ref);
band_lane!(tanh_band_matches_mpfr, "tanh", tanh, tanh_ref);
