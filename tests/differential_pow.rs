//! MPFR differential: `BigFloat::pow` against `rug::Float::pow`,
//! **bit-exact** across all five IEEE rounding modes.
//!
//! Slice 7c routes `pow` through the Ziv interval-test driver
//! (ADR-0022): an exact integer exponent takes a square-and-multiply
//! fast path, every other case evaluates `exp(y · ln(x))`, and both
//! are correctly rounded to the target under the caller's mode.
//! `pow` is therefore the first transcendental off the
//! NearestEven-only differential tier — this lane compares for exact
//! equality (`assert_eq!`, the `differential_ei` idiom) under
//! [`BIT_EXACT_ROUNDING_MODES`].
//!
//! **Oracle note (NearestAway).** MPFR has no roundTiesToAway mode:
//! `MPFR_RNDA` (rug `Round::AwayZero`) is *directed* round-away (it
//! takes the farther neighbor of any inexact value), and `MPFR_RNDN`
//! is ties-to-*even*. An integer base raised to a small integer
//! power is frequently an exact tie at the target precision (e.g.
//! `99⁸ = 9227446944279201` is exactly between two p=53 values), so
//! neither `RNDA` nor `RNDN` is a valid roundTiesToAway oracle. This
//! lane therefore synthesizes the roundTiesToAway value from a
//! high-precision MPFR result and rounds it itself (nearest, ties to
//! the larger magnitude — all results here are positive). The shared
//! `differential::mpfr_round_of` still maps `NearestAway → AwayZero`;
//! that mapping is unused elsewhere (every other lane is
//! NearestEven-only) and is filed as separate cleanup, not widened
//! here.
//!
//! Restricted to positive bases and small finite exponents. The
//! full IEEE 754-2019 §9.2.1 table (zero base, infinity base,
//! negative base with integer / non-integer exponent) is covered
//! by the Kani harnesses in `src/verify/pow.rs`.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, rug_from_i64, sweep_size,
    BIT_EXACT_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};
use pfloat::RoundingMode;
use rug::float::Round;
use rug::ops::Pow;
use rug::Float;

/// Exact (for these inputs) high-precision `base^exp`, then the IEEE
/// roundTiesToAway value at precision `p`. `base ∈ [1,100]`,
/// `exp ∈ [-10,10]`: positive exponents give an exact integer (a
/// genuine tie source), negative give `1/integer` (non-terminating
/// binary, so never an exact `p`-bit tie). `p + 128` captures the
/// integer exactly and resolves every non-tie unambiguously, so the
/// tie test below is exact.
fn pow_ties_to_away(base: i64, exp: i64, p: u32) -> Float {
    let guard = p + 128;
    let b = Float::with_val(guard, base);
    let e = Float::with_val(guard, exp);
    let hp = Float::with_val(guard, Pow::pow(&b, &e));

    let (lo, _) = Float::with_val_round(p, &hp, Round::Zero);
    let (hi, _) = Float::with_val_round(p, &hp, Round::AwayZero);
    if lo == hi {
        return lo; // exactly representable at p
    }
    let d_lo = Float::with_val(guard, &hp - &lo).abs();
    let d_hi = Float::with_val(guard, &hi - &hp).abs();
    if d_hi < d_lo {
        hi
    } else if d_lo < d_hi {
        lo
    } else {
        hi // exact tie → away from zero; all results positive
    }
}

#[test]
fn pow_matches_mpfr_on_positive_base_small_exponent() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat7c");
    let cases = sweep_size().min(500);

    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            // Base in [1, 100], exponent in [-10, 10] keeps the
            // result inside i64 range without testing the
            // overflow / underflow paths.
            let base = next_i64_in(&mut state, 1, 100);
            let exp = next_i64_in(&mut state, -10, 10);
            for &mode in BIT_EXACT_ROUNDING_MODES {
                let bf_r = {
                    let b_bf = bigfloat_from_i64(base, p);
                    let e_bf = bigfloat_from_i64(exp, p);
                    let (r, _status) = b_bf.pow(&e_bf, mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = if matches!(mode, RoundingMode::NearestAway) {
                    pow_ties_to_away(base, exp, p)
                } else {
                    let b_rg = rug_from_i64(base, p);
                    let e_rg = rug_from_i64(exp, p);
                    let (r, _ord) =
                        Float::with_val_round(p, Pow::pow(&b_rg, &e_rg), mpfr_round_of(mode));
                    r
                };
                assert_eq!(bf_r, rug_r, "pow({base}, {exp}) at p={p}, mode={mode:?}");
            }
        }
    }
}
