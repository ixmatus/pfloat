//! Cross-check `BigFloat::to_f32_round` / `to_f64_round` (the pure-Rust
//! `no_std` conversion in `src/convert.rs`) against the rug/MPFR
//! reference over a broad sweep of off-grid values and all rounding
//! modes.
//!
//! The unit tests in `src/convert.rs` cover the round-trip identity and
//! hand-derived corner cases. This lane is the independent oracle: it
//! builds values that do not sit on the target grid (so a real rounding
//! decision is made) and asserts pfloat's pure-Rust result is bit-exact
//! with MPFR's for the same value under every mode MPFR provides a
//! primitive for. The `f32` lane reuses the harness's
//! `certified_round_bf_to_f32` (which itself routes through MPFR and
//! synthesizes `NearestAway`), so it covers all five modes; the `f64`
//! lane checks the four MPFR-native modes directly, with `NearestAway`
//! covered by the width-generic `f32` lane.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod oracle;

use oracle::convert::{bigfloat_to_rug, certified_round_bf_to_f32};
use pfloat::{BigFloat, RoundingMode};
use rug::float::Round;

const MODES: [RoundingMode; 5] = [
    RoundingMode::NearestEven,
    RoundingMode::NearestAway,
    RoundingMode::TowardZero,
    RoundingMode::TowardPositive,
    RoundingMode::TowardNegative,
];

const ITERS: usize = 200_000;

/// Deterministic 64-bit LCG (Knuth MMIX constants). Pure-function
/// generator so the sweep is reproducible without a `rand` dependency.
fn next(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state
}

/// A high-precision off-grid value: the 128-bit quotient of two random
/// finite f64s, which lands off both the f32 and the f64 grid.
fn ratio(a_bits: u64, b_bits: u64) -> Option<BigFloat> {
    let a = f64::from_bits(a_bits);
    let b = f64::from_bits(b_bits);
    if !a.is_finite() || !b.is_finite() || b == 0.0 {
        return None;
    }
    let a128 = BigFloat::from_f64(a)
        .round_to_precision(128, RoundingMode::NearestEven)
        .unwrap()
        .0;
    let b128 = BigFloat::from_f64(b)
        .round_to_precision(128, RoundingMode::NearestEven)
        .unwrap()
        .0;
    Some(a128.div(&b128, RoundingMode::NearestEven).0)
}

#[test]
fn to_f32_round_matches_mpfr() {
    let mut state = 0x1234_5678_9abc_def0u64;
    for _ in 0..ITERS {
        // Two value families: p=53 from a random f64 (off the f32 grid),
        // and a p=128 ratio (off every grid).
        let direct = {
            let x = f64::from_bits(next(&mut state));
            if x.is_finite() {
                Some(BigFloat::from_f64(x))
            } else {
                None
            }
        };
        let ratio = ratio(next(&mut state), next(&mut state));

        for bf in [direct, ratio].into_iter().flatten() {
            for mode in MODES {
                let mine = bf.to_f32_round(mode).0.to_bits();
                let reference =
                    certified_round_bf_to_f32(&bf, mode).expect("finite value rounds to Some");
                assert_eq!(
                    mine, reference,
                    "to_f32_round {mode:?}: pfloat 0x{mine:08x} != mpfr 0x{reference:08x}"
                );
            }
        }
    }
}

#[test]
fn to_f64_round_matches_mpfr() {
    let mut state = 0x0fed_cba9_8765_4321u64;
    let native = [
        (RoundingMode::NearestEven, Round::Nearest),
        (RoundingMode::TowardZero, Round::Zero),
        (RoundingMode::TowardPositive, Round::Up),
        (RoundingMode::TowardNegative, Round::Down),
    ];
    for _ in 0..ITERS {
        let Some(bf) = ratio(next(&mut state), next(&mut state)) else {
            continue;
        };
        let rug_val = bigfloat_to_rug(&bf);
        for (mode, round) in native {
            let mine = bf.to_f64_round(mode).0.to_bits();
            let reference = rug_val.to_f64_round(round).to_bits();
            assert_eq!(
                mine, reference,
                "to_f64_round {mode:?}: pfloat 0x{mine:016x} != mpfr 0x{reference:016x}"
            );
        }
    }
}
