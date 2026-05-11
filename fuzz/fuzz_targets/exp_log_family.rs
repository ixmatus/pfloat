//! Fuzz target: exp/log family dispatch.
//!
//! Feeds an `(Op, x, y)` triple through each of `exp`, `expm1`,
//! `exp2`, `exp10`, `ln`, `log1p`, `log2`, `log10`, and `pow`.
//! Panic-freedom only. Accuracy is delegated to MPFR differential
//! per ADR-0013.

#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;

use pfloat::{BigFloat, RoundingMode};

#[derive(Arbitrary, Debug)]
struct Input {
    op: u8,
    x: i64,
    y: i64,
    mode: u8,
}

fn rounding_mode(byte: u8) -> RoundingMode {
    match byte % 5 {
        0 => RoundingMode::NearestEven,
        1 => RoundingMode::NearestAway,
        2 => RoundingMode::TowardZero,
        3 => RoundingMode::TowardPositive,
        4 => RoundingMode::TowardNegative,
        _ => unreachable!(),
    }
}

fuzz_target!(|input: Input| {
    let prec: u32 = 113;
    let mode = rounding_mode(input.mode);

    let Ok(x) = BigFloat::try_from_i64_exact(input.x, prec) else {
        return;
    };
    let Ok(y) = BigFloat::try_from_i64_exact(input.y, prec) else {
        return;
    };

    match input.op % 9 {
        0 => {
            let _ = x.exp(mode);
        }
        1 => {
            let _ = x.expm1(mode);
        }
        2 => {
            let _ = x.exp2(mode);
        }
        3 => {
            let _ = x.exp10(mode);
        }
        4 => {
            let _ = x.ln(mode);
        }
        5 => {
            let _ = x.log1p(mode);
        }
        6 => {
            let _ = x.log2(mode);
        }
        7 => {
            let _ = x.log10(mode);
        }
        8 => {
            let _ = x.pow(&y, mode);
        }
        _ => unreachable!(),
    }
});
