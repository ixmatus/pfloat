//! Fuzz target: forward + inverse hyperbolic dispatch.

#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;

use pfloat::{BigFloat, RoundingMode};

#[derive(Arbitrary, Debug)]
struct Input {
    op: u8,
    x: i64,
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

    match input.op % 6 {
        0 => {
            let _ = x.sinh(mode);
        }
        1 => {
            let _ = x.cosh(mode);
        }
        2 => {
            let _ = x.tanh(mode);
        }
        3 => {
            let _ = x.asinh(mode);
        }
        4 => {
            let _ = x.acosh(mode);
        }
        5 => {
            let _ = x.atanh(mode);
        }
        _ => unreachable!(),
    }
});
