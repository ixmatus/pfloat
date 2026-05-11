//! Fuzz target: forward + inverse circular trig dispatch.

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

    match input.op % 7 {
        0 => {
            let _ = x.sin(mode);
        }
        1 => {
            let _ = x.cos(mode);
        }
        2 => {
            let _ = x.tan(mode);
        }
        3 => {
            let _ = x.asin(mode);
        }
        4 => {
            let _ = x.acos(mode);
        }
        5 => {
            let _ = x.atan(mode);
        }
        6 => {
            let _ = y.atan2(&x, mode);
        }
        _ => unreachable!(),
    }
});
