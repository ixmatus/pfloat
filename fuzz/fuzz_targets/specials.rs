//! Fuzz target: tier-1 special functions.

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

    match input.op % 28 {
        0 => {
            let _ = x.erf(mode);
        }
        1 => {
            let _ = x.erfc(mode);
        }
        2 => {
            let _ = x.gamma(mode);
        }
        3 => {
            let _ = x.lgamma(mode);
        }
        4 => {
            let _ = x.digamma(mode);
        }
        5 => {
            let _ = x.beta(&y, mode);
        }
        6 => {
            let _ = x.agm(&y, mode);
        }
        7 => {
            let _ = x.ei(mode);
        }
        8 => {
            let _ = x.si(mode);
        }
        9 => {
            let _ = x.ci(mode);
        }
        10 => {
            let _ = x.li(mode);
        }
        11 => {
            let _ = x.ai(mode);
        }
        12 => {
            let _ = x.bi(mode);
        }
        13 => {
            let _ = x.ai_prime(mode);
        }
        14 => {
            let _ = x.bi_prime(mode);
        }
        15 => {
            let _ = x.j0(mode);
        }
        16 => {
            let _ = x.j1(mode);
        }
        17 => {
            // Fixed small order; `input.y` selects a bounded order
            // including the negative-order parity path.
            let _ = x.jn((input.y % 6) as i32, mode);
        }
        18 => {
            let _ = x.y0(mode);
        }
        19 => {
            let _ = x.y1(mode);
        }
        20 => {
            // Bounded order incl. the negative-order parity path;
            // `Y` is real only for `x > 0` (negative / zero `x`
            // exercise the INVALID / pole arms).
            let _ = x.yn((input.y % 6) as i32, mode);
        }
        21 => {
            let _ = x.i0(mode);
        }
        22 => {
            let _ = x.i1(mode);
        }
        23 => {
            // `I` is entire; bounded order incl. the negative-order
            // parity path (`I₋ₙ = Iₙ`, even, no sign) and the
            // negative-argument argument parity.
            let _ = x.in_((input.y % 6) as i32, mode);
        }
        24 => {
            let _ = x.k0(mode);
        }
        25 => {
            let _ = x.k1(mode);
        }
        26 => {
            // `K` is real only for `x > 0` (negative / zero / −0 `x`
            // exercise the INVALID / positive-pole arms); bounded
            // order incl. the negative-order parity path
            // (`K₋ₙ = Kₙ`, even, no sign).
            let _ = x.kn((input.y % 6) as i32, mode);
        }
        27 => {
            // ζ exercises every path: the pole `s = 1`
            // (DIV_BY_ZERO), `ζ(0) = −1/2` and the trivial zeros
            // `ζ(−2n) = 0` (special-cased), `s > 0` (Borwein), and
            // `s < 0` (the functional equation through Γ/sin/pow).
            let _ = x.zeta(mode);
        }
        _ => unreachable!(),
    }
});
