//! Differential lane for binary64.
//!
//! The 2^64 input space cannot be enumerated, so f64 rests on a
//! structured random sample plus an edge battery, each input driven
//! through the shell and certified against the MPFR oracle under all
//! five rounding modes via `verify_input`. Domain errors and poles are
//! not filtered out: the oracle and the flag gate handle them (a
//! negative `ln` input certifies NaN + `INVALID`), so feeding any
//! finite value is sound. The two binary functions are exercised at
//! both widths. Mirrors pfloat's `differential_*.rs` lanes, with the
//! shell's hardware output as the comparison subject.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod harness;

use harness::{
    next_f64_banded, next_u64, sweep_size, verify_input, Hw, LibmArg, LibmFnId, StatusGate,
    Verdict, ALL_MODES,
};

const GATE: StatusGate = StatusGate::ValueAndDomainHard;

/// A verdict the differential lane accepts: a clean match, or the
/// measure-zero hard-to-round straddle (recorded by the exhaustive
/// sweep, not a failure here).
fn accept<H: Hw>(f: LibmFnId, input: H::Bits, arg: LibmArg, fails: &mut Vec<String>, label: &str) {
    for &mode in ALL_MODES {
        let v = verify_input::<H>(f, input, arg, mode, GATE);
        if !matches!(v, Verdict::Ok | Verdict::OracleInconclusive { .. }) {
            fails.push(format!("{label} {mode:?}: {v:?}"));
            if fails.len() > 50 {
                return;
            }
        }
    }
}

/// Domain-relevant f64 input for `f`. Bounded-domain functions
/// concentrate the sample in their interesting range (with a little
/// spill past the boundary to exercise `INVALID`); the rest span a wide
/// exponent band.
fn gen(f: LibmFnId, st: &mut u64) -> f64 {
    let unit = |st: &mut u64| (next_u64(st) >> 11) as f64 / (1u64 << 53) as f64; // [0, 1)
    let sign = |st: &mut u64| if next_u64(st) & 1 == 0 { 1.0 } else { -1.0 };
    match f {
        LibmFnId::Asin | LibmFnId::Acos => sign(st) * unit(st) * 1.1,
        LibmFnId::Atanh => sign(st) * unit(st) * 1.05,
        LibmFnId::Acosh => 1.0 + next_f64_banded(st, -40, 40).abs(),
        _ => next_f64_banded(st, -60, 60),
    }
}

#[test]
fn edge_battery_all_functions() {
    let edges: [f64; 13] = [
        0.0,
        -0.0,
        f64::from_bits(0x0000_0000_0000_0001), // smallest subnormal
        f64::from_bits(0x000F_FFFF_FFFF_FFFF), // largest subnormal
        f64::from_bits(0x0010_0000_0000_0000), // smallest normal
        1.0,
        -1.0,
        2.0,
        0.5,
        f64::MAX,
        -f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ];
    let mut fails = Vec::new();
    for &f in LibmFnId::UNARY {
        for &x in &edges {
            accept::<f64>(
                f,
                x.to_bits(),
                LibmArg::None,
                &mut fails,
                &format!("{} edge {x}", f.name()),
            );
        }
        // NaN propagation.
        accept::<f64>(
            f,
            f64::NAN.to_bits(),
            LibmArg::None,
            &mut fails,
            &format!("{} NaN", f.name()),
        );
    }
    assert!(fails.is_empty(), "edge battery failures: {fails:#?}");
}

#[test]
fn random_sample_unary() {
    let n = sweep_size().min(500);
    let mut fails = Vec::new();
    for &f in LibmFnId::UNARY {
        let mut st = 0x6c69_626d_5f64_6966u64 ^ (f.name().len() as u64).wrapping_mul(0x9E37);
        for _ in 0..n {
            let x = gen(f, &mut st);
            accept::<f64>(
                f,
                x.to_bits(),
                LibmArg::None,
                &mut fails,
                &format!("{} rand {x:e}", f.name()),
            );
            if fails.len() > 50 {
                break;
            }
        }
    }
    assert!(fails.is_empty(), "random unary failures: {fails:#?}");
}

#[test]
fn hypot_both_widths() {
    let n = sweep_size().min(400);
    let mut fails = Vec::new();
    let mut st = 0x6879_706f_745f_6c6du64;
    for _ in 0..n {
        let x = next_f64_banded(&mut st, -60, 60);
        // Every fourth pair uses near-equal operands to stress the
        // sum-of-squares with no cancellation; otherwise an independent y.
        let y = if next_u64(&mut st) & 3 == 0 {
            f64::from_bits(x.to_bits().wrapping_add(1))
        } else {
            next_f64_banded(&mut st, -60, 60)
        };
        accept::<f64>(
            LibmFnId::Hypot,
            x.to_bits(),
            LibmArg::HypotY(y.to_bits()),
            &mut fails,
            &format!("f64 hypot({x:e},{y:e})"),
        );
        // f32 sibling: the §9.2.1 worst case wants both widths.
        let xf = x as f32;
        let yf = y as f32;
        accept::<f32>(
            LibmFnId::Hypot,
            xf.to_bits(),
            LibmArg::HypotY(u64::from(yf.to_bits())),
            &mut fails,
            &format!("f32 hypot({xf:e},{yf:e})"),
        );
        if fails.len() > 50 {
            break;
        }
    }
    assert!(fails.is_empty(), "hypot failures: {fails:#?}");
}

#[test]
fn rootn_both_widths() {
    let n = sweep_size().min(120);
    let mut fails = Vec::new();
    let mut st = 0x726f_6f74_6e5f_6c6du64;
    for order in [-8i32, -5, -3, -2, 2, 3, 4, 5, 8] {
        for _ in 0..n {
            // Positive base across magnitudes (negative base + even
            // order is a domain error, covered by the smoke battery).
            let base = next_f64_banded(&mut st, -40, 40).abs();
            accept::<f64>(
                LibmFnId::Rootn(order),
                base.to_bits(),
                LibmArg::None,
                &mut fails,
                &format!("f64 rootn({base:e},{order})"),
            );
            let bf = base as f32;
            accept::<f32>(
                LibmFnId::Rootn(order),
                bf.to_bits(),
                LibmArg::None,
                &mut fails,
                &format!("f32 rootn({bf:e},{order})"),
            );
            if fails.len() > 50 {
                break;
            }
        }
    }
    assert!(fails.is_empty(), "rootn failures: {fails:#?}");
}
