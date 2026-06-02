//! Fast per-push correctness gate for the MPFR harness.
//!
//! Sweeps a small strided neighbourhood (`N` consecutive ULPs from each
//! function's domain anchor) across all 25 unary functions at both
//! widths and the 2 binary functions with fixed partners, under all
//! five rounding modes, against the MPFR oracle. Asserts zero value
//! mismatches, zero `INVALID`/`DIV_BY_ZERO` mismatches, zero panics,
//! zero inconclusives. This is the lane the libm CI job runs on every
//! push; the exhaustive 2^32 grid is the EC2 deliverable
//! (`examples/libm_sweep.rs`). Mirrors pfloat's `oracle_smoke_gate.rs`,
//! widened to both widths and all modes.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod harness;

use harness::{run_function, DriverOutcome, LibmArg, LibmFnId, StatusGate, ALL_MODES};

/// Inputs per (function, width). Small enough that the whole gate runs
/// in seconds (release); dense enough to catch near-boundary rounding.
const N: u32 = 64;
const GATE: StatusGate = StatusGate::ValueAndDomainHard;

/// f32 domain anchor: the bit pattern the `N`-ULP sweep starts at. Most
/// functions anchor at `0.5` (in every domain, away from the reciprocal
/// poles); `acosh` needs `x >= 1`.
fn anchor_f32(f: LibmFnId) -> u32 {
    match f {
        LibmFnId::Acosh => 0x3fc0_0000, // 1.5
        _ => 0x3f00_0000,               // 0.5
    }
}

/// f64 domain anchor, the f64 analogue of [`anchor_f32`].
fn anchor_f64(f: LibmFnId) -> u64 {
    match f {
        LibmFnId::Acosh => 0x3ff8_0000_0000_0000, // 1.5
        _ => 0x3fe0_0000_0000_0000,               // 0.5
    }
}

/// Record a failing function for the final assertion message.
struct Failures(Vec<String>);

impl Failures {
    fn check(&mut self, label: &str, out: &DriverOutcome) {
        if out.has_errors() || !out.inconclusive.is_empty() {
            for (i, &(input, mode, expected, got)) in out.value_mismatch.iter().take(3).enumerate()
            {
                eprintln!(
                    "[smoke] {label} value #{i}: input={input:#018x} mode={mode:?} \
                     expected={expected:#018x} got={got:#018x}"
                );
            }
            for (i, &(input, mode, flag, expected, got)) in
                out.flag_mismatch.iter().take(3).enumerate()
            {
                eprintln!(
                    "[smoke] {label} flag #{i}: input={input:#018x} mode={mode:?} \
                     {} expected={expected} got={got}",
                    flag.name()
                );
            }
            for (i, &(input, mode)) in out.inconclusive.iter().take(3).enumerate() {
                eprintln!("[smoke] {label} inconclusive #{i}: input={input:#018x} mode={mode:?}");
            }
            for (i, (input, mode, msg)) in out.panic.iter().take(3).enumerate() {
                eprintln!("[smoke] {label} panic #{i}: input={input:#018x} mode={mode:?} {msg}");
            }
            self.0.push(format!(
                "{label} (value={} flag={} inconclusive={} panic={})",
                out.value_mismatch.len(),
                out.flag_mismatch.len(),
                out.inconclusive.len(),
                out.panic.len()
            ));
        }
    }
}

#[test]
fn smoke_gate_all_functions_clean() {
    let mut fails = Failures(Vec::new());
    let mut verdicts: u64 = 0;

    // Unary, both widths.
    for &f in LibmFnId::UNARY {
        let a32 = anchor_f32(f);
        let out = run_function::<f32, _>(f, a32..(a32 + N), LibmArg::None, ALL_MODES, GATE);
        verdicts += out.total();
        fails.check(&format!("f32 {}", f.name()), &out);

        let a64 = anchor_f64(f);
        let out =
            run_function::<f64, _>(f, a64..(a64 + u64::from(N)), LibmArg::None, ALL_MODES, GATE);
        verdicts += out.total();
        fails.check(&format!("f64 {}", f.name()), &out);
    }

    // hypot(x, 0.5) over the anchor neighbourhood, both widths.
    let a32 = 0x3f00_0000u32;
    let out = run_function::<f32, _>(
        LibmFnId::Hypot,
        a32..(a32 + N),
        LibmArg::HypotY(u64::from(0.5f32.to_bits())),
        ALL_MODES,
        GATE,
    );
    verdicts += out.total();
    fails.check("f32 hypot", &out);
    let a64 = 0x3fe0_0000_0000_0000u64;
    let out = run_function::<f64, _>(
        LibmFnId::Hypot,
        a64..(a64 + u64::from(N)),
        LibmArg::HypotY(0.5f64.to_bits()),
        ALL_MODES,
        GATE,
    );
    verdicts += out.total();
    fails.check("f64 hypot", &out);

    // rootn(x, n) over a few orders (positive base neighbourhood).
    for n in [2i32, 3, -2] {
        let out = run_function::<f32, _>(
            LibmFnId::Rootn(n),
            a32..(a32 + N),
            LibmArg::None,
            ALL_MODES,
            GATE,
        );
        verdicts += out.total();
        fails.check(&format!("f32 rootn:{n}"), &out);
        let out = run_function::<f64, _>(
            LibmFnId::Rootn(n),
            a64..(a64 + u64::from(N)),
            LibmArg::None,
            ALL_MODES,
            GATE,
        );
        verdicts += out.total();
        fails.check(&format!("f64 rootn:{n}"), &out);
    }

    eprintln!(
        "[smoke] swept {verdicts} verdicts; {} failing groups",
        fails.0.len()
    );
    assert!(fails.0.is_empty(), "smoke gate failures: {:?}", fails.0);
}
