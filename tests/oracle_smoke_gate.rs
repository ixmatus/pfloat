//! Per-push smoke gate for the Phase 1 Oracle harness.
//!
//! Runs `verify_input` on every MPFR-primary `FnId` across a small
//! representative input set (64 inputs per function) under
//! `NearestEven`, asserting zero mismatches, zero panics, zero
//! oracle-inconclusives. This is the lane the
//! `differential-mpfr` CI job runs on every push to surface a
//! regression in either pfloat's kernel or the oracle harness
//! itself.
//!
//! Runtime budget. Debug-mode CI completes in around 20-30 seconds
//! across the 33 MPFR-primary functions; release-mode runs in
//! around 5 seconds. The exhaustive sweep ships through
//! `examples/oracle_sweep.rs` (next commit) for the per-release
//! lane that pays many CPU-minutes per function.
//!
//! The full-surface exhaustive sweep is a separate runner; see
//! `examples/oracle_sweep.rs` (next commit). Per-release CI invokes
//! that runner with the `--exhaustive` flag.
//!
//! Per-function domain. The smoke gate picks 64 f32 bit patterns
//! starting at the function's domain anchor (typically `0.5`, but
//! `acosh` needs `x >= 1` so it anchors at `1.5`). The inputs walk
//! consecutive f32 ULPs from there; the range is dense enough to
//! cover near-rounding-boundary cases but small enough that the
//! whole gate runs in seconds.
//!
//! Findings expectation. Slice p1.2 closed the five-finding
//! has-errors class on the v1.0 surface (the L-M corpus had
//! everything passing by slice p1.2.4). The smoke gate is biased
//! toward NORMAL-range inputs near `0.5` or `1.5`, so the L-M-style
//! hard-to-round cases at subnormal boundaries are not directly
//! exercised here; the standalone runner is where exhaustive sweeps
//! catch those. If a function unexpectedly has-errors at this
//! gate, the finding is captured in the driver's regression corpus,
//! a defect bead gets filed, and the kernel fix lands in a follow-up
//! slice unless cheap enough to absorb in-slice.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod oracle;

use oracle::{pfloat_kernel, run_function, FnId, Kernel, MpfrOracle, RoundingStatus};
use pfloat::RoundingMode;

const NE: RoundingMode = RoundingMode::NearestEven;

/// Per-function domain anchor (the bit pattern the 64-element ULP
/// sweep starts at). Most functions anchor at f32 `0.5`
/// (`0x3f000000`), which is in every transcendental function's
/// domain. Functions requiring `x >= 1` anchor at `1.5`.
fn domain_anchor(f: FnId) -> u32 {
    match f {
        FnId::Acosh => 0x3fc0_0000u32, // 1.5
        _ => 0x3f00_0000u32,           // 0.5
    }
}

/// 64 f32 ULPs starting at the function's anchor. The slice is
/// fixed-size so per-push runtime is predictable across the
/// 33-function set.
fn smoke_inputs(f: FnId) -> std::ops::Range<u32> {
    let anchor = domain_anchor(f);
    anchor..(anchor + 64)
}

const MPFR_PRIMARY_FNIDS: &[FnId] = &[
    // Elementary.
    FnId::Sqrt,
    FnId::Exp,
    FnId::Exp2,
    FnId::Exp10,
    FnId::Expm1,
    FnId::Ln,
    FnId::Log1p,
    FnId::Log2,
    FnId::Log10,
    FnId::Sin,
    FnId::Cos,
    FnId::Tan,
    FnId::Asin,
    FnId::Acos,
    FnId::Atan,
    FnId::Sinh,
    FnId::Cosh,
    FnId::Tanh,
    FnId::Asinh,
    FnId::Acosh,
    FnId::Atanh,
    // Specials with MPFR primitive.
    FnId::Erf,
    FnId::Erfc,
    FnId::Gamma,
    FnId::Lgamma,
    FnId::Digamma,
    FnId::Zeta,
    FnId::Ei,
    // Airy: only Ai has an MPFR primitive.
    FnId::Ai,
    // Bessel J / Y fixed-order; parametric Jn / Yn deferred to the
    // runner where the order list is configurable.
    FnId::BesselJ0,
    FnId::BesselJ1,
    FnId::BesselY0,
    FnId::BesselY1,
];

#[test]
fn smoke_gate_all_mpfr_primary_functions_clean() {
    let oracle = MpfrOracle;
    let kernel: &Kernel = &pfloat_kernel;
    let mut failures: Vec<(FnId, RoundingStatus, usize, usize)> = Vec::new();
    let mut total_verdicts: u32 = 0;
    for &f in MPFR_PRIMARY_FNIDS {
        let outcome = run_function(&oracle, kernel, f, smoke_inputs(f), &[NE]);
        total_verdicts += outcome.total();
        let status = outcome.rounding_status();
        let inconclusive = outcome.inconclusive.len();
        let panic_count = outcome.panic.len();
        if status != RoundingStatus::CorrectlyRounded || inconclusive > 0 || panic_count > 0 {
            failures.push((f, status, inconclusive, panic_count));
            // Also print the first few mismatches for diagnostic
            // context; the captured ones land in the regression
            // corpus when the standalone runner runs.
            for (i, &(input, mode, expected, got)) in outcome.mismatch.iter().take(3).enumerate() {
                eprintln!(
                    "[smoke] {f:?} mismatch #{i}: input={input:#010x} \
                     mode={mode:?} expected={expected:#010x} got={got:#010x}"
                );
            }
            for (i, &(input, mode)) in outcome.inconclusive.iter().take(3).enumerate() {
                eprintln!("[smoke] {f:?} inconclusive #{i}: input={input:#010x} mode={mode:?}");
            }
            for (i, (input, mode, msg)) in outcome.panic.iter().take(3).enumerate() {
                eprintln!("[smoke] {f:?} panic #{i}: input={input:#010x} mode={mode:?} msg={msg}");
            }
        }
    }
    eprintln!(
        "[smoke] swept {} verdicts across {} functions; {} failures",
        total_verdicts,
        MPFR_PRIMARY_FNIDS.len(),
        failures.len()
    );
    assert!(failures.is_empty(), "smoke gate failures: {failures:?}");
}
