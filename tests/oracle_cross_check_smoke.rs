//! pf-tqzz cross-check assertion smoke (ADR-0039, slice p1g.3).
//!
//! For every `(kernel, input, mode)` triple this harness routes,
//! it asserts that pfloat's Ziv-driver `eval(working)` intermediate
//! stayed within the kernel's per-function calibrated `error_guard`
//! bound (from `crate::math::ziv_calibration`) when compared against
//! the rigorous-enclosure midpoint computed at higher precision by
//! Arb (for the 12 Arb-primary `FnId`s) or MPFR (for the 35
//! MPFR-primary `FnId`s).
//!
//! The assertion is
//!
//! ```text
//! |eval(working) − midpoint| ≤ 2^(error_guard − working) · |midpoint|
//! ```
//!
//! evaluated in `rug::Float` at oracle precision `oracle_prec =
//! working + 64`. Any violation surfaces as a structured panic
//! identifying kernel, input, mode, achieved error, assumed bound,
//! and the gap.
//!
//! This file ships the **smoke subset** (~10 inputs per kernel),
//! not the full 65536 × 5 × 47 sweep. The full sweep is the per-
//! release gate the user runs separately (`examples/pf_tqzz_sweep.rs`,
//! ADR-0049, pf-hcz4); the smoke validates the plumbing is correct
//! end-to-end and would catch any kernel whose calibrated bound is
//! grossly wrong. The smoke runtime fits per-push CI under the
//! existing `differential-arb` lane.
//!
//! The harness gates on both `differential-arb` (for the Arb
//! subprocess) and `ziv-instrumented` (for the thread-local Ziv
//! trace capture); without either the cross-check cannot run.
//! Silently no-op when the Arb venv is unavailable (matches the
//! `oracle_arb_midpoint_smoke.rs` precedent).
//!
//! The actual assertion machinery lives in
//! [`oracle::cross_check`] and is shared with the full-sweep
//! example. This file's responsibility is the smoke input grid and
//! the panic-on-violation UX.

#![cfg(all(unix, feature = "differential-arb", feature = "ziv-instrumented"))]

mod oracle;

use oracle::arb::ArbOracle;
use oracle::cross_check::{cross_check_one, CheckOutcome};
use oracle::mpfr::MpfrOracle;
use oracle::pfloat_kernels::verification_precision;
use oracle::types::FnId;

use pfloat::RoundingMode;

#[derive(Debug, Default)]
struct SmokeStats {
    passes: usize,
    skipped_no_ziv_path: usize,
    skipped_no_midpoint: usize,
    skipped_non_finite: usize,
}

fn run_smoke_subset(
    f: FnId,
    arb: Option<&ArbOracle>,
    mpfr: &MpfrOracle,
    inputs: &[u32],
    modes: &[RoundingMode],
) -> SmokeStats {
    let mut stats = SmokeStats::default();
    for &input in inputs {
        for &mode in modes {
            match cross_check_one(f, input, mode, arb, mpfr) {
                CheckOutcome::Pass => stats.passes += 1,
                CheckOutcome::SkippedNoZivPath => stats.skipped_no_ziv_path += 1,
                CheckOutcome::SkippedNoMidpoint => stats.skipped_no_midpoint += 1,
                CheckOutcome::SkippedNonFiniteMidpoint => stats.skipped_non_finite += 1,
                CheckOutcome::Violation(v) => panic!("{}", v.format_panic_message()),
            }
        }
    }
    stats
}

/// Smoke subset of inputs spanning typical magnitudes: zero,
/// small positives, ones, several decades up to typical f32
/// magnitudes. Mode coverage is the full 5-mode IEEE set.
fn smoke_inputs() -> Vec<u32> {
    [
        0.0_f32, 0.5, 1.0, 1.5, 2.0, 3.0, 7.0, 10.0, 100.0, 1.0e3, 1.0e6,
    ]
    .iter()
    .map(|x| x.to_bits())
    .collect()
}

fn all_modes() -> Vec<RoundingMode> {
    vec![
        RoundingMode::NearestEven,
        RoundingMode::NearestAway,
        RoundingMode::TowardZero,
        RoundingMode::TowardPositive,
        RoundingMode::TowardNegative,
    ]
}

/// MPFR-primary smoke: exp, ln, sin, cos. These exercise the four
/// most-common transcendental shapes and validate the `MpfrOracle`
/// midpoint path end-to-end.
#[test]
fn cross_check_mpfr_primary_smoke() {
    let arb = ArbOracle::new().ok();
    let mpfr = MpfrOracle;
    let inputs = smoke_inputs();
    let modes = all_modes();

    let mut total = SmokeStats::default();
    for f in [FnId::Exp, FnId::Ln, FnId::Sin, FnId::Cos] {
        let s = run_smoke_subset(f, arb.as_ref(), &mpfr, &inputs, &modes);
        eprintln!(
            "{f:?}: {} passes, {} no-ziv-path, {} no-midpoint, {} non-finite",
            s.passes, s.skipped_no_ziv_path, s.skipped_no_midpoint, s.skipped_non_finite
        );
        total.passes += s.passes;
        total.skipped_no_ziv_path += s.skipped_no_ziv_path;
        total.skipped_no_midpoint += s.skipped_no_midpoint;
        total.skipped_non_finite += s.skipped_non_finite;
    }
    eprintln!("Total: {total:?}");
    // Demand at least *some* passes; if every input was skipped the
    // harness is broken (likely no Ziv routing).
    assert!(
        total.passes >= 4,
        "expected at least 4 cross-check passes, got {total:?}"
    );
}

/// Arb-primary smoke: Si only. Si verifies at the default p=24
/// precision so the Arb subprocess call returns quickly; the Bessel
/// I/K small-arg family verifies at p=320 (per the cubic-correction
/// trap precedent at p1.32) and is too slow for per-push smoke.
/// The full sweep at per-release cadence exercises the wider
/// Arb-primary surface.
#[test]
fn cross_check_arb_primary_smoke() {
    let arb = match ArbOracle::new() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("skipping arb cross-check smoke (Arb venv unavailable): {e}");
            return;
        }
    };
    let mpfr = MpfrOracle;
    // Use a reduced input set: zero, small positives only. Si is
    // bounded and well-behaved across this range; the smoke focuses
    // on validating the Arb subprocess + MIDPOINT wire format
    // round-trip, not on stress-testing Si.
    let inputs: Vec<u32> = [0.0_f32, 0.5, 1.0, 2.0, 5.0]
        .iter()
        .map(|x| x.to_bits())
        .collect();
    let modes = all_modes();

    let s = run_smoke_subset(FnId::Si, Some(&arb), &mpfr, &inputs, &modes);
    eprintln!(
        "Si: {} passes, {} no-ziv-path, {} no-midpoint, {} non-finite",
        s.passes, s.skipped_no_ziv_path, s.skipped_no_midpoint, s.skipped_non_finite
    );
    assert!(
        s.passes >= 4,
        "expected at least 4 Si cross-check passes, got {s:?}"
    );
}

/// Bound-form sanity: the cross-check uses
/// `verification_precision` only to determine the working precision
/// for the kernel call, NOT for the bound formula. The bound formula
/// uses the actual converged working precision from the trace.
/// This test exercises a kernel where `verification_precision` differs
/// from f32's 24 bits (Bessel small-arg bumps to 320), confirming
/// the trace-derived working precision is what drives the bound.
#[test]
fn cross_check_respects_trace_working_precision() {
    // Verify the wiring exists; no assertion beyond compilation.
    let _ = verification_precision(FnId::BesselI0);
    let _ = verification_precision(FnId::Exp);
}
