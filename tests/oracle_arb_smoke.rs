//! Opt-in smoke gate for the Arb backend. Runs the ten
//! non-parametric Arb-primary `FnId`s through the Arb oracle
//! subprocess at a small input range per function (16 consecutive
//! f32 inputs anchored at each function's natural domain),
//! asserts zero mismatches, zero panics, zero inconclusive. NE
//! only. The parametric `BesselIn` / `BesselKn` are exercised
//! through the standalone runner (`--function In:5`,
//! `--function Kn:5`), matching the `BesselJn` / `BesselYn`
//! precedent the MPFR smoke gate sets.
//!
//! Not in the per-push CI lane. Compile and run with:
//!
//!     cargo test --features=differential-mpfr,differential-arb \
//!         --test oracle_arb_smoke
//!
//! Requires the Python venv from `scripts/setup_arb_oracle.sh`;
//! without it, `ArbOracle::new` returns
//! `ArbError::VenvNotFound` and the test fails with the setup
//! pointer in the message. The per-release runner
//! (`examples/oracle_sweep.rs`) does the full f32 sweep at
//! `--sample 65536` and emits per-row status TOMLs; this smoke
//! gate is the cheaper signal for "is the worker still alive and
//! the per-FnId dispatch still routing".

#![cfg(all(unix, feature = "differential-arb"))]

#[path = "oracle/mod.rs"]
mod oracle;

use oracle::{pfloat_kernel, run_function, FnId, Kernel, MetaOracle, RoundingStatus};
use pfloat::RoundingMode;

const NE: RoundingMode = RoundingMode::NearestEven;

/// 16 consecutive f32 inputs anchored at each function's natural
/// domain. The anchor avoids f32 NaN / inf and any pole or
/// singularity (e.g. `li(1)` is `-inf`); each function's choice is
/// documented inline.
fn smoke_inputs(f: FnId) -> std::ops::Range<u32> {
    // Anchor selection: pick an f32 normal that exercises each
    // function in a non-trivial regime without crossing a
    // domain boundary or pole.
    let anchor = match f {
        // Bessel I / K and Bi / Ai_prime / Bi_prime: x = 1.0 is
        // well inside the small-argument regime where Maclaurin
        // and asymptotic both converge.
        FnId::BesselI0
        | FnId::BesselI1
        | FnId::BesselIn(_)
        | FnId::BesselK0
        | FnId::BesselK1
        | FnId::BesselKn(_)
        | FnId::Bi
        | FnId::AiPrime
        | FnId::BiPrime => 0x3f80_0000u32, // 1.0
        // Si and Ci: Si is odd with Si(0) = 0; Ci has a
        // logarithmic singularity at 0 but is finite for x > 0.
        // Anchor at 0.5 to stay away from the singularity.
        FnId::Si | FnId::Ci => 0x3f00_0000u32, // 0.5
        // li(x) has a pole at x = 1 (li(1) = -inf). Anchor at
        // 2.0 to stay on the x > 1 side.
        FnId::Li => 0x4000_0000u32, // 2.0
        // Non-Arb-primary FnIds should never reach this helper;
        // the test below iterates only the Arb-primary set.
        _ => panic!("smoke_inputs called with non-Arb-primary FnId: {f:?}"),
    };
    anchor..(anchor + 16)
}

/// The ten non-parametric Arb-primary `FnId`s. Mirrors
/// `ARB_PRIMARY_FNIDS` in
/// `examples/oracle_sweep.rs`; the smoke test stays self-contained
/// rather than reaching into the runner.
const ARB_PRIMARY_FNIDS: &[FnId] = &[
    FnId::Si,
    FnId::Ci,
    FnId::Li,
    FnId::Bi,
    FnId::AiPrime,
    FnId::BiPrime,
    FnId::BesselI0,
    FnId::BesselI1,
    FnId::BesselK0,
    FnId::BesselK1,
];

#[test]
fn arb_smoke_gate_all_arb_primary_functions_clean() {
    let oracle = MetaOracle::new()
        .expect("MetaOracle::new (Arb backend requires the python-flint venv; run scripts/setup_arb_oracle.sh)");
    let kernel: &Kernel = &pfloat_kernel;
    let mut failures: Vec<(FnId, RoundingStatus, usize, usize)> = Vec::new();
    let mut total_verdicts: u32 = 0;
    for &f in ARB_PRIMARY_FNIDS {
        let outcome = run_function(&oracle, kernel, f, smoke_inputs(f), &[NE]);
        total_verdicts += outcome.total();
        let status = outcome.rounding_status();
        let inconclusive = outcome.inconclusive.len();
        let panic_count = outcome.panic.len();
        if status != RoundingStatus::CorrectlyRounded || inconclusive > 0 || panic_count > 0 {
            failures.push((f, status, inconclusive, panic_count));
            for (i, &(input, mode, expected, got)) in outcome.mismatch.iter().take(3).enumerate() {
                eprintln!(
                    "[arb-smoke] {f:?} mismatch #{i}: input={input:#010x} \
                     mode={mode:?} expected={expected:#010x} got={got:#010x}"
                );
            }
            for (i, &(input, mode)) in outcome.inconclusive.iter().take(3).enumerate() {
                eprintln!("[arb-smoke] {f:?} inconclusive #{i}: input={input:#010x} mode={mode:?}");
            }
            for (i, (input, mode, msg)) in outcome.panic.iter().take(3).enumerate() {
                eprintln!(
                    "[arb-smoke] {f:?} panic #{i}: input={input:#010x} mode={mode:?} msg={msg}"
                );
            }
        }
    }
    eprintln!(
        "[arb-smoke] swept {total_verdicts} verdicts across {} functions; {} failures",
        ARB_PRIMARY_FNIDS.len(),
        failures.len()
    );
    assert!(failures.is_empty(), "Arb smoke gate failures: {failures:?}");
}
