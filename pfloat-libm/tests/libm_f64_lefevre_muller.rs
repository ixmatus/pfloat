//! Lefevre-Muller hard-to-round adversarial lane (binary64).
//!
//! For each of the 20 corpus-covered functions this lane runs two
//! checks on every `(input_bits, ne_output_bits)` case:
//!
//! 1. The pre-pinned NE assertion: the shell's nearest-even result
//!    equals the corpus's mpmath-computed binary64 output bit-for-bit.
//!    This needs no live oracle (the expected value is pinned).
//! 2. The live MPFR cross-check under all five rounding modes via
//!    `verify_input`, which extends the NE-only corpus to the directed
//!    modes against an independent oracle.
//!
//! Mirrors pfloat's `differential_lefevre_muller.rs`.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod harness;

use harness::lm::COVERED;
use harness::{lm_seeds_for, verify_input, Hw, LibmArg, LibmFnId, StatusGate, Verdict, ALL_MODES};
use pfloat_libm::RoundingMode;

const GATE: StatusGate = StatusGate::ValueAndDomainHard;

/// Nearest-even shell result bits for `f(from_bits(input))` at f64.
fn shell_ne_bits(f: LibmFnId, input: u64) -> u64 {
    <f64 as Hw>::shell(f, input, LibmArg::None, RoundingMode::NearestEven).0
}

#[test]
fn lefevre_muller_corpus_pinned_and_cross_checked() {
    let mut pinned = 0u64;
    let mut crosschecked = 0u64;
    let mut failures: Vec<String> = Vec::new();

    for &f in COVERED {
        let cases = lm_seeds_for(f);
        assert!(
            !cases.is_empty(),
            "{} should have Lefevre-Muller cases",
            f.name()
        );
        for &(input, want_ne) in cases {
            // (1) Pinned NE assertion (NaN-aware).
            let got_ne = shell_ne_bits(f, input);
            let nan_match = f64::from_bits(got_ne).is_nan() && f64::from_bits(want_ne).is_nan();
            if !(nan_match || got_ne == want_ne) {
                failures.push(format!(
                    "{} NE pin: input={input:#018x} want={want_ne:#018x} got={got_ne:#018x}",
                    f.name()
                ));
            }
            pinned += 1;

            // (2) Live MPFR cross-check, all five modes.
            for &mode in ALL_MODES {
                let v = verify_input::<f64>(f, input, LibmArg::None, mode, GATE);
                if !matches!(v, Verdict::Ok | Verdict::OracleInconclusive { .. }) {
                    failures.push(format!(
                        "{} {mode:?} cross-check: input={input:#018x}: {v:?}",
                        f.name()
                    ));
                }
                crosschecked += 1;
            }
        }
    }

    eprintln!(
        "[lm] pinned {pinned} NE outputs, cross-checked {crosschecked} (fn x mode x case); \
         {} failures",
        failures.len()
    );
    assert!(
        failures.is_empty(),
        "Lefevre-Muller failures: {failures:#?}"
    );
}
