//! Smoke tests for the per-function driver and status table
//! emitter.
//!
//! The driver runs pfloat's sqrt under the MPFR oracle across a
//! small dense sweep; the resulting `DriverOutcome` and
//! `StatusRow` shapes are asserted. A separate test verifies the
//! TOML emission matches the documented schema field-for-field.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod oracle;

use oracle::{
    outcome_to_status_row, pfloat_kernel, run_function, write_mismatch_corpus, DomainCoverage,
    FnId, Kernel, MpfrOracle, RoundingStatus, StatusRow,
};
use pfloat::RoundingMode;

const NE: RoundingMode = RoundingMode::NearestEven;

#[test]
fn run_function_sqrt_dense_1024_ok() {
    let oracle = MpfrOracle;
    let kernel: &Kernel = &pfloat_kernel;
    let inputs = 0x3f80_0000u32..0x3f80_0400; // 1024 f32 normals at 1.0+.
    let outcome = run_function(&oracle, kernel, FnId::Sqrt, inputs, &[NE]);
    assert_eq!(outcome.total(), 1024, "input cap was 1024 per mode");
    assert_eq!(outcome.ok, 1024);
    assert!(outcome.mismatch.is_empty());
    assert!(outcome.inconclusive.is_empty());
    assert!(outcome.panic.is_empty());
    assert_eq!(outcome.rounding_status(), RoundingStatus::CorrectlyRounded);
}

#[test]
fn run_function_sqrt_all_modes_ok() {
    let oracle = MpfrOracle;
    let kernel: &Kernel = &pfloat_kernel;
    let modes = [
        RoundingMode::NearestEven,
        RoundingMode::NearestAway,
        RoundingMode::TowardZero,
        RoundingMode::TowardPositive,
        RoundingMode::TowardNegative,
    ];
    let inputs = 0x3f80_0000u32..0x3f80_0040; // 64 inputs × 5 modes = 320 verdicts
    let outcome = run_function(&oracle, kernel, FnId::Sqrt, inputs, &modes);
    assert_eq!(outcome.total(), 64 * 5);
    assert_eq!(outcome.ok, 64 * 5);
}

#[test]
fn run_function_catches_panic_from_kernel() {
    let oracle = MpfrOracle;
    let kernel: &Kernel = &|_f, _input, _mode| panic!("synthetic kernel panic");
    let outcome = run_function(&oracle, kernel, FnId::Sqrt, 0u32..1u32, &[NE]);
    assert_eq!(outcome.total(), 1);
    assert_eq!(outcome.panic.len(), 1);
    let (input, _mode, msg) = &outcome.panic[0];
    assert_eq!(*input, 0);
    assert!(
        msg.contains("synthetic kernel panic"),
        "panic message {msg:?} did not contain expected substring"
    );
    assert_eq!(outcome.rounding_status(), RoundingStatus::HasErrors);
}

#[test]
fn status_row_toml_emission_matches_schema() {
    let row = StatusRow {
        function: "sqrt",
        order: String::new(),
        kernel_kind: "primary",
        domain_coverage: DomainCoverage::Sampled(1024),
        oracle: "MPFR",
        oracle_independence: "independent",
        rounding_modes: vec![NE],
        rounding_status: RoundingStatus::CorrectlyRounded,
        worst_ulp: 0.0,
        mismatch_count: 0,
        inconclusive_count: 0,
        panic_count: 0,
        vectors: String::new(),
    };
    let toml = row.to_toml();
    // Spot-check the schema fields named in ADR-0034.
    assert!(toml.contains("function           = \"sqrt\""));
    assert!(toml.contains("order              = \"\""));
    assert!(toml.contains("kernel_kind        = \"primary\""));
    assert!(toml.contains("domain_coverage    = \"sampled(1024)\""));
    assert!(toml.contains("oracle             = \"MPFR\""));
    assert!(toml.contains("oracle_independence = \"independent\""));
    assert!(toml.contains("rounding_modes     = \"RNE\""));
    assert!(toml.contains("rounding_status    = \"correctly-rounded\""));
    assert!(toml.contains("worst_ulp          = 0"));
    assert!(toml.contains("mismatch_count     = 0"));
    assert!(toml.contains("inconclusive_count = 0"));
    assert!(toml.contains("panic_count        = 0"));
    assert!(toml.contains("vectors            = \"\""));
}

#[test]
fn status_row_for_bessel_carries_order() {
    let row = outcome_to_status_row(
        FnId::BesselJn(7),
        &oracle::DriverOutcome::default(),
        DomainCoverage::Sampled(0),
        "MPFR",
        &[NE],
        "",
    );
    assert_eq!(row.function, "Jn");
    assert_eq!(row.order, "7");
}

#[test]
fn write_mismatch_corpus_emits_expected_byte_count() {
    let mut outcome = oracle::DriverOutcome::default();
    outcome
        .mismatch
        .push((0x3f80_0000, NE, 0x3f80_0000, 0x3f80_0001));
    outcome
        .mismatch
        .push((0x4000_0000, NE, 0x4000_0000, 0x3fff_ffff));
    // Record is 13 bytes: u32 input + u8 mode + u32 expected + u32 got.
    let tmp = std::env::temp_dir().join("pfloat-oracle-test-corpus.bin");
    write_mismatch_corpus(&outcome, &tmp).unwrap();
    let bytes = std::fs::read(&tmp).unwrap();
    assert_eq!(bytes.len(), 2 * 13);
    std::fs::remove_file(&tmp).ok();
}
