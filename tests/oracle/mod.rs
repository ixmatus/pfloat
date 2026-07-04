//! Phase 1 Oracle harness: shared types and dispatch.
//!
//! Each `tests/oracle_*.rs` test crate (and the standalone
//! `examples/oracle_sweep.rs` runner, when it lands) does
//! `mod oracle;` to pull this module in. Cargo treats the file as a
//! shared submodule of every consumer (not as a standalone test
//! binary). The pattern mirrors `tests/differential/mod.rs`.
//!
//! ADR-0034 records the design. The oracle's correctness lift over
//! the existing differential lane is the enclosure posture: the
//! oracle returns a proven bracket `[lo, hi]` of the true value at
//! a requested working precision, not a rounded scalar. The
//! verifier asks whether both bracket endpoints round to the same
//! `f32` under the caller's mode; if they do, every point in the
//! bracket (including the true value) rounds there too and the
//! correctly-rounded `f32` is determined. If they straddle, the
//! working precision doubles and the bracket tightens; this is the
//! Ziv-at-oracle loop from ADR-0022 applied to the oracle's
//! evaluation rather than to a pfloat working-precision pass.

#![cfg(all(unix, feature = "differential-mpfr"))]
// Consumers use different subsets; suppress unused-warnings under
// any single crate's compilation. The `unused_imports` allow
// applies to the `pub use` re-exports below, which any one
// consumer may or may not reach into.
#![allow(dead_code)]
#![allow(unused_imports)]

#[cfg(feature = "differential-arb")]
pub mod arb;
pub mod convert;
#[cfg(all(feature = "differential-arb", feature = "ziv-instrumented"))]
pub mod cross_check;
pub mod driver;
#[cfg(feature = "differential-arb")]
pub mod maxima;
pub mod meta;
pub mod mpfr;
#[cfg(feature = "differential-arb")]
pub mod mpmath;
pub mod pfloat_kernels;
pub mod status;
pub mod types;
pub mod verify;

#[cfg(feature = "differential-arb")]
pub use arb::{ArbError, ArbOracle};
pub use convert::{bf24_of_bits, bf_to_f32_bits, certified_round_bf_to_f32, round_f32};
pub use driver::{outcome_to_status_row, run_function, write_mismatch_corpus, DriverOutcome};
#[cfg(feature = "differential-arb")]
pub use maxima::MaximaOracle;
pub use meta::{is_arb_primary, oracle_name_for, MetaError, MetaOracle};
pub use mpfr::MpfrOracle;
#[cfg(feature = "differential-arb")]
pub use mpmath::MpmathOracle;
pub use pfloat_kernels::pfloat_kernel;
pub use status::{DomainCoverage, PerModeStatus, RoundingStatus, StatusRow};
pub use types::{Enclosed, Enclosure, FnId, OracleBackend, Verdict};
pub use verify::{certified_round_f32, verify_input, Kernel, MAX_PREC, START_PREC};
