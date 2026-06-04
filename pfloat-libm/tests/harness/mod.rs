//! Shared MPFR verification harness for pfloat-libm.
//!
//! Each consumer (`tests/libm_*.rs`, `examples/libm_sweep.rs`) pulls
//! this in via `mod harness;`. Cargo compiles it as a shared submodule
//! of each consumer crate, not a standalone target. The harness
//! reimplements pfloat's oracle/differential patterns (the dependency
//! direction forbids importing pfloat's test code) over an MPFR-only
//! oracle (`rug`). See ADR-0058.

#![cfg(all(unix, feature = "differential-mpfr"))]
// Each consumer (test/example crate) pulls in the whole harness via
// `mod harness;` but uses a different subset of its items and
// re-exports, so unused-item and unused-re-export warnings are expected
// here and not a signal. Mirrors pfloat's `tests/differential/mod.rs`.
#![allow(dead_code, unused_imports)]

pub mod convert;
pub mod driver;
pub mod hw;
pub mod lefevre_muller_data;
pub mod lm;
pub mod oracle;
pub mod rng;
pub mod sharding;
pub mod status;
pub mod status_gate;
pub mod types;
pub mod verify;

pub use convert::{round_f32, round_f64};
pub use driver::{outcome_to_status_row, run_function, write_mismatch_corpus, DriverOutcome};
pub use hw::Hw;
pub use lm::{lm_seeds_for, Case};
pub use oracle::{enclose, ORACLE_NAME};
pub use rng::{next_f64_banded, next_i64_in, next_u64, sweep_size, ALL_MODES};
pub use sharding::shard_range;
pub use status::{fnid_to_status_fields, DomainCoverage, PerModeStatus, RoundingStatus, StatusRow};
pub use status_gate::{check_flags, expected_div_by_zero, expected_invalid, StatusGate};
pub use types::{Enclosure, FlagKind, LibmArg, LibmFnId, Verdict, Width};
pub use verify::{certified_round, verify_input, MAX_PREC, START_PREC};
