//! Status table schema and TOML emission.
//!
//! One `StatusRow` per `(function, rounding mode set)` pair the
//! sweep runs. The schema mirrors ADR-0034 and the Phase 1 plan's
//! "Status table output schema" section verbatim. Emission is by a
//! hand-written TOML writer; the row's fields are flat scalars and
//! enums, so pulling in serde + toml as test-only deps is not
//! warranted under the frugality posture.

#![cfg(all(unix, feature = "differential-mpfr"))]

use core::fmt::Write;

use pfloat::RoundingMode;

use super::types::FnId;

/// How thoroughly the sweep covered the function's binary32 input
/// space.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DomainCoverage {
    /// Every binary32 bit pattern was verified (`Exhaustive`) or a
    /// subset of that many inputs (`Sampled(n)`).
    Exhaustive,
    Sampled(u32),
}

/// The Phase 1 verdict for the (function, modes) pair the row
/// represents. ADR-0033 makes `HasErrors` a v1.0 blocker.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RoundingStatus {
    CorrectlyRounded,
    Faithful,
    HasErrors,
}

/// Per-(function, mode-set) sweep result; one row per pair lands in
/// `tests/oracle/status/<fn>.toml`. The README's per-function table
/// at v1.0 publishes from this schema directly.
#[derive(Clone, Debug)]
pub struct StatusRow {
    pub function: &'static str,
    /// Empty `""` for non-Bessel rows; the integer order otherwise.
    pub order: String,
    /// `"primary"` for direct kernels; `"derived_alias"` for
    /// composition wrappers (none today per ADR-0032).
    pub kernel_kind: &'static str,
    pub domain_coverage: DomainCoverage,
    /// Backend's `name()`: `"MPFR"`, `"Arb"`, or `"Arb+table"`.
    pub oracle: &'static str,
    /// Whether the oracle is independent of pfloat's algorithm
    /// class. `"independent"` for MPFR/Arb where the underlying
    /// algorithms differ; `"shared_algorithm_class"` flagged when
    /// the oracle and pfloat happen to share a series or recurrence
    /// shape.
    pub oracle_independence: &'static str,
    pub rounding_modes: Vec<RoundingMode>,
    pub rounding_status: RoundingStatus,
    /// Worst observed error in ULP across the (sub)swept space.
    /// `0.0` for `CorrectlyRounded` rows.
    pub worst_ulp: f64,
    pub mismatch_count: u32,
    pub inconclusive_count: u32,
    pub panic_count: u32,
    /// Path (relative to the repo root) of the regression corpus
    /// file capturing mismatch and panic inputs. Empty `""` when
    /// the file does not exist (no mismatches captured).
    pub vectors: String,
}

impl StatusRow {
    /// Hand-emit a TOML representation. The schema is flat enough
    /// that a serde+toml dep would not earn its keep.
    pub fn to_toml(&self) -> String {
        let mut out = String::new();
        writeln!(out, "function           = \"{}\"", self.function).unwrap();
        writeln!(out, "order              = \"{}\"", self.order).unwrap();
        writeln!(out, "kernel_kind        = \"{}\"", self.kernel_kind).unwrap();
        match self.domain_coverage {
            DomainCoverage::Exhaustive => {
                writeln!(out, "domain_coverage    = \"exhaustive\"").unwrap();
            }
            DomainCoverage::Sampled(n) => {
                writeln!(out, "domain_coverage    = \"sampled({n})\"").unwrap();
            }
        }
        writeln!(out, "oracle             = \"{}\"", self.oracle).unwrap();
        writeln!(
            out,
            "oracle_independence = \"{}\"",
            self.oracle_independence
        )
        .unwrap();
        let modes: Vec<&'static str> = self
            .rounding_modes
            .iter()
            .map(|m| match m {
                RoundingMode::NearestEven => "RNE",
                RoundingMode::NearestAway => "RNA",
                RoundingMode::TowardZero => "RZ",
                RoundingMode::TowardPositive => "RP",
                RoundingMode::TowardNegative => "RM",
            })
            .collect();
        writeln!(out, "rounding_modes     = \"{}\"", modes.join(",")).unwrap();
        let status = match self.rounding_status {
            RoundingStatus::CorrectlyRounded => "correctly-rounded",
            RoundingStatus::Faithful => "faithful",
            RoundingStatus::HasErrors => "has-errors",
        };
        writeln!(out, "rounding_status    = \"{status}\"").unwrap();
        writeln!(out, "worst_ulp          = {}", self.worst_ulp).unwrap();
        writeln!(out, "mismatch_count     = {}", self.mismatch_count).unwrap();
        writeln!(out, "inconclusive_count = {}", self.inconclusive_count).unwrap();
        writeln!(out, "panic_count        = {}", self.panic_count).unwrap();
        writeln!(out, "vectors            = \"{}\"", self.vectors).unwrap();
        out
    }
}

/// Format an `FnId` for the `function` and `order` columns. The
/// non-parametric variants emit just the function name; the
/// parametric Bessel variants emit the family name and the order
/// separately.
pub fn fnid_to_status_fields(f: FnId) -> (&'static str, String) {
    match f {
        FnId::BesselJn(n) | FnId::BesselYn(n) | FnId::BesselIn(n) | FnId::BesselKn(n) => {
            (f.name(), n.to_string())
        }
        _ => (f.name(), String::new()),
    }
}
