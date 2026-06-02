//! Status table schema and TOML emission.
//!
//! One [`StatusRow`] per function, written to a `<fn>.toml` file: the
//! libm analogue of `pfloat/tests/oracle/status/*.toml`, and the
//! credibility document the v0.1 README will publish from. The schema
//! mirrors pfloat's byte-for-byte so the two crates' tables are
//! directly comparable (a hand-written TOML writer; the row is flat
//! scalars plus one nested table, so serde+toml would not earn its
//! keep). See ADR-0058.

#![cfg(all(unix, feature = "differential-mpfr"))]

use core::fmt::Write;

use super::types::LibmFnId;

/// How thoroughly the sweep covered the input space.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DomainCoverage {
    Exhaustive,
    Sampled(u32),
}

/// The verdict for one (function, mode) pair. `HasErrors` is a release
/// blocker; `Unswept` records a mode the sweep did not certify.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RoundingStatus {
    CorrectlyRounded,
    Faithful,
    HasErrors,
    Unswept,
}

/// Per-mode verdict table; emits as the `[rounding_status]` TOML table.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PerModeStatus {
    pub ne: RoundingStatus,
    pub na: RoundingStatus,
    pub tz: RoundingStatus,
    pub tp: RoundingStatus,
    pub tn: RoundingStatus,
}

impl PerModeStatus {
    pub const fn all(status: RoundingStatus) -> Self {
        Self {
            ne: status,
            na: status,
            tz: status,
            tp: status,
            tn: status,
        }
    }
}

/// Per-function sweep result; one row per `<fn>.toml`.
#[derive(Clone, Debug)]
pub struct StatusRow {
    pub function: &'static str,
    /// Empty `""` except `rootn`, which emits its order.
    pub order: String,
    pub kernel_kind: &'static str,
    pub domain_coverage: DomainCoverage,
    pub oracle: &'static str,
    pub oracle_independence: &'static str,
    pub rounding_status: PerModeStatus,
    pub worst_ulp: f64,
    pub mismatch_count: u32,
    pub inconclusive_count: u32,
    pub panic_count: u32,
    pub vectors: String,
    pub lm_seeds_run: u32,
}

impl StatusRow {
    /// Hand-emit a TOML representation matching pfloat's schema.
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
        writeln!(out, "worst_ulp          = {}", self.worst_ulp).unwrap();
        writeln!(out, "mismatch_count     = {}", self.mismatch_count).unwrap();
        writeln!(out, "inconclusive_count = {}", self.inconclusive_count).unwrap();
        writeln!(out, "panic_count        = {}", self.panic_count).unwrap();
        writeln!(out, "vectors            = \"{}\"", self.vectors).unwrap();
        writeln!(out, "lm_seeds_run       = {}", self.lm_seeds_run).unwrap();
        // [rounding_status] LAST so the table scope does not capture the
        // scalar fields above.
        writeln!(out).unwrap();
        writeln!(out, "[rounding_status]").unwrap();
        writeln!(out, "NE = \"{}\"", status_str(self.rounding_status.ne)).unwrap();
        writeln!(out, "NA = \"{}\"", status_str(self.rounding_status.na)).unwrap();
        writeln!(out, "TZ = \"{}\"", status_str(self.rounding_status.tz)).unwrap();
        writeln!(out, "TP = \"{}\"", status_str(self.rounding_status.tp)).unwrap();
        writeln!(out, "TN = \"{}\"", status_str(self.rounding_status.tn)).unwrap();
        out
    }
}

fn status_str(status: RoundingStatus) -> &'static str {
    match status {
        RoundingStatus::CorrectlyRounded => "correctly-rounded",
        RoundingStatus::Faithful => "faithful",
        RoundingStatus::HasErrors => "has-errors",
        RoundingStatus::Unswept => "unswept",
    }
}

/// The `function` and `order` columns for a `LibmFnId`. Only `rootn`
/// carries an order.
pub fn fnid_to_status_fields(f: LibmFnId) -> (&'static str, String) {
    match f {
        LibmFnId::Rootn(n) => (f.name(), n.to_string()),
        _ => (f.name(), String::new()),
    }
}
