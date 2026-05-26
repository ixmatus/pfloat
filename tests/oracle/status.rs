//! Status table schema and TOML emission.
//!
//! One `StatusRow` per function the sweep runs; the row records a
//! per-mode `RoundingStatus` verdict via the `PerModeStatus` table
//! emitted as the `[rounding_status]` TOML section at the end of
//! each row file. The schema mirrors ADR-0034 (the initial Phase 1
//! design) as refined by ADR-0038 (Phase 1f's per-mode migration);
//! the per-row file format documents the v2 schema verbatim.
//! Emission is by a hand-written TOML writer; the row's fields are
//! flat scalars and enums plus one nested table, so pulling in
//! serde + toml as test-only deps is not warranted under the
//! frugality posture.

#![cfg(all(unix, feature = "differential-mpfr"))]

use core::fmt::Write;

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

/// The Phase 1 verdict for the (function, mode) pair the row's
/// per-mode entry represents. ADR-0033 makes `HasErrors` a v1.0
/// blocker. ADR-0038 (Phase 1f) introduces the `Unswept`
/// transitional state: a per-mode entry is `Unswept` while the
/// per-mode sweep has not yet certified that mode. No row exits
/// Phase 1f reading `Unswept` for any mode; the v1.0 ship criterion
/// is `CorrectlyRounded` across all five modes uniformly.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RoundingStatus {
    CorrectlyRounded,
    Faithful,
    HasErrors,
    /// Phase-1f-internal: the per-mode sweep has not yet certified
    /// this mode for this function. Never the resting state at
    /// v1.0 per ADR-0038's no-narrowing principle.
    Unswept,
}

/// Per-mode `RoundingStatus` table for a status row. Emits as the
/// `[rounding_status]` TOML table at the end of each per-row file.
/// The schema was migrated from a single-verdict scalar at Phase 1f
/// slice p1.23 (ADR-0038); the row now records one verdict per IEEE
/// 754-2019 rounding mode rather than collapsing the five modes
/// into one.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PerModeStatus {
    pub ne: RoundingStatus,
    pub na: RoundingStatus,
    pub tz: RoundingStatus,
    pub tp: RoundingStatus,
    pub tn: RoundingStatus,
}

impl PerModeStatus {
    /// All five modes record the same verdict. Convenience for the
    /// pre-migration data (every existing TOML row's data carried
    /// the same verdict for every mode it covered) and for fully-
    /// migrated rows (post Phase 1f, every row reads
    /// `CorrectlyRounded` uniformly).
    pub const fn all(status: RoundingStatus) -> Self {
        Self {
            ne: status,
            na: status,
            tz: status,
            tp: status,
            tn: status,
        }
    }

    /// NE-only sweep: NE carries `status`; the four directed modes
    /// read `Unswept`. The transitional shape slice p1.23 emits
    /// after the schema migration but before the per-family slices
    /// sweep the directed modes.
    pub const fn ne_only(status: RoundingStatus) -> Self {
        Self {
            ne: status,
            na: RoundingStatus::Unswept,
            tz: RoundingStatus::Unswept,
            tp: RoundingStatus::Unswept,
            tn: RoundingStatus::Unswept,
        }
    }
}

/// Per-function sweep result; one row per function lands in
/// `tests/oracle/status/<fn>.toml`. The README's per-function table
/// at v1.0 publishes from this schema directly.
///
/// Schema v2 (Phase 1f slice p1.23, ADR-0038): a single row records
/// one `RoundingStatus` per IEEE 754-2019 rounding mode via the
/// `[rounding_status]` TOML table. The pre-migration schema had a
/// `rounding_modes` field listing the modes the row covered plus a
/// single `rounding_status` scalar; the new schema treats every row
/// as covering all five modes implicitly (each mode entry carries
/// its own verdict; `Unswept` records modes the sweep has not yet
/// certified).
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
    /// Per-mode verdict table. Emits as the `[rounding_status]`
    /// TOML section at the end of the row file. Phase 1f's
    /// no-narrowing principle (ADR-0038) bars `Unswept` from any
    /// row at v1.0; the per-family slices p1.24 through p1.34
    /// migrate each kernel's directed-mode entries from `Unswept`
    /// to `CorrectlyRounded`.
    pub rounding_status: PerModeStatus,
    /// Worst observed error in ULP across the (sub)swept space.
    /// `0.0` for fully `CorrectlyRounded` rows.
    pub worst_ulp: f64,
    pub mismatch_count: u32,
    pub inconclusive_count: u32,
    pub panic_count: u32,
    /// Path (relative to the repo root) of the regression corpus
    /// file capturing mismatch and panic inputs. Empty `""` when
    /// the file does not exist (no mismatches captured).
    pub vectors: String,
    /// Number of Lefèvre-Muller hard-to-round inputs the sweep
    /// runner prepended to the linear input range for this row.
    /// `0` when the runner is the smoke gate (L-M seeds are
    /// runner-only) or when this `FnId` is outside the L-M
    /// corpus's 24-function coverage. The field documents the
    /// verification posture: a non-zero value records that the
    /// adversarial-seed lane ran in addition to the linear sweep.
    pub lm_seeds_run: u32,
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
        writeln!(out, "worst_ulp          = {}", self.worst_ulp).unwrap();
        writeln!(out, "mismatch_count     = {}", self.mismatch_count).unwrap();
        writeln!(out, "inconclusive_count = {}", self.inconclusive_count).unwrap();
        writeln!(out, "panic_count        = {}", self.panic_count).unwrap();
        writeln!(out, "vectors            = \"{}\"", self.vectors).unwrap();
        writeln!(out, "lm_seeds_run       = {}", self.lm_seeds_run).unwrap();
        // [rounding_status] section LAST so the named-table scope
        // doesn't capture the row's scalar fields above. The five
        // keys map 1:1 to the IEEE 754-2019 modes; the values are
        // the same vocabulary the pre-migration scalar carried with
        // the addition of `unswept` for Phase 1f's transitional
        // state (ADR-0038).
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
