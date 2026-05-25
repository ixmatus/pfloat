//! Per-function driver: drive `verify_input` across an input
//! iterator under each rounding mode, accumulate verdicts, capture
//! regression corpus, emit a `StatusRow`.
//!
//! The driver is single-threaded; sharding across cores is a
//! follow-up slice's concern. Per the Phase 1 plan: "Single-
//! threaded driver ships first; the shard coordinator and
//! rebalancer ... is a follow-up slice if release-CI runtime
//! motivates it."

#![cfg(all(unix, feature = "differential-mpfr"))]

use std::path::Path;

use pfloat::RoundingMode;

use super::status::{
    fnid_to_status_fields, DomainCoverage, PerModeStatus, RoundingStatus, StatusRow,
};
use super::types::{FnId, OracleBackend, Verdict};
use super::verify::{verify_input, Kernel};

/// Result of a per-function driver run before status-row emission.
/// The driver returns this so callers can decide where the
/// regression corpus and status row are written (the smoke gate
/// writes to a temp dir; the standalone runner writes to the
/// canonical `tests/vectors/` and `tests/oracle/status/` paths).
#[derive(Debug, Default)]
pub struct DriverOutcome {
    pub ok: u32,
    pub mismatch: Vec<(u32, RoundingMode, u32, u32)>,
    pub inconclusive: Vec<(u32, RoundingMode)>,
    pub panic: Vec<(u32, RoundingMode, String)>,
}

impl DriverOutcome {
    pub fn total(&self) -> u32 {
        self.ok
            + self.mismatch.len() as u32
            + self.inconclusive.len() as u32
            + self.panic.len() as u32
    }

    /// Compose the aggregate Phase 1 verdict from the verdict
    /// counts. `HasErrors` is any mismatch or panic; `Faithful`
    /// requires the worst observed error to be at most 1 ULP (the
    /// driver does not currently measure ULP and treats any
    /// mismatch as `HasErrors`, so a follow-up slice's
    /// per-mismatch ULP measurement is what unlocks the `Faithful`
    /// rung).
    ///
    /// The Phase 1f per-mode schema (ADR-0038) prefers
    /// [`Self::rounding_status_per_mode`] over this method;
    /// `rounding_status` is retained for back-compat with
    /// pre-migration callers (the smoke gate's assertions, the
    /// L-M corpus's success criterion).
    pub fn rounding_status(&self) -> RoundingStatus {
        if !self.mismatch.is_empty() || !self.panic.is_empty() {
            RoundingStatus::HasErrors
        } else {
            RoundingStatus::CorrectlyRounded
        }
    }

    /// Compose the Phase 1f per-mode verdict table. For each of
    /// the five IEEE 754-2019 rounding modes the driver swept
    /// (`swept_modes`), the corresponding entry records
    /// `CorrectlyRounded` when no mismatch or panic in this
    /// outcome was tagged with that mode, `HasErrors` when at
    /// least one was. Modes the driver did NOT sweep read
    /// `Unswept` (the Phase 1f transitional state; never the
    /// resting state at v1.0 per ADR-0038's no-narrowing
    /// principle).
    pub fn rounding_status_per_mode(&self, swept_modes: &[RoundingMode]) -> PerModeStatus {
        let mut status = PerModeStatus::all(RoundingStatus::Unswept);
        let entries: [(RoundingMode, &mut RoundingStatus); 5] = [
            (RoundingMode::NearestEven, &mut status.ne),
            (RoundingMode::NearestAway, &mut status.na),
            (RoundingMode::TowardZero, &mut status.tz),
            (RoundingMode::TowardPositive, &mut status.tp),
            (RoundingMode::TowardNegative, &mut status.tn),
        ];
        for (mode, slot) in entries {
            if !swept_modes.contains(&mode) {
                // Already initialized to Unswept.
                continue;
            }
            let any_failure = self
                .mismatch
                .iter()
                .any(|(_, m, _, _)| *m == mode)
                || self.panic.iter().any(|(_, m, _)| *m == mode);
            *slot = if any_failure {
                RoundingStatus::HasErrors
            } else {
                RoundingStatus::CorrectlyRounded
            };
        }
        status
    }
}

/// Drive `verify_input` for `f` across `inputs` under every mode in
/// `modes`. Catches panics from the pfloat kernel via
/// `std::panic::catch_unwind`; the catcher unwinds the
/// `verify_input` call directly so a kernel panic is captured at
/// the driver boundary rather than aborting the entire sweep.
pub fn run_function(
    oracle: &dyn OracleBackend,
    kernel: &Kernel<'_>,
    f: FnId,
    inputs: impl Iterator<Item = u32>,
    modes: &[RoundingMode],
) -> DriverOutcome {
    let mut out = DriverOutcome::default();
    for input in inputs {
        for &mode in modes {
            let verdict = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                verify_input(oracle, f, input, mode, kernel)
            }));
            match verdict {
                Ok(Verdict::Ok) => out.ok += 1,
                Ok(Verdict::Mismatch {
                    input,
                    mode,
                    expected,
                    got,
                }) => {
                    out.mismatch.push((input, mode, expected, got));
                }
                Ok(Verdict::OracleInconclusive { input, mode }) => {
                    out.inconclusive.push((input, mode));
                }
                Ok(Verdict::Panic { .. }) => {
                    // verify_input does not produce Verdict::Panic
                    // itself (the kernel can panic but the verifier
                    // does not synthesize it); this arm exists for
                    // future expansion.
                    unreachable!("verify_input does not emit Verdict::Panic")
                }
                Err(payload) => {
                    let msg = panic_message(payload.as_ref());
                    out.panic.push((input, mode, msg));
                }
            }
        }
    }
    out
}

/// Build a `StatusRow` from a `DriverOutcome` plus the per-run
/// metadata. `vectors_path` is the relative path where the
/// regression corpus was written, or empty when no corpus was
/// captured.
pub fn outcome_to_status_row(
    f: FnId,
    outcome: &DriverOutcome,
    domain_coverage: DomainCoverage,
    oracle_name: &'static str,
    modes: &[RoundingMode],
    vectors_path: &str,
) -> StatusRow {
    let (function, order) = fnid_to_status_fields(f);
    StatusRow {
        function,
        order,
        kernel_kind: "primary",
        domain_coverage,
        oracle: oracle_name,
        oracle_independence: "independent",
        rounding_status: outcome.rounding_status_per_mode(modes),
        worst_ulp: 0.0,
        mismatch_count: outcome.mismatch.len() as u32,
        inconclusive_count: outcome.inconclusive.len() as u32,
        panic_count: outcome.panic.len() as u32,
        vectors: vectors_path.to_string(),
        lm_seeds_run: 0,
    }
}

/// Serialize a `DriverOutcome`'s mismatch and panic entries to a
/// binary regression corpus file at `path`. The format is one
/// record per failure: `[input: u32 LE][mode_tag: u8][expected: u32
/// LE][got: u32 LE]` for mismatches; panics use a separate file
/// (the message is variable-length and benefits from human
/// inspection rather than binary regression replay).
pub fn write_mismatch_corpus(outcome: &DriverOutcome, path: &Path) -> std::io::Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    for &(input, mode, expected, got) in &outcome.mismatch {
        f.write_all(&input.to_le_bytes())?;
        f.write_all(&[mode_tag(mode)])?;
        f.write_all(&expected.to_le_bytes())?;
        f.write_all(&got.to_le_bytes())?;
    }
    Ok(())
}

fn mode_tag(m: RoundingMode) -> u8 {
    match m {
        RoundingMode::NearestEven => 0,
        RoundingMode::NearestAway => 1,
        RoundingMode::TowardZero => 2,
        RoundingMode::TowardPositive => 3,
        RoundingMode::TowardNegative => 4,
    }
}

fn panic_message(payload: &(dyn std::any::Any + Send)) -> String {
    if let Some(s) = payload.downcast_ref::<&'static str>() {
        (*s).to_string()
    } else if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else {
        "<non-string panic payload>".to_string()
    }
}
