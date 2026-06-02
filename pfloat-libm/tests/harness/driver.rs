//! Per-function driver: run [`verify_input`] across an input iterator
//! under each mode, accumulate verdicts, capture a regression corpus,
//! emit a [`StatusRow`].
//!
//! Width-generic over [`Hw`]; one driver serves the f32 exhaustive
//! sweep and the f64 differential. A function `HasErrors` iff it has a
//! value mismatch, a flag mismatch, or a panic. `inconclusive` (a
//! measure-zero hard-to-round straddle at the precision cap) does not
//! flip the verdict, mirroring pfloat's posture. Mirrors
//! `pfloat/tests/oracle/driver.rs`.

#![cfg(all(unix, feature = "differential-mpfr"))]

use std::path::Path;

use pfloat_libm::RoundingMode;

use super::hw::Hw;
use super::status::{
    fnid_to_status_fields, DomainCoverage, PerModeStatus, RoundingStatus, StatusRow,
};
use super::status_gate::StatusGate;
use super::types::{FlagKind, LibmArg, LibmFnId, Verdict};
use super::verify::verify_input;

/// Accumulated per-function result before status-row emission.
#[derive(Debug, Default)]
pub struct DriverOutcome {
    pub ok: u32,
    /// `(input, mode, expected_bits, got_bits)`.
    pub value_mismatch: Vec<(u64, RoundingMode, u64, u64)>,
    /// `(input, mode, flag, expected, got)`.
    pub flag_mismatch: Vec<(u64, RoundingMode, FlagKind, bool, bool)>,
    pub inconclusive: Vec<(u64, RoundingMode)>,
    pub panic: Vec<(u64, RoundingMode, String)>,
}

impl DriverOutcome {
    pub fn total(&self) -> u64 {
        u64::from(self.ok)
            + self.value_mismatch.len() as u64
            + self.flag_mismatch.len() as u64
            + self.inconclusive.len() as u64
            + self.panic.len() as u64
    }

    /// Any value mismatch, flag mismatch, or panic is `HasErrors`.
    pub fn has_errors(&self) -> bool {
        !self.value_mismatch.is_empty() || !self.flag_mismatch.is_empty() || !self.panic.is_empty()
    }

    /// Aggregate verdict (over all modes). `CorrectlyRounded` when clean.
    pub fn rounding_status(&self) -> RoundingStatus {
        if self.has_errors() {
            RoundingStatus::HasErrors
        } else {
            RoundingStatus::CorrectlyRounded
        }
    }

    /// Per-mode verdict table: `CorrectlyRounded` for a swept mode with
    /// no failure tagged to it, `HasErrors` otherwise; modes not in
    /// `swept_modes` read `Unswept`.
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
                continue;
            }
            let any_failure = self.value_mismatch.iter().any(|(_, m, _, _)| *m == mode)
                || self.flag_mismatch.iter().any(|(_, m, _, _, _)| *m == mode)
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

/// Drive `verify_input` for `f` across `inputs` (the swept axis) with a
/// fixed `arg` (the binary partner; [`LibmArg::None`] for unary) under
/// every mode in `modes`. Kernel panics are caught at the driver
/// boundary.
pub fn run_function<H, I>(
    f: LibmFnId,
    inputs: I,
    arg: LibmArg,
    modes: &[RoundingMode],
    gate: StatusGate,
) -> DriverOutcome
where
    H: Hw,
    I: Iterator<Item = H::Bits>,
{
    let mut out = DriverOutcome::default();
    for input in inputs {
        for &mode in modes {
            let verdict = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                verify_input::<H>(f, input, arg, mode, gate)
            }));
            match verdict {
                Ok(Verdict::Ok) => out.ok += 1,
                Ok(Verdict::ValueMismatch {
                    input,
                    mode,
                    expected,
                    got,
                }) => out.value_mismatch.push((input, mode, expected, got)),
                Ok(Verdict::FlagMismatch {
                    input,
                    mode,
                    flag,
                    expected,
                    got,
                }) => out.flag_mismatch.push((input, mode, flag, expected, got)),
                Ok(Verdict::OracleInconclusive { input, mode }) => {
                    out.inconclusive.push((input, mode))
                }
                Err(payload) => {
                    out.panic
                        .push((H::bits_to_u64(input), mode, panic_message(payload.as_ref())))
                }
            }
        }
    }
    out
}

/// Build a `StatusRow` from an outcome plus run metadata.
#[allow(clippy::too_many_arguments)]
pub fn outcome_to_status_row(
    f: LibmFnId,
    outcome: &DriverOutcome,
    domain_coverage: DomainCoverage,
    oracle_name: &'static str,
    modes: &[RoundingMode],
    vectors_path: &str,
    lm_seeds_run: u32,
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
        mismatch_count: (outcome.value_mismatch.len() + outcome.flag_mismatch.len()) as u32,
        inconclusive_count: outcome.inconclusive.len() as u32,
        panic_count: outcome.panic.len() as u32,
        vectors: vectors_path.to_string(),
        lm_seeds_run,
    }
}

/// Serialize value mismatches to a binary regression corpus:
/// `[input u64 LE][mode u8][expected u64 LE][got u64 LE]` per record
/// (f32 bits occupy the low 32 bits). Flag mismatches and panics are
/// human-inspected, not replayed, so they are not written here.
pub fn write_mismatch_corpus(outcome: &DriverOutcome, path: &Path) -> std::io::Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;
    for &(input, mode, expected, got) in &outcome.value_mismatch {
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
