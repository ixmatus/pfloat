//! Standalone Phase 1 Oracle sweep runner.
//!
//! Invoke this binary by hand or via release CI to run the
//! exhaustive (or sampled) `f32` sweep across pfloat's MPFR-primary
//! surface. The same `verify_input` + `run_function` machinery the
//! per-push smoke gate uses; this binary just drives it with a
//! larger input budget and emits one TOML status row plus a binary
//! regression corpus per function.
//!
//! Build (requires `differential-mpfr`):
//!
//!     cargo build --release --features differential-mpfr \
//!         --example oracle_sweep
//!
//! Run examples:
//!
//!     # Sweep one function at 2^16 = 65536 sampled inputs under NE.
//!     ./target/release/examples/oracle_sweep --function sqrt \
//!         --sample 65536 --mode RNE
//!
//!     # Exhaustive 2^32 sweep of one function (long; minutes to
//!     # hours per the kernel cost).
//!     ./target/release/examples/oracle_sweep --function sqrt \
//!         --exhaustive --mode RNE
//!
//!     # Sweep every MPFR-primary function at 2^20 inputs under NE.
//!     ./target/release/examples/oracle_sweep --sample 1048576
//!
//! Output. One TOML status row per function in
//! `tests/oracle/status/<fn>.toml`; a binary regression corpus per
//! function with any mismatches in `tests/vectors/<fn>_regression.bin`.
//! Per ADR-0034 these are the v1.0 status table artifact and the
//! regression-replay store.

#![cfg(all(unix, feature = "differential-mpfr"))]

// The runner depends on the same Oracle harness modules as the
// integration tests. Include them via the same `mod oracle;`
// pattern (the harness lives under tests/oracle/ so an example
// target picks it up by a relative path include).

#[path = "../tests/oracle/mod.rs"]
mod oracle;

/// L-M hard-to-round corpus, included via relative path so the
/// runner can prepend per-function adversarial seeds without
/// re-exporting through the oracle harness module tree.
#[path = "../tests/differential/lefevre_muller_data.rs"]
mod lm_data;

use std::path::PathBuf;
use std::process::ExitCode;

use oracle::{
    outcome_to_status_row, pfloat_kernel, run_function, write_mismatch_corpus, DomainCoverage,
    FnId, Kernel, MpfrOracle, RoundingStatus,
};
use pfloat::RoundingMode;

/// MPFR-primary surface the runner covers by default. Mirrors the
/// smoke gate's list; the runner also accepts Bessel parametric
/// orders via `--function Jn:7` style syntax (next slice).
const MPFR_PRIMARY_FNIDS: &[FnId] = &[
    FnId::Sqrt,
    FnId::Exp,
    FnId::Exp2,
    FnId::Exp10,
    FnId::Expm1,
    FnId::Ln,
    FnId::Log1p,
    FnId::Log2,
    FnId::Log10,
    FnId::Sin,
    FnId::Cos,
    FnId::Tan,
    FnId::Asin,
    FnId::Acos,
    FnId::Atan,
    FnId::Sinh,
    FnId::Cosh,
    FnId::Tanh,
    FnId::Asinh,
    FnId::Acosh,
    FnId::Atanh,
    FnId::Erf,
    FnId::Erfc,
    FnId::Gamma,
    FnId::Lgamma,
    FnId::Digamma,
    FnId::Zeta,
    FnId::Ei,
    FnId::Ai,
    FnId::BesselJ0,
    FnId::BesselJ1,
    FnId::BesselY0,
    FnId::BesselY1,
];

struct Args {
    function: Option<FnId>,
    exhaustive: bool,
    sample: u32,
    modes: Vec<RoundingMode>,
    output_status: PathBuf,
    output_vectors: PathBuf,
}

impl Args {
    fn parse() -> Result<Self, String> {
        let mut function: Option<FnId> = None;
        let mut exhaustive = false;
        let mut sample: u32 = 1 << 20; // 1M default
        let mut modes: Vec<RoundingMode> = vec![RoundingMode::NearestEven];
        let mut output_status = PathBuf::from("tests/oracle/status");
        let mut output_vectors = PathBuf::from("tests/vectors");

        let mut args = std::env::args().skip(1);
        while let Some(a) = args.next() {
            match a.as_str() {
                "--function" => {
                    let name = args.next().ok_or("--function needs an argument")?;
                    function = Some(parse_fn_id(&name)?);
                }
                "--exhaustive" => exhaustive = true,
                "--sample" => {
                    let n = args.next().ok_or("--sample needs an argument")?;
                    sample = n.parse().map_err(|_| format!("bad --sample {n}"))?;
                }
                "--mode" => {
                    let s = args.next().ok_or("--mode needs an argument")?;
                    modes = s
                        .split(',')
                        .map(parse_mode)
                        .collect::<Result<Vec<_>, _>>()?;
                }
                "--output-status" => {
                    output_status =
                        PathBuf::from(args.next().ok_or("--output-status needs an argument")?);
                }
                "--output-vectors" => {
                    output_vectors =
                        PathBuf::from(args.next().ok_or("--output-vectors needs an argument")?);
                }
                "--help" | "-h" => {
                    print_help();
                    std::process::exit(0);
                }
                other => return Err(format!("unknown flag {other}")),
            }
        }

        Ok(Args {
            function,
            exhaustive,
            sample,
            modes,
            output_status,
            output_vectors,
        })
    }
}

fn parse_fn_id(name: &str) -> Result<FnId, String> {
    // Accept the canonical FnId::name() values plus the parametric
    // Bessel forms ("Jn:5" etc.).
    if let Some((family, n_str)) = name.split_once(':') {
        let n: i32 = n_str.parse().map_err(|_| format!("bad order in {name}"))?;
        return match family {
            "Jn" => Ok(FnId::BesselJn(n)),
            "Yn" => Ok(FnId::BesselYn(n)),
            "In" => Ok(FnId::BesselIn(n)),
            "Kn" => Ok(FnId::BesselKn(n)),
            other => Err(format!("unknown parametric family {other}")),
        };
    }
    for &f in MPFR_PRIMARY_FNIDS {
        if f.name() == name {
            return Ok(f);
        }
    }
    Err(format!("unknown --function {name}"))
}

fn parse_mode(s: &str) -> Result<RoundingMode, String> {
    match s {
        "RNE" | "NE" => Ok(RoundingMode::NearestEven),
        "RNA" | "NA" => Ok(RoundingMode::NearestAway),
        "RZ" => Ok(RoundingMode::TowardZero),
        "RP" => Ok(RoundingMode::TowardPositive),
        "RM" => Ok(RoundingMode::TowardNegative),
        other => Err(format!("unknown mode {other}")),
    }
}

fn print_help() {
    eprintln!(
        "Phase 1 Oracle sweep runner

Usage:
    oracle_sweep [OPTIONS]

Options:
    --function NAME   One of: sqrt, exp, exp2, exp10, expm1, ln,
                      log1p, log2, log10, sin, cos, tan, asin,
                      acos, atan, sinh, cosh, tanh, asinh, acosh,
                      atanh, erf, erfc, gamma, lgamma, digamma,
                      zeta, Ei, Ai, J0, J1, Y0, Y1, or parametric
                      forms Jn:N / Yn:N / In:N / Kn:N. Default:
                      every MPFR-primary function.
    --exhaustive      Sweep all 2^32 f32 bit patterns. Implies the
                      DomainCoverage::Exhaustive status row.
    --sample N        Sample N consecutive f32 bit patterns
                      starting at 0. Default: 2^20 (1048576).
    --mode MODES      Comma-separated rounding modes (RNE, RNA,
                      RZ, RP, RM). Default: RNE.
    --output-status DIR  Write per-function TOML status rows here.
                         Default: tests/oracle/status.
    --output-vectors DIR Write per-function regression corpus here.
                         Default: tests/vectors.
    --help            Print this message."
    );
}

fn run() -> Result<u32, String> {
    let args = Args::parse().inspect_err(|_e| {
        print_help();
    })?;

    std::fs::create_dir_all(&args.output_status)
        .map_err(|e| format!("create {}: {e}", args.output_status.display()))?;
    std::fs::create_dir_all(&args.output_vectors)
        .map_err(|e| format!("create {}: {e}", args.output_vectors.display()))?;

    let oracle = MpfrOracle;
    let kernel: &Kernel = &pfloat_kernel;
    let functions: Vec<FnId> = match args.function {
        Some(f) => vec![f],
        None => MPFR_PRIMARY_FNIDS.to_vec(),
    };

    let (input_iter_count, domain_coverage) = if args.exhaustive {
        (u32::MAX, DomainCoverage::Exhaustive)
    } else {
        (args.sample, DomainCoverage::Sampled(args.sample))
    };

    let mut has_errors_count: u32 = 0;
    for f in functions {
        let lm_seeds = lm_seeds_for(f);
        let lm_seeds_run = lm_seeds.len() as u32;
        let inputs = lm_seeds
            .into_iter()
            .chain((0u32..input_iter_count).take(input_iter_count as usize));
        eprint!("[oracle_sweep] {} ", f.name());
        let start = std::time::Instant::now();
        let outcome = run_function(&oracle, kernel, f, inputs, &args.modes);
        let elapsed = start.elapsed();
        eprintln!(
            "({}s): ok={}, mismatch={}, inconclusive={}, panic={}",
            elapsed.as_secs(),
            outcome.ok,
            outcome.mismatch.len(),
            outcome.inconclusive.len(),
            outcome.panic.len()
        );

        let vectors_path = if outcome.mismatch.is_empty() {
            String::new()
        } else {
            let rel = format!("tests/vectors/{}_regression.bin", f.name());
            let abs = args
                .output_vectors
                .join(format!("{}_regression.bin", f.name()));
            write_mismatch_corpus(&outcome, &abs)
                .map_err(|e| format!("write {}: {e}", abs.display()))?;
            rel
        };
        let mut row = outcome_to_status_row(
            f,
            &outcome,
            domain_coverage,
            "MPFR",
            &args.modes,
            &vectors_path,
        );
        row.lm_seeds_run = lm_seeds_run;
        let status_path = args.output_status.join(format!("{}.toml", row_filename(f)));
        std::fs::write(&status_path, row.to_toml())
            .map_err(|e| format!("write {}: {e}", status_path.display()))?;

        if row.rounding_status == RoundingStatus::HasErrors {
            has_errors_count += 1;
        }
    }

    eprintln!("[oracle_sweep] done. has-errors functions: {has_errors_count}");
    Ok(has_errors_count)
}

fn row_filename(f: FnId) -> String {
    match f {
        FnId::BesselJn(n) | FnId::BesselYn(n) | FnId::BesselIn(n) | FnId::BesselKn(n) => {
            format!("{}_{n}", f.name())
        }
        _ => f.name().to_string(),
    }
}

/// Per-`FnId` Lefèvre-Muller hard-to-round inputs, cast from
/// binary64 to binary32 bit patterns and deduplicated. The L-M
/// corpus covers 24 elementary and special-function entries (no
/// Bessel / Airy / Ei / Si / Ci / Li / digamma / zeta); other
/// `FnId`s return an empty vector. Inputs that map to `f32` NaN or
/// `f32` ±∞ after the cast are dropped so the linear sweep stays
/// the only carrier of those classes.
fn lm_seeds_for(f: FnId) -> Vec<u32> {
    use lm_data::{
        ACOSH_CASES, ACOS_CASES, ASINH_CASES, ASIN_CASES, ATANH_CASES, ATAN_CASES, COSH_CASES,
        COS_CASES, ERFC_CASES, ERF_CASES, EXP10_CASES, EXP2_CASES, EXPM1_CASES, EXP_CASES,
        GAMMA_CASES, LGAMMA_CASES, LN_CASES, LOG10_CASES, LOG1P_CASES, LOG2_CASES, SINH_CASES,
        SIN_CASES, TANH_CASES, TAN_CASES,
    };
    let cases: &[(u64, u64)] = match f {
        FnId::Exp => EXP_CASES,
        FnId::Ln => LN_CASES,
        FnId::Sin => SIN_CASES,
        FnId::Cos => COS_CASES,
        FnId::Tan => TAN_CASES,
        FnId::Atan => ATAN_CASES,
        FnId::Asin => ASIN_CASES,
        FnId::Acos => ACOS_CASES,
        FnId::Exp2 => EXP2_CASES,
        FnId::Exp10 => EXP10_CASES,
        FnId::Expm1 => EXPM1_CASES,
        FnId::Log1p => LOG1P_CASES,
        FnId::Log2 => LOG2_CASES,
        FnId::Log10 => LOG10_CASES,
        FnId::Sinh => SINH_CASES,
        FnId::Cosh => COSH_CASES,
        FnId::Tanh => TANH_CASES,
        FnId::Asinh => ASINH_CASES,
        FnId::Acosh => ACOSH_CASES,
        FnId::Atanh => ATANH_CASES,
        FnId::Erf => ERF_CASES,
        FnId::Erfc => ERFC_CASES,
        FnId::Gamma => GAMMA_CASES,
        FnId::Lgamma => LGAMMA_CASES,
        _ => return Vec::new(),
    };
    let mut seeds: Vec<u32> = cases
        .iter()
        .filter_map(|(input_bits, _)| {
            let x_f32 = f64::from_bits(*input_bits) as f32;
            if x_f32.is_finite() {
                Some(x_f32.to_bits())
            } else {
                None
            }
        })
        .collect();
    seeds.sort_unstable();
    seeds.dedup();
    seeds
}

fn main() -> ExitCode {
    match run() {
        Ok(0) => ExitCode::SUCCESS,
        // Any has-errors function means a v1.0 blocker. Exit code
        // reflects the count so CI can gate on it.
        Ok(_n) => ExitCode::from(1),
        Err(msg) => {
            eprintln!("error: {msg}");
            ExitCode::from(2)
        }
    }
}
