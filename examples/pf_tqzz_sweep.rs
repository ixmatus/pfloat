//! pf-tqzz full per-release cross-check sweep (pf-hcz4, ADR-0049).
//!
//! Runs the assertion machinery in [`oracle::cross_check`] across
//! one `FnId` shard at full input × mode granularity. Defaults to
//! 65 536 inputs × 5 IEEE rounding modes = 327 680 assertions per
//! invocation; the 47-FnId / 63-status-TOML surface gets the full
//! 15.4M-assertion sweep when each shard runs once.
//!
//! Failure handling: violations DO NOT panic. Each violating triple
//! lands in a sidecar list; the binary runs to completion across
//! every input × mode for the given `FnId` and emits a single
//! result JSON. The aggregator (`scripts/pf-hcz4-aggregate.py`)
//! merges the 63 shard files into the per-release v1.0 baseline.
//!
//! Build (requires `differential-arb` and `ziv-instrumented`):
//!
//!     cargo build --release --features differential-arb,ziv-instrumented \
//!         --example pf_tqzz_sweep
//!
//! Run examples:
//!
//!     # Sweep one FnId at the default 65536 × 5 grid.
//!     ./target/release/examples/pf_tqzz_sweep --fn-id Exp \
//!         --output /tmp/pf_tqzz_Exp.json
//!
//!     # Parametric Bessel: order 5.
//!     ./target/release/examples/pf_tqzz_sweep --fn-id Yn:5 \
//!         --output /tmp/pf_tqzz_Yn_5.json
//!
//!     # Single mode (smoke; the EC2 shard runs --modes all).
//!     ./target/release/examples/pf_tqzz_sweep --fn-id Exp \
//!         --modes RNE --output /tmp/pf_tqzz_Exp_RNE.json
//!
//!     # Smaller sample for local sanity (under a minute).
//!     ./target/release/examples/pf_tqzz_sweep --fn-id Exp \
//!         --sample 1024 --output /tmp/pf_tqzz_Exp_smoke.json
//!
//! Output JSON schema is documented in
//! `docs/decisions/0049-pf-hcz4-full-cross-check-sweep.md`.

#![cfg(all(unix, feature = "differential-arb", feature = "ziv-instrumented"))]

#[path = "../tests/oracle/mod.rs"]
mod oracle;

#[path = "../tests/differential/lefevre_muller_data.rs"]
mod lm_data;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::ExitCode;

use oracle::arb::ArbOracle;
use oracle::cross_check::{cross_check_one, CheckOutcome, ViolationRecord};
use oracle::mpfr::MpfrOracle;
use oracle::types::FnId;
use pfloat::RoundingMode;

// ===== CLI =====

struct Args {
    fn_id: FnId,
    output: PathBuf,
    modes: Vec<RoundingMode>,
    sample: u32,
    skip_lm_seeds: bool,
    instance_type: Option<String>,
}

impl Args {
    fn parse() -> Result<Self, String> {
        let mut fn_id: Option<FnId> = None;
        let mut output: Option<PathBuf> = None;
        let mut modes: Vec<RoundingMode> = all_modes();
        let mut sample: u32 = 65_536;
        let mut skip_lm_seeds = false;
        let mut instance_type: Option<String> = None;

        let mut args = std::env::args().skip(1);
        while let Some(a) = args.next() {
            match a.as_str() {
                "--fn-id" => {
                    let name = args.next().ok_or("--fn-id needs an argument")?;
                    fn_id = Some(parse_fn_id(&name)?);
                }
                "--output" => {
                    let p = args.next().ok_or("--output needs an argument")?;
                    output = Some(PathBuf::from(p));
                }
                "--modes" => {
                    let s = args.next().ok_or("--modes needs an argument")?;
                    modes = parse_modes(&s)?;
                }
                "--sample" => {
                    let n = args.next().ok_or("--sample needs an argument")?;
                    sample = n.parse().map_err(|_| format!("bad --sample {n}"))?;
                }
                "--skip-lm-seeds" => skip_lm_seeds = true,
                "--instance-type" => {
                    instance_type = Some(args.next().ok_or("--instance-type needs an argument")?);
                }
                "--help" | "-h" => {
                    print_help();
                    std::process::exit(0);
                }
                other => return Err(format!("unknown flag {other}")),
            }
        }

        let fn_id = fn_id.ok_or("missing required --fn-id")?;
        let output = output
            .unwrap_or_else(|| PathBuf::from(format!("/tmp/pf_tqzz_{}.json", row_filename(fn_id))));

        Ok(Args {
            fn_id,
            output,
            modes,
            sample,
            skip_lm_seeds,
            instance_type,
        })
    }
}

fn parse_fn_id(name: &str) -> Result<FnId, String> {
    // Mirrors examples/oracle_sweep.rs::parse_fn_id; accept canonical
    // names plus the parametric Bessel forms "Jn:5", "Yn:5", "In:5",
    // "Kn:5".
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
    for &f in MPFR_PRIMARY_FNIDS.iter().chain(ARB_PRIMARY_FNIDS.iter()) {
        if f.name() == name {
            return Ok(f);
        }
    }
    Err(format!("unknown --fn-id {name}"))
}

fn parse_modes(s: &str) -> Result<Vec<RoundingMode>, String> {
    if s == "all" {
        return Ok(all_modes());
    }
    s.split(',')
        .map(|tok| match tok {
            "RNE" | "NE" => Ok(RoundingMode::NearestEven),
            "RNA" | "NA" => Ok(RoundingMode::NearestAway),
            "RZ" => Ok(RoundingMode::TowardZero),
            "RP" => Ok(RoundingMode::TowardPositive),
            "RM" => Ok(RoundingMode::TowardNegative),
            other => Err(format!("unknown mode {other}")),
        })
        .collect()
}

fn all_modes() -> Vec<RoundingMode> {
    vec![
        RoundingMode::NearestEven,
        RoundingMode::NearestAway,
        RoundingMode::TowardZero,
        RoundingMode::TowardPositive,
        RoundingMode::TowardNegative,
    ]
}

fn print_help() {
    eprintln!(
        "pf-tqzz full per-release cross-check sweep (pf-hcz4, ADR-0049).

Usage:
    pf_tqzz_sweep --fn-id NAME [OPTIONS]

Required:
    --fn-id NAME      One canonical FnId name (Exp, Ln, Sin, ...) or
                      parametric Bessel form (Jn:5, Yn:5, In:5, Kn:5).

Options:
    --output PATH     Result JSON output path. Default:
                      /tmp/pf_tqzz_<fn_id>.json.
    --modes MODES     `all` for every IEEE 754-2019 mode, or
                      comma-separated subset (RNE, RNA, RZ, RP, RM).
                      Default: all.
    --sample N        Sweep first N f32 bit patterns (0u32..N).
                      Default: 65536 (matches pf-hcz4 bead spec).
    --skip-lm-seeds   Skip the Lefèvre-Muller hard-to-round corpus
                      (per-FnId; off by default — LM seeds add stress
                      to the Ziv-bound check).
    --instance-type S Optional informational tag recorded in JSON
                      (e.g. `c8g.large`).
    --help            Print this message."
    );
}

// ===== FnId surface =====

/// MPFR-primary surface mirrors `examples/oracle_sweep.rs`; the
/// parametric Bessel families (Jn / Yn / In / Kn) come in via the
/// `--fn-id Jn:5` syntax so this list does not need to enumerate
/// individual orders.
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

const ARB_PRIMARY_FNIDS: &[FnId] = &[
    FnId::Si,
    FnId::Ci,
    FnId::Li,
    FnId::Bi,
    FnId::AiPrime,
    FnId::BiPrime,
    FnId::BesselI0,
    FnId::BesselI1,
    FnId::BesselK0,
    FnId::BesselK1,
];

fn row_filename(f: FnId) -> String {
    match f {
        FnId::BesselJn(n) | FnId::BesselYn(n) | FnId::BesselIn(n) | FnId::BesselKn(n) => {
            format!("{}_{n}", f.name())
        }
        _ => f.name().to_string(),
    }
}

// ===== Per-mode result accumulator =====

#[derive(Default, Debug, Clone)]
struct ModeStats {
    passes: u64,
    skipped_no_ziv_path: u64,
    skipped_no_midpoint: u64,
    skipped_non_finite: u64,
    skipped_trace_not_final: u64,
    violations: u64,
}

impl ModeStats {
    fn record(&mut self, outcome: &CheckOutcome) {
        match outcome {
            CheckOutcome::Pass => self.passes += 1,
            CheckOutcome::SkippedNoZivPath => self.skipped_no_ziv_path += 1,
            CheckOutcome::SkippedNoMidpoint => self.skipped_no_midpoint += 1,
            CheckOutcome::SkippedNonFiniteMidpoint => self.skipped_non_finite += 1,
            CheckOutcome::SkippedTraceNotFinal => self.skipped_trace_not_final += 1,
            CheckOutcome::Violation(_) => self.violations += 1,
        }
    }
}

// ===== LM seeds =====

fn lm_seeds_for(f: FnId) -> Vec<u32> {
    use lm_data::{
        ACOSH_CASES, ACOS_CASES, ASINH_CASES, ASIN_CASES, ATANH_CASES, ATAN_CASES, COSH_CASES,
        COS_CASES, ERFC_CASES, ERF_CASES, EXP10_CASES, EXP2_CASES, EXPM1_CASES, EXP_CASES,
        GAMMA_CASES, LGAMMA_CASES, LN_CASES, LOG10_CASES, LOG1P_CASES, LOG2_CASES, SINH_CASES,
        SIN_CASES, TANH_CASES, TAN_CASES,
    };
    let cases: &[(u64, u64)] = match f {
        FnId::Exp => EXP_CASES,
        FnId::Exp2 => EXP2_CASES,
        FnId::Exp10 => EXP10_CASES,
        FnId::Expm1 => EXPM1_CASES,
        FnId::Ln => LN_CASES,
        FnId::Log1p => LOG1P_CASES,
        FnId::Log2 => LOG2_CASES,
        FnId::Log10 => LOG10_CASES,
        FnId::Sin => SIN_CASES,
        FnId::Cos => COS_CASES,
        FnId::Tan => TAN_CASES,
        FnId::Asin => ASIN_CASES,
        FnId::Acos => ACOS_CASES,
        FnId::Atan => ATAN_CASES,
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
    // Cast binary64 → binary32 bit pattern. Drop entries that fall
    // on NaN / ±∞ after the cast; the linear sweep covers those.
    let mut out: Vec<u32> = cases
        .iter()
        .map(|&(input, _)| f32::from_bits((input >> 32) as u32))
        .filter(|x| x.is_finite())
        .map(f32::to_bits)
        .collect();
    out.sort_unstable();
    out.dedup();
    out
}

// ===== Sweep loop (single-threaded per mode) =====

fn sweep_one_mode(
    fn_id: FnId,
    mode: RoundingMode,
    inputs: &[u32],
    arb: Option<&ArbOracle>,
    mpfr: &MpfrOracle,
) -> (ModeStats, Vec<ViolationRecord>) {
    let mut stats = ModeStats::default();
    let mut violations: Vec<ViolationRecord> = Vec::new();
    for &input in inputs {
        let outcome = cross_check_one(fn_id, input, mode, arb, mpfr);
        stats.record(&outcome);
        if let CheckOutcome::Violation(v) = outcome {
            violations.push(v);
        }
    }
    (stats, violations)
}

// ===== JSON emit =====

fn write_result_json(
    path: &std::path::Path,
    fn_id: FnId,
    sample: u32,
    used_lm_seeds: usize,
    modes: &[RoundingMode],
    per_mode: &[(RoundingMode, ModeStats, Vec<ViolationRecord>)],
    wall_clock_seconds: f64,
    instance_type: Option<&str>,
) -> std::io::Result<()> {
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);

    let order = match fn_id {
        FnId::BesselJn(n) | FnId::BesselYn(n) | FnId::BesselIn(n) | FnId::BesselKn(n) => Some(n),
        _ => None,
    };

    let mut totals = ModeStats::default();
    let mut arb_calls: u64 = 0;
    let mut mpfr_calls: u64 = 0;
    let is_arb = matches!(
        fn_id,
        FnId::Si
            | FnId::Ci
            | FnId::Li
            | FnId::Bi
            | FnId::AiPrime
            | FnId::BiPrime
            | FnId::BesselI0
            | FnId::BesselI1
            | FnId::BesselIn(_)
            | FnId::BesselK0
            | FnId::BesselK1
            | FnId::BesselKn(_)
    );
    for (_m, s, _v) in per_mode {
        let attempted = s.passes
            + s.skipped_no_ziv_path
            + s.skipped_no_midpoint
            + s.skipped_non_finite
            + s.skipped_trace_not_final
            + s.violations;
        // A midpoint is fetched only after the no-ziv-path and
        // trace-not-final guards pass; subtract both to count actual
        // backend MIDPOINT calls.
        let midpoint_calls = attempted - s.skipped_no_ziv_path - s.skipped_trace_not_final;
        if is_arb {
            arb_calls += midpoint_calls;
        } else {
            mpfr_calls += midpoint_calls;
        }
        totals.passes += s.passes;
        totals.skipped_no_ziv_path += s.skipped_no_ziv_path;
        totals.skipped_no_midpoint += s.skipped_no_midpoint;
        totals.skipped_non_finite += s.skipped_non_finite;
        totals.skipped_trace_not_final += s.skipped_trace_not_final;
        totals.violations += s.violations;
    }

    let error_guard = oracle::cross_check::error_guard_for(fn_id);
    let git_sha = std::env::var("PFLOAT_GIT_SHA")
        .unwrap_or_else(|_| read_git_sha().unwrap_or_else(|| "unknown".to_string()));
    let instance_arch = std::env::consts::ARCH;
    let instance_type_str = instance_type.unwrap_or("local");
    let pfloat_version = env!("CARGO_PKG_VERSION");

    writeln!(w, "{{")?;
    writeln!(w, "  \"schema_version\": 1,")?;
    writeln!(w, "  \"fn_id\": \"{}\",", fn_id.name())?;
    if let Some(n) = order {
        writeln!(w, "  \"order\": {n},")?;
    } else {
        writeln!(w, "  \"order\": null,")?;
    }
    writeln!(w, "  \"sample\": {sample},")?;
    writeln!(w, "  \"lm_seeds_used\": {used_lm_seeds},")?;
    write!(w, "  \"modes\": [")?;
    for (i, m) in modes.iter().enumerate() {
        if i > 0 {
            write!(w, ", ")?;
        }
        write!(w, "\"{}\"", mode_name(*m))?;
    }
    writeln!(w, "],")?;

    writeln!(w, "  \"totals\": {{")?;
    writeln!(w, "    \"passes\": {},", totals.passes)?;
    writeln!(
        w,
        "    \"skipped_no_ziv_path\": {},",
        totals.skipped_no_ziv_path
    )?;
    writeln!(
        w,
        "    \"skipped_no_midpoint\": {},",
        totals.skipped_no_midpoint
    )?;
    writeln!(
        w,
        "    \"skipped_non_finite\": {},",
        totals.skipped_non_finite
    )?;
    writeln!(
        w,
        "    \"skipped_trace_not_final\": {},",
        totals.skipped_trace_not_final
    )?;
    writeln!(w, "    \"violations\": {}", totals.violations)?;
    writeln!(w, "  }},")?;

    writeln!(w, "  \"per_mode\": {{")?;
    for (i, (m, s, _v)) in per_mode.iter().enumerate() {
        let comma = if i + 1 < per_mode.len() { "," } else { "" };
        writeln!(w, "    \"{}\": {{", mode_name(*m))?;
        writeln!(w, "      \"passes\": {},", s.passes)?;
        writeln!(
            w,
            "      \"skipped_no_ziv_path\": {},",
            s.skipped_no_ziv_path
        )?;
        writeln!(
            w,
            "      \"skipped_no_midpoint\": {},",
            s.skipped_no_midpoint
        )?;
        writeln!(w, "      \"skipped_non_finite\": {},", s.skipped_non_finite)?;
        writeln!(
            w,
            "      \"skipped_trace_not_final\": {},",
            s.skipped_trace_not_final
        )?;
        writeln!(w, "      \"violations\": {}", s.violations)?;
        writeln!(w, "    }}{comma}")?;
    }
    writeln!(w, "  }},")?;

    writeln!(w, "  \"violations\": [")?;
    // Flatten per-mode violations, sort by ratio_log2 descending.
    let mut all_violations: Vec<&ViolationRecord> =
        per_mode.iter().flat_map(|(_, _, v)| v.iter()).collect();
    all_violations.sort_by(|a, b| {
        ratio_log2(b)
            .partial_cmp(&ratio_log2(a))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for (i, v) in all_violations.iter().enumerate() {
        let comma = if i + 1 < all_violations.len() {
            ","
        } else {
            ""
        };
        write_violation(&mut w, v)?;
        writeln!(w, "{comma}")?;
    }
    writeln!(w, "  ],")?;

    writeln!(w, "  \"error_guard_const\": {error_guard},")?;
    writeln!(w, "  \"wall_clock_seconds\": {wall_clock_seconds:.3},")?;
    writeln!(w, "  \"arb_midpoint_calls\": {arb_calls},")?;
    writeln!(w, "  \"mpfr_midpoint_calls\": {mpfr_calls},")?;
    writeln!(w, "  \"git_sha\": \"{git_sha}\",")?;
    writeln!(w, "  \"instance_arch\": \"{instance_arch}\",")?;
    writeln!(w, "  \"instance_type\": \"{instance_type_str}\",")?;
    writeln!(w, "  \"pfloat_version\": \"{pfloat_version}\"")?;
    writeln!(w, "}}")?;
    w.flush()?;
    Ok(())
}

fn write_violation(w: &mut BufWriter<File>, v: &ViolationRecord) -> std::io::Result<()> {
    let input_f32 = f32::from_bits(v.input);
    writeln!(w, "    {{")?;
    writeln!(w, "      \"fn_id\": \"{}\",", v.fn_id.name())?;
    writeln!(w, "      \"input_u32\": {},", v.input)?;
    writeln!(w, "      \"input_hex\": \"0x{:08x}\",", v.input)?;
    writeln!(w, "      \"input_f32_repr\": \"{input_f32}\",")?;
    writeln!(w, "      \"mode\": \"{}\",", mode_name(v.mode))?;
    writeln!(w, "      \"working_prec\": {},", v.working_prec)?;
    writeln!(w, "      \"error_guard_const\": {},", v.error_guard)?;
    writeln!(w, "      \"eval_w_str\": \"{}\",", v.eval_w)?;
    writeln!(w, "      \"midpoint_str\": \"{}\",", v.midpoint)?;
    writeln!(w, "      \"abs_diff_str\": \"{}\",", v.abs_diff)?;
    writeln!(w, "      \"bound_str\": \"{}\",", v.bound)?;
    writeln!(w, "      \"gap_str\": \"{}\",", v.gap)?;
    writeln!(w, "      \"ratio_log2\": {:.6}", ratio_log2(v))?;
    write!(w, "    }}")?;
    Ok(())
}

/// `log2(abs_diff / bound)`. Negative = pass by margin, positive =
/// violation severity. Capped at ±256 to keep JSON output bounded.
fn ratio_log2(v: &ViolationRecord) -> f64 {
    // bound > 0 by construction (|midpoint| > 0 reached this path);
    // if abs_diff > bound we are in the violation branch by the
    // caller's choice.
    let lo_abs = v.abs_diff.clone().abs().to_f64();
    let lo_bnd = v.bound.clone().abs().to_f64();
    if lo_bnd == 0.0 {
        return 256.0;
    }
    (lo_abs / lo_bnd).log2().clamp(-256.0, 256.0)
}

fn mode_name(m: RoundingMode) -> &'static str {
    match m {
        RoundingMode::NearestEven => "NearestEven",
        RoundingMode::NearestAway => "NearestAway",
        RoundingMode::TowardZero => "TowardZero",
        RoundingMode::TowardPositive => "TowardPositive",
        RoundingMode::TowardNegative => "TowardNegative",
    }
}

fn read_git_sha() -> Option<String> {
    let out = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()?;
    if out.status.success() {
        Some(String::from_utf8_lossy(&out.stdout).trim().to_string())
    } else {
        None
    }
}

// ===== main =====

fn run() -> Result<(), String> {
    let args = Args::parse().inspect_err(|_e| print_help())?;

    let mpfr = MpfrOracle;
    let arb = ArbOracle::new().ok();
    if arb.is_none() && oracle::meta::is_arb_primary(args.fn_id) {
        return Err(format!(
            "fn_id {} requires Arb but the venv is unavailable; run scripts/setup_arb_oracle.sh",
            args.fn_id.name()
        ));
    }

    // Build the input grid: LM seeds (deduplicated) + (0..sample) bit patterns.
    let lm = if args.skip_lm_seeds {
        Vec::new()
    } else {
        lm_seeds_for(args.fn_id)
    };
    let used_lm = lm.len();
    let mut inputs: Vec<u32> = lm;
    inputs.extend(0u32..args.sample);
    inputs.sort_unstable();
    inputs.dedup();

    eprintln!(
        "[pf_tqzz_sweep] fn_id={} inputs={} modes={} → {} assertions",
        args.fn_id.name(),
        inputs.len(),
        args.modes.len(),
        inputs.len() * args.modes.len()
    );

    let start = std::time::Instant::now();
    let mut per_mode: Vec<(RoundingMode, ModeStats, Vec<ViolationRecord>)> = Vec::new();
    for mode in &args.modes {
        let mode_start = std::time::Instant::now();
        let (stats, violations) = sweep_one_mode(args.fn_id, *mode, &inputs, arb.as_ref(), &mpfr);
        let elapsed = mode_start.elapsed();
        eprintln!(
            "[pf_tqzz_sweep] mode={} ({}s): passes={} skipped={} violations={}",
            mode_name(*mode),
            elapsed.as_secs(),
            stats.passes,
            stats.skipped_no_ziv_path
                + stats.skipped_no_midpoint
                + stats.skipped_non_finite
                + stats.skipped_trace_not_final,
            stats.violations,
        );
        per_mode.push((*mode, stats, violations));
    }
    let wall = start.elapsed().as_secs_f64();

    write_result_json(
        &args.output,
        args.fn_id,
        args.sample,
        used_lm,
        &args.modes,
        &per_mode,
        wall,
        args.instance_type.as_deref(),
    )
    .map_err(|e| format!("write {}: {e}", args.output.display()))?;

    let total_violations: u64 = per_mode.iter().map(|(_, s, _)| s.violations).sum();
    eprintln!(
        "[pf_tqzz_sweep] done. wall={wall:.1}s violations={total_violations} output={}",
        args.output.display()
    );
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
