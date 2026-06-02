//! The pfloat-libm correctness sweep runner.
//!
//! Sweeps one function over a (sharded) input range, certifying the
//! shell's output against the MPFR oracle. The headline mode is the
//! exhaustive `f32` unary sweep: split the 2^32 `binary32` inputs across
//! shards with `--shard-index`/`--shard-count` (the capability pfloat's
//! pf-hcz4 runner lacked, which sharded only one function per instance),
//! verify NearestEven over the full shard range, and verify the four
//! directed modes over a strided subsample. Binary functions (`hypot`,
//! `rootn`) cannot be exhausted, so they take `--sample`. Emits a status
//! TOML row (mirroring pfloat's schema) and a JSON sidecar the
//! `pf-lm3-aggregate.py` script merges across shards.
//!
//! Examples:
//!   libm_sweep --function sin --width f32 --sample 65536 --mode all
//!   libm_sweep --function exp --exhaustive --shard-index 0 --shard-count 16 \
//!              --directed-sample 1048576 --output-json /tmp/exp_0of16.json
//!   libm_sweep --function hypot --sample 65536 --hypot-partner 0x3f000000

#![cfg(all(unix, feature = "differential-mpfr"))]

#[path = "../tests/harness/mod.rs"]
mod harness;

use std::fmt::Write as _;
use std::process::ExitCode;
use std::time::Instant;

use harness::status::DomainCoverage;
use harness::{
    lm_seeds_for, outcome_to_status_row, run_function, write_mismatch_corpus, DriverOutcome, Hw,
    LibmArg, LibmFnId, StatusGate, Width,
};
use pfloat_libm::RoundingMode;

const NE: RoundingMode = RoundingMode::NearestEven;
const DIRECTED: [RoundingMode; 4] = [
    RoundingMode::NearestAway,
    RoundingMode::TowardZero,
    RoundingMode::TowardPositive,
    RoundingMode::TowardNegative,
];
const ALL5: [RoundingMode; 5] = [
    RoundingMode::NearestEven,
    RoundingMode::NearestAway,
    RoundingMode::TowardZero,
    RoundingMode::TowardPositive,
    RoundingMode::TowardNegative,
];

const TWO_POW_32: u64 = 1 << 32;

struct Args {
    function: String,
    width: Width,
    exhaustive: bool,
    sample: u64,
    shard_index: u64,
    shard_count: u64,
    shard_start: Option<u64>,
    modes: Vec<RoundingMode>,
    directed_sample: u64,
    gate: StatusGate,
    hypot_partner: Option<u64>,
    output_status: Option<String>,
    output_vectors: Option<String>,
    output_json: Option<String>,
    instance_type: String,
}

fn main() -> ExitCode {
    match run() {
        Ok(false) => ExitCode::SUCCESS,
        Ok(true) => ExitCode::from(1), // has-errors
        Err(msg) => {
            eprintln!("error: {msg}");
            ExitCode::from(2)
        }
    }
}

fn run() -> Result<bool, String> {
    let args = parse_args()?;
    let (f, is_binary) = resolve_function(&args.function)?;

    if args.exhaustive && (is_binary || args.width == Width::F64) {
        return Err("--exhaustive applies only to unary f32 (binary and f64 use --sample)".into());
    }

    match args.width {
        Width::F32 => sweep::<f32>(f, is_binary, &args),
        Width::F64 => sweep::<f64>(f, is_binary, &args),
    }
}

/// Per-mode tallies for the JSON sidecar.
#[derive(Default, Clone, Copy)]
struct ModeTally {
    value_mismatch: u32,
    flag_mismatch: u32,
    inconclusive: u32,
    panic: u32,
}

fn merge(a: &mut DriverOutcome, b: DriverOutcome) {
    a.ok += b.ok;
    a.value_mismatch.extend(b.value_mismatch);
    a.flag_mismatch.extend(b.flag_mismatch);
    a.inconclusive.extend(b.inconclusive);
    a.panic.extend(b.panic);
}

fn sweep<H: Hw>(f: LibmFnId, is_binary: bool, args: &Args) -> Result<bool, String> {
    let arg = if matches!(f, LibmFnId::Hypot) {
        let p = args.hypot_partner.unwrap_or(default_partner(args.width));
        LibmArg::HypotY(p)
    } else {
        LibmArg::None
    };

    // Resolve the swept input range [start, end) over the u64 counter.
    let (start, end, exhaustive) = if args.exhaustive {
        let (s, e) = harness::shard_range(args.shard_index, args.shard_count, TWO_POW_32);
        if s >= e {
            return Err(format!(
                "shard-index {} is empty for shard-count {}",
                args.shard_index, args.shard_count
            ));
        }
        (s, e, true)
    } else {
        let s = args.shard_start.unwrap_or(0);
        let span = match args.width {
            Width::F32 => TWO_POW_32,
            Width::F64 => u64::MAX,
        };
        let e = s.saturating_add(args.sample).min(span);
        (s, e, false)
    };

    // Lefevre-Muller seeds fold in on shard 0 only (avoid cross-shard
    // duplication); cast to this width.
    let on_first_shard = start == 0;
    let seeds: Vec<H::Bits> = if on_first_shard {
        lm_seeds_for(f)
            .iter()
            .map(|&(input, _)| H::seed_from_f64_bits(input))
            .collect()
    } else {
        Vec::new()
    };
    let lm_seeds_run = seeds.len() as u32;

    let t0 = Instant::now();

    // Build the outcome. Exhaustive: NE over the full range + directed
    // over a strided subsample. Sample: the requested modes over the
    // sample.
    let mut outcome;
    let swept_modes: &[RoundingMode];
    let ne_count: u64;
    if exhaustive {
        let ne_iter = seeds
            .iter()
            .copied()
            .chain((start..end).map(H::bits_from_u64));
        outcome = run_function::<H, _>(f, ne_iter, arg, &[NE], args.gate);
        ne_count = (end - start) + u64::from(lm_seeds_run);

        if args.directed_sample > 0 {
            let span = end - start;
            let stride = (span / args.directed_sample.max(1)).max(1);
            let dir_iter = seeds
                .iter()
                .copied()
                .chain((start..end).step_by(stride as usize).map(H::bits_from_u64));
            let dir = run_function::<H, _>(f, dir_iter, arg, &DIRECTED, args.gate);
            merge(&mut outcome, dir);
        }
        swept_modes = &ALL5;
    } else {
        let iter = seeds
            .iter()
            .copied()
            .chain((start..end).map(H::bits_from_u64));
        outcome = run_function::<H, _>(f, iter, arg, &args.modes, args.gate);
        ne_count = (end - start) + u64::from(lm_seeds_run);
        swept_modes = &args.modes;
    }

    let wall = t0.elapsed().as_secs_f64();
    let has_errors = outcome.has_errors();

    // Status TOML row.
    let domain_coverage = if exhaustive {
        DomainCoverage::Exhaustive
    } else {
        DomainCoverage::Sampled(ne_count.min(u64::from(u32::MAX)) as u32)
    };
    let vectors_path = match (&args.output_vectors, outcome.value_mismatch.is_empty()) {
        (Some(dir), false) => {
            let p = format!("{dir}/{}_regression.bin", file_stem(f));
            write_mismatch_corpus(&outcome, std::path::Path::new(&p))
                .map_err(|e| format!("writing corpus: {e}"))?;
            p
        }
        _ => String::new(),
    };
    let row = outcome_to_status_row(
        f,
        &outcome,
        domain_coverage,
        "MPFR",
        swept_modes,
        &vectors_path,
        lm_seeds_run,
    );
    let toml = row.to_toml();

    if let Some(dir) = &args.output_status {
        let p = format!("{dir}/{}.toml", file_stem(f));
        std::fs::create_dir_all(dir).map_err(|e| format!("mkdir {dir}: {e}"))?;
        std::fs::write(&p, &toml).map_err(|e| format!("writing {p}: {e}"))?;
    }

    // JSON sidecar for the aggregator.
    let json = build_json(f, args, start, end, ne_count, lm_seeds_run, &outcome, wall);
    if let Some(p) = &args.output_json {
        std::fs::write(p, &json).map_err(|e| format!("writing {p}: {e}"))?;
    }

    // Human summary to stderr; the TOML row to stdout.
    eprintln!(
        "[sweep] {} {} range=[{:#x},{:#x}) ne={} directed_sample={} \
         ok={} value_mismatch={} flag_mismatch={} inconclusive={} panic={} {:.2}s {}",
        file_stem(f),
        args.width.name(),
        start,
        end,
        ne_count,
        if exhaustive { args.directed_sample } else { 0 },
        outcome.ok,
        outcome.value_mismatch.len(),
        outcome.flag_mismatch.len(),
        outcome.inconclusive.len(),
        outcome.panic.len(),
        wall,
        if has_errors { "HAS-ERRORS" } else { "clean" },
    );
    print!("{toml}");
    let _ = is_binary;
    Ok(has_errors)
}

fn per_mode_tally(outcome: &DriverOutcome, mode: RoundingMode) -> ModeTally {
    ModeTally {
        value_mismatch: outcome
            .value_mismatch
            .iter()
            .filter(|(_, m, _, _)| *m == mode)
            .count() as u32,
        flag_mismatch: outcome
            .flag_mismatch
            .iter()
            .filter(|(_, m, _, _, _)| *m == mode)
            .count() as u32,
        inconclusive: outcome
            .inconclusive
            .iter()
            .filter(|(_, m)| *m == mode)
            .count() as u32,
        panic: outcome.panic.iter().filter(|(_, m, _)| *m == mode).count() as u32,
    }
}

#[allow(clippy::too_many_arguments)]
fn build_json(
    f: LibmFnId,
    args: &Args,
    start: u64,
    end: u64,
    ne_count: u64,
    lm_seeds_run: u32,
    outcome: &DriverOutcome,
    wall: f64,
) -> String {
    let git_sha = std::env::var("PFLOAT_GIT_SHA").unwrap_or_else(|_| "unknown".into());
    let order = match f {
        LibmFnId::Rootn(n) => n.to_string(),
        _ => String::new(),
    };
    let directed = if args.exhaustive {
        args.directed_sample
    } else {
        0
    };
    let mut s = String::new();
    let _ = writeln!(s, "{{");
    let _ = writeln!(s, "  \"schema_version\": 1,");
    let _ = writeln!(s, "  \"function\": \"{}\",", f.name());
    let _ = writeln!(s, "  \"order\": \"{order}\",");
    let _ = writeln!(s, "  \"width\": \"{}\",", args.width.name());
    let _ = writeln!(s, "  \"exhaustive\": {},", args.exhaustive);
    let _ = writeln!(s, "  \"shard_index\": {},", args.shard_index);
    let _ = writeln!(s, "  \"shard_count\": {},", args.shard_count);
    let _ = writeln!(s, "  \"range_start\": {start},");
    let _ = writeln!(s, "  \"range_end\": {end},");
    let _ = writeln!(s, "  \"ne_inputs\": {ne_count},");
    let _ = writeln!(s, "  \"directed_sample\": {directed},");
    let _ = writeln!(s, "  \"lm_seeds_run\": {lm_seeds_run},");
    let _ = writeln!(s, "  \"ok\": {},", outcome.ok);
    let _ = writeln!(s, "  \"value_mismatch\": {},", outcome.value_mismatch.len());
    let _ = writeln!(s, "  \"flag_mismatch\": {},", outcome.flag_mismatch.len());
    let _ = writeln!(s, "  \"inconclusive\": {},", outcome.inconclusive.len());
    let _ = writeln!(s, "  \"panic\": {},", outcome.panic.len());
    // Per-mode breakdown.
    let _ = write!(s, "  \"per_mode\": {{");
    for (i, &m) in ALL5.iter().enumerate() {
        let t = per_mode_tally(outcome, m);
        let _ = write!(
            s,
            "{}\n    \"{}\": {{\"value_mismatch\": {}, \"flag_mismatch\": {}, \"inconclusive\": {}, \"panic\": {}}}",
            if i == 0 { "" } else { "," },
            mode_name(m),
            t.value_mismatch,
            t.flag_mismatch,
            t.inconclusive,
            t.panic
        );
    }
    let _ = writeln!(s, "\n  }},");
    // A bounded sample of value mismatches for triage.
    let _ = write!(s, "  \"sample_value_mismatches\": [");
    for (i, &(input, mode, expected, got)) in outcome.value_mismatch.iter().take(16).enumerate() {
        let _ = write!(
            s,
            "{}\n    {{\"input\": \"{input:#018x}\", \"mode\": \"{}\", \"expected\": \"{expected:#018x}\", \"got\": \"{got:#018x}\"}}",
            if i == 0 { "" } else { "," },
            mode_name(mode)
        );
    }
    let close = if outcome.value_mismatch.is_empty() {
        ""
    } else {
        "\n  "
    };
    let _ = writeln!(s, "{close}],");
    let _ = writeln!(s, "  \"wall_clock_seconds\": {wall:.3},");
    let _ = writeln!(s, "  \"instance_type\": \"{}\",", args.instance_type);
    let _ = writeln!(s, "  \"git_sha\": \"{git_sha}\",");
    let _ = writeln!(
        s,
        "  \"pfloat_libm_version\": \"{}\"",
        env!("CARGO_PKG_VERSION")
    );
    let _ = writeln!(s, "}}");
    s
}

fn mode_name(m: RoundingMode) -> &'static str {
    match m {
        RoundingMode::NearestEven => "NE",
        RoundingMode::NearestAway => "NA",
        RoundingMode::TowardZero => "TZ",
        RoundingMode::TowardPositive => "TP",
        RoundingMode::TowardNegative => "TN",
    }
}

/// File stem for the status/corpus files: the function name, with
/// `rootn`'s order appended (e.g. `rootn_-2`).
fn file_stem(f: LibmFnId) -> String {
    match f {
        LibmFnId::Rootn(n) => format!("rootn_{n}"),
        _ => f.name().to_string(),
    }
}

fn default_partner(width: Width) -> u64 {
    match width {
        Width::F32 => u64::from(0.5f32.to_bits()),
        Width::F64 => 0.5f64.to_bits(),
    }
}

fn resolve_function(name: &str) -> Result<(LibmFnId, bool), String> {
    if let Some(rest) = name.strip_prefix("rootn:") {
        let n: i32 = rest
            .parse()
            .map_err(|_| format!("bad rootn order: {rest}"))?;
        return Ok((LibmFnId::Rootn(n), true));
    }
    let f = match name {
        "exp" => LibmFnId::Exp,
        "exp2" => LibmFnId::Exp2,
        "exp10" => LibmFnId::Exp10,
        "expm1" => LibmFnId::Expm1,
        "ln" => LibmFnId::Ln,
        "log2" => LibmFnId::Log2,
        "log10" => LibmFnId::Log10,
        "log1p" => LibmFnId::Log1p,
        "sqrt" => LibmFnId::Sqrt,
        "cbrt" => LibmFnId::Cbrt,
        "sin" => LibmFnId::Sin,
        "cos" => LibmFnId::Cos,
        "tan" => LibmFnId::Tan,
        "cot" => LibmFnId::Cot,
        "sec" => LibmFnId::Sec,
        "csc" => LibmFnId::Csc,
        "asin" => LibmFnId::Asin,
        "acos" => LibmFnId::Acos,
        "atan" => LibmFnId::Atan,
        "sinh" => LibmFnId::Sinh,
        "cosh" => LibmFnId::Cosh,
        "tanh" => LibmFnId::Tanh,
        "asinh" => LibmFnId::Asinh,
        "acosh" => LibmFnId::Acosh,
        "atanh" => LibmFnId::Atanh,
        "hypot" => return Ok((LibmFnId::Hypot, true)),
        "rootn" => return Err("rootn needs an order: use rootn:N".into()),
        other => return Err(format!("unknown function: {other}")),
    };
    Ok((f, false))
}

fn parse_modes(spec: &str) -> Result<Vec<RoundingMode>, String> {
    if spec == "all" {
        return Ok(ALL5.to_vec());
    }
    spec.split(',')
        .map(|t| match t.trim() {
            "RNE" | "NE" => Ok(RoundingMode::NearestEven),
            "RNA" | "NA" => Ok(RoundingMode::NearestAway),
            "RZ" | "TZ" => Ok(RoundingMode::TowardZero),
            "RP" | "TP" => Ok(RoundingMode::TowardPositive),
            "RM" | "TN" => Ok(RoundingMode::TowardNegative),
            other => Err(format!("unknown mode: {other}")),
        })
        .collect()
}

fn parse_u64(s: &str) -> Result<u64, String> {
    let s = s.trim();
    if let Some(hex) = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X")) {
        u64::from_str_radix(hex, 16).map_err(|e| format!("bad hex {s}: {e}"))
    } else {
        s.parse::<u64>().map_err(|e| format!("bad number {s}: {e}"))
    }
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        function: String::new(),
        width: Width::F32,
        exhaustive: false,
        sample: 1 << 20,
        shard_index: 0,
        shard_count: 1,
        shard_start: None,
        modes: ALL5.to_vec(),
        directed_sample: 1 << 20,
        gate: StatusGate::ValueAndDomainHard,
        hypot_partner: None,
        output_status: None,
        output_vectors: None,
        output_json: None,
        instance_type: String::new(),
    };
    let mut it = std::env::args().skip(1);
    while let Some(flag) = it.next() {
        let mut next = || it.next().ok_or_else(|| format!("{flag} needs a value"));
        match flag.as_str() {
            "--function" => a.function = next()?,
            "--width" => {
                a.width = match next()?.as_str() {
                    "f32" => Width::F32,
                    "f64" => Width::F64,
                    w => return Err(format!("bad width: {w}")),
                }
            }
            "--exhaustive" => a.exhaustive = true,
            "--sample" => a.sample = parse_u64(&next()?)?,
            "--shard-index" => a.shard_index = parse_u64(&next()?)?,
            "--shard-count" => a.shard_count = parse_u64(&next()?)?,
            "--shard-start" => a.shard_start = Some(parse_u64(&next()?)?),
            "--mode" => a.modes = parse_modes(&next()?)?,
            "--directed-sample" => a.directed_sample = parse_u64(&next()?)?,
            "--status-gate" => {
                a.gate = match next()?.as_str() {
                    "hard" => StatusGate::ValueAndDomainHard,
                    "value-only" => StatusGate::ValueOnly,
                    g => return Err(format!("bad status-gate: {g}")),
                }
            }
            "--hypot-partner" => a.hypot_partner = Some(parse_u64(&next()?)?),
            "--output-status" => a.output_status = Some(next()?),
            "--output-vectors" => a.output_vectors = Some(next()?),
            "--output-json" => a.output_json = Some(next()?),
            "--instance-type" => a.instance_type = next()?,
            "--help" | "-h" => return Err(USAGE.into()),
            other => return Err(format!("unknown flag: {other}\n{USAGE}")),
        }
    }
    if a.function.is_empty() {
        return Err(format!("--function is required\n{USAGE}"));
    }
    Ok(a)
}

const USAGE: &str = "\
libm_sweep --function NAME [options]

  --function NAME       a unary name (exp, ln, sin, ...), `hypot`, or `rootn:N`
  --width {f32|f64}     default f32 (exhaustive is f32-unary only)
  --exhaustive          sweep all 2^32 binary32 inputs (sharded)
  --sample N            sample N inputs from --shard-start (default 2^20)
  --shard-index K       this shard's index (default 0)
  --shard-count M       number of shards splitting 2^32 (default 1)
  --shard-start BITS    explicit sample start (default 0)
  --mode {all|RNE,...}  modes for the sample run (default all)
  --directed-sample N   directed-mode subsample size for --exhaustive (default 2^20)
  --status-gate {hard|value-only}   default hard (value+INVALID+DIV_BY_ZERO)
  --hypot-partner BITS  fixed y bits for hypot (default 0.5 at the width)
  --output-status DIR   write <fn>.toml here
  --output-vectors DIR  write <fn>_regression.bin here on value mismatch
  --output-json PATH    write the JSON shard sidecar here
  --instance-type S     informational tag recorded in the JSON";
