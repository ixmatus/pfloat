//! Arb oracle backend: the [`OracleBackend`] implementation for the
//! twelve `FnId`s the MPFR backend cannot cover (`Si`, `Ci`, `Li`,
//! `Bi`, `Ai_prime`, `Bi_prime`, `BesselI{0,1,n}`,
//! `BesselK{0,1,n}`).
//!
//! Evaluations happen out-of-process in a long-lived `python-flint`
//! worker; the worker reads one request per line on its stdin and
//! emits one certified `f32` bit pattern per line on its stdout.
//! This module owns the subprocess, the request / response protocol,
//! and the construction of the single-point [`Enclosure`] the
//! verifier rounds at the caller's mode.
//!
//! See ADR-0034 for the framing (two-backend layering, LGPL
//! isolation via subprocess) and ADR-0035 for the refined protocol
//! (worker reports certified f32 directly, no decimal bridge). The
//! slice p1.7 reclassification of pf-6a4e showed the decimal
//! bracket protocol was correctness-load-bearing in two ways that
//! were silently broken; ADR-0035 records the cure and slice p1.8
//! implements it.
//!
//! ## Venv resolution
//!
//! The Arb backend needs a Python venv with `python-flint`
//! installed. The default path is `${HOME}/.cache/pfloat-arb-oracle/venv`;
//! override via the `PFLOAT_ARB_ORACLE_VENV` env var. The
//! [`scripts/setup_arb_oracle.sh`](../../scripts/setup_arb_oracle.sh)
//! helper creates and verifies the venv idempotently.
//!
//! ## Worker protocol (ADR-0035)
//!
//! Request: `<fn_id> <order_or_dash> <input_bits_hex> <mode>`.
//!
//! Response: `OK <f32_bits_hex>` (the certified f32 bit pattern as
//! 8 lowercase hex chars), `INC` (the worker's internal Ziv loop
//! could not certify a unique f32 at its maximum precision), or
//! `ERR <message>` (an error processing the request).
//!
//! The worker runs the Ziv-at-oracle loop in-process (ball
//! arithmetic stays in binary, no decimal bridge); the
//! [`ArbOracle::enclose`] response is wrapped as a single-point
//! [`Enclosure`] at the certified `f32` (or NaN endpoints for
//! `INC`, which the verifier surfaces as `OracleInconclusive`).
//! [`ArbOracle::is_authoritative`] returns `true` so the verifier
//! short-circuits its outer Ziv loop and accepts the worker's
//! single answer.

#![cfg(all(unix, feature = "differential-arb"))]

use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};
use std::sync::Mutex;

use rug::float::Special;
use rug::Float;

use pfloat::RoundingMode;

use super::types::{Enclosure, FnId, OracleBackend};

/// Errors the Arb oracle can surface at construction time. Per-call
/// failures (worker died mid-request, malformed response) escalate
/// to panics inside [`OracleBackend::enclose`] because the trait
/// signature is infallible; the driver's `catch_unwind` records
/// those as `Verdict::Panic`.
#[derive(Debug)]
pub enum ArbError {
    /// The configured venv path does not contain a Python
    /// interpreter; carries the path that was checked.
    VenvNotFound(PathBuf),
    /// The worker script is missing or not readable; carries the
    /// path that was checked.
    WorkerScriptNotFound(PathBuf),
    /// Failure spawning the worker subprocess.
    Spawn(std::io::Error),
    /// The worker's response to the initial `ready?` ping was
    /// unexpected.
    HandshakeFailed(String),
}

impl std::fmt::Display for ArbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::VenvNotFound(p) => write!(
                f,
                "Arb oracle venv not found at {}; run scripts/setup_arb_oracle.sh \
                 (or set PFLOAT_ARB_ORACLE_VENV)",
                p.display()
            ),
            Self::WorkerScriptNotFound(p) => {
                write!(f, "Arb oracle worker script not found at {}", p.display())
            }
            Self::Spawn(e) => write!(f, "Failed to spawn Arb oracle worker: {e}"),
            Self::HandshakeFailed(msg) => write!(f, "Arb oracle worker handshake failed: {msg}"),
        }
    }
}

impl std::error::Error for ArbError {}

/// Resolve the venv path: env var override, otherwise
/// `${HOME}/.cache/pfloat-arb-oracle/venv`.
fn default_venv_path() -> PathBuf {
    if let Ok(p) = std::env::var("PFLOAT_ARB_ORACLE_VENV") {
        return PathBuf::from(p);
    }
    let home =
        std::env::var("HOME").expect("HOME must be set for the default Arb oracle venv path");
    PathBuf::from(home).join(".cache/pfloat-arb-oracle/venv")
}

/// Path to the in-tree worker script. Resolved from `CARGO_MANIFEST_DIR`
/// so the runner / smoke tests find it regardless of the cwd they
/// were launched from.
fn worker_script_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scripts/arb_oracle_worker.py")
}

/// Owned handle to the long-lived worker subprocess. The single
/// instance is wrapped in a Mutex inside [`ArbOracle`] so the
/// `&self`-receiver of [`OracleBackend::enclose`] can mutate the
/// stdin and stdout pipes.
struct ArbWorker {
    child: Child,
    stdin: BufWriter<ChildStdin>,
    stdout: BufReader<ChildStdout>,
}

impl ArbWorker {
    /// Spawn a fresh worker process and ping-pong to confirm it is
    /// alive.
    fn spawn(venv_python: &PathBuf, script: &PathBuf) -> Result<Self, ArbError> {
        let mut child = Command::new(venv_python)
            .arg(script)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            // stderr inherited so any traceback from the worker
            // surfaces in the test log.
            .spawn()
            .map_err(ArbError::Spawn)?;
        let stdin = BufWriter::new(child.stdin.take().expect("stdin pipe requested above"));
        let stdout = BufReader::new(child.stdout.take().expect("stdout pipe requested above"));
        let mut worker = Self {
            child,
            stdin,
            stdout,
        };
        // Sanity check: the worker should answer "OK ready" to the
        // initial `ready?` ping. Any other answer means the worker
        // is broken (wrong Python, missing module, ...).
        let response = worker
            .request_raw("ready?")
            .map_err(|e| ArbError::HandshakeFailed(format!("ping: {e}")))?;
        if response != "OK ready" {
            return Err(ArbError::HandshakeFailed(format!(
                "expected `OK ready`, got `{response}`"
            )));
        }
        Ok(worker)
    }

    /// Write a request line, flush, read one response line. Does
    /// not retry on failure (the caller in [`ArbOracle`] handles
    /// the one-shot restart policy).
    fn request_raw(&mut self, line: &str) -> std::io::Result<String> {
        writeln!(self.stdin, "{line}")?;
        self.stdin.flush()?;
        let mut response = String::new();
        let n = self.stdout.read_line(&mut response)?;
        if n == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "worker closed stdout",
            ));
        }
        Ok(response.trim_end().to_string())
    }
}

impl Drop for ArbWorker {
    fn drop(&mut self) {
        // Close stdin to signal EOF; the worker exits cleanly. If
        // the child has already exited (broken pipe path), `kill`
        // is a no-op and `wait` reaps the zombie either way.
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

/// The Arb oracle backend. Owns one worker subprocess; restarts it
/// once on a failed read before propagating the failure.
pub struct ArbOracle {
    worker: Mutex<ArbWorker>,
    venv_python: PathBuf,
    script: PathBuf,
}

impl ArbOracle {
    /// Construct a new `ArbOracle`, spawning the worker subprocess.
    /// Returns `Err(ArbError::VenvNotFound)` if the configured venv
    /// is missing.
    pub fn new() -> Result<Self, ArbError> {
        let venv = default_venv_path();
        let venv_python = venv.join("bin/python3");
        if !venv_python.exists() {
            return Err(ArbError::VenvNotFound(venv));
        }
        let script = worker_script_path();
        if !script.exists() {
            return Err(ArbError::WorkerScriptNotFound(script));
        }
        let worker = ArbWorker::spawn(&venv_python, &script)?;
        Ok(Self {
            worker: Mutex::new(worker),
            venv_python,
            script,
        })
    }

    /// Issue a request line, retrying once with a fresh worker if
    /// the first attempt fails.
    fn request(&self, line: &str) -> Result<String, String> {
        let mut worker = self.worker.lock().expect("ArbOracle mutex poisoned");
        match worker.request_raw(line) {
            Ok(r) => Ok(r),
            Err(first_err) => {
                // Try once to restart and re-issue.
                match ArbWorker::spawn(&self.venv_python, &self.script) {
                    Ok(new_worker) => {
                        *worker = new_worker;
                        worker
                            .request_raw(line)
                            .map_err(|e| format!("restart succeeded but request failed: {e}"))
                    }
                    Err(restart_err) => Err(format!(
                        "request failed ({first_err}); restart also failed ({restart_err})"
                    )),
                }
            }
        }
    }
}

impl OracleBackend for ArbOracle {
    fn enclose(&self, f: FnId, input: u32, mode: RoundingMode, working_prec: u32) -> Enclosure {
        let (fn_id, order) = fnid_to_worker_args(f);
        let mode_str = mode_to_str(mode);
        let request = format!("{fn_id} {order} {input:08x} {mode_str}");
        let response = self
            .request(&request)
            .unwrap_or_else(|e| panic!("Arb oracle request `{request}` failed: {e}"));
        parse_response(&response, working_prec).unwrap_or_else(|e| {
            panic!(
                "Arb oracle response parse failed: request=`{request}` \
                 response=`{response}` error={e}"
            )
        })
    }

    fn name(&self) -> &'static str {
        "Arb"
    }

    /// `true`: the worker's internal Ziv loop produces a single
    /// certified `f32` per call, and further calls at higher
    /// `working_prec` would not change the answer. The verifier
    /// short-circuits its outer Ziv-at-oracle loop and accepts the
    /// worker's first response. ADR-0035.
    fn is_authoritative(&self) -> bool {
        true
    }
}

/// Wire form of a [`RoundingMode`] for the worker protocol.
fn mode_to_str(mode: RoundingMode) -> &'static str {
    match mode {
        RoundingMode::NearestEven => "NE",
        RoundingMode::NearestAway => "RNA",
        RoundingMode::TowardZero => "RZ",
        RoundingMode::TowardPositive => "RP",
        RoundingMode::TowardNegative => "RM",
    }
}

/// Map a `FnId` to the worker's `(fn_id, order_or_dash)` tuple.
/// Panics on a non-Arb-primary `FnId` because the dispatcher
/// (`MetaOracle` in slice p1.5.3) is responsible for routing those
/// to `MpfrOracle` instead.
fn fnid_to_worker_args(f: FnId) -> (&'static str, String) {
    match f {
        FnId::Si => ("si", "-".to_string()),
        FnId::Ci => ("ci", "-".to_string()),
        FnId::Li => ("li", "-".to_string()),
        FnId::Bi => ("bi", "-".to_string()),
        FnId::AiPrime => ("ai_prime", "-".to_string()),
        FnId::BiPrime => ("bi_prime", "-".to_string()),
        FnId::BesselI0 => ("i", "0".to_string()),
        FnId::BesselI1 => ("i", "1".to_string()),
        FnId::BesselIn(n) => ("i", n.to_string()),
        FnId::BesselK0 => ("k", "0".to_string()),
        FnId::BesselK1 => ("k", "1".to_string()),
        FnId::BesselKn(n) => ("k", n.to_string()),
        _ => panic!("ArbOracle::enclose called with non-Arb-primary FnId: {f:?}"),
    }
}

/// Parse one worker response line into an [`Enclosure`].
///
/// Under the ADR-0035 protocol the response is one of:
///
/// - `OK <f32_bits_hex>`: the certified `f32` bit pattern. The
///   enclosure returned is a single-point bracket at that `f32`
///   value (both endpoints equal). The verifier's
///   `certified_round_f32` on the single point under any rounding
///   mode returns the same `f32`.
/// - `INC`: the worker's Ziv loop could not certify a unique `f32`
///   at its max precision. We return an Enclosure with NaN
///   endpoints; the verifier's NaN-aware certified-rounding then
///   reports `OracleInconclusive` because the pfloat kernel's
///   non-NaN output cannot match a NaN certified answer.
/// - `ERR <message>`: an error processing the request. Propagated
///   as a `Result::Err` for the caller to panic.
fn parse_response(line: &str, working_prec: u32) -> Result<Enclosure, String> {
    let mut parts = line.split_whitespace();
    let tag = parts.next().ok_or_else(|| "empty response".to_string())?;
    if tag == "ERR" {
        let msg: String = parts.collect::<Vec<_>>().join(" ");
        return Err(format!("worker reported: {msg}"));
    }
    if tag == "INC" {
        let nan = Float::with_val(working_prec, Special::Nan);
        return Ok(Enclosure {
            lo: nan.clone(),
            hi: nan,
        });
    }
    if tag != "OK" {
        return Err(format!("unexpected response tag: {tag}"));
    }
    let bits_str = parts.next().ok_or_else(|| "missing f32 bits".to_string())?;
    if parts.next().is_some() {
        return Err(format!(
            "trailing data after f32 bits in response: `{line}`"
        ));
    }
    let bits = u32::from_str_radix(bits_str, 16)
        .map_err(|e| format!("parse `{bits_str}` as u32 hex: {e}"))?;
    let value = f32_to_float_endpoint(bits, working_prec);
    Ok(Enclosure {
        lo: value.clone(),
        hi: value,
    })
}

/// Construct a rug `Float` at the requested working precision
/// holding the exact value of an f32 bit pattern. Handles `+0` /
/// `-0` / `+inf` / `-inf` / NaN as the Float special variants;
/// finite values lift through `f32::from_bits` and `Float::with_val`
/// which preserves the exact value at any `working_prec >= 24`.
fn f32_to_float_endpoint(bits: u32, working_prec: u32) -> Float {
    let exp_field = (bits >> 23) & 0xFF;
    let mant = bits & 0x007F_FFFF;
    if exp_field == 0xFF {
        if mant == 0 {
            // Infinity.
            if (bits >> 31) & 1 == 1 {
                return Float::with_val(working_prec, Special::NegInfinity);
            }
            return Float::with_val(working_prec, Special::Infinity);
        }
        // NaN.
        return Float::with_val(working_prec, Special::Nan);
    }
    // Finite: f32 -> f32 round-trip via from_bits is exact, and
    // Float::with_val at any precision >= 24 holds the value
    // without rounding (f32 has at most 24 bits of precision).
    Float::with_val(working_prec, f32::from_bits(bits))
}
