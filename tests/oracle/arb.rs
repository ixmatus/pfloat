//! Arb oracle backend: the second [`OracleBackend`] implementation,
//! sized for the twelve `FnId`s the MPFR backend cannot cover
//! (`Si`, `Ci`, `Li`, `Bi`, `Ai_prime`, `Bi_prime`,
//! `BesselI{0,1,n}`, `BesselK{0,1,n}`).
//! Evaluations happen out-of-process in a long-lived `python-flint`
//! worker; the worker reads one request per line on its stdin and
//! emits one enclosure per line on its stdout. This module owns the
//! subprocess, the request / response protocol, and the conversion
//! between Arb's decimal enclosure mantissas and rug's [`Float`] for
//! the [`Enclosure`] type.
//!
//! See ADR-0034 for the design (the "Arb backend posture (next
//! slice)" section is exactly this slice). FLINT and Arb are both
//! LGPL; keeping them in a Python subprocess means they never
//! enter the shipped Rust crate's link graph.
//!
//! ## Venv resolution
//!
//! The Arb backend needs a Python venv with `python-flint`
//! installed. The default path is `${HOME}/.cache/pfloat-arb-oracle/venv`;
//! override via the `PFLOAT_ARB_ORACLE_VENV` env var. The
//! [`scripts/setup_arb_oracle.sh`](../../scripts/setup_arb_oracle.sh)
//! helper creates and verifies the venv idempotently.
//!
//! ## Worker protocol
//!
//! Request: `<fn_id> <order_or_dash> <input_bits_hex> <working_prec>`.
//! Response: `OK <lo_decimal> <hi_decimal>` or `ERR <message>`.
//! Endpoints come back as `<mantissa>e<exp>` decimals (parsed by
//! rug's `Float::parse` accepting scientific notation) or as
//! `nan` / `inf` / `-inf` for non-finite cases. The Python worker
//! adds a `+/-1` absorption to the decimal mantissas so the
//! resulting rug `Float` endpoints still rigorously bracket the
//! true value after the binary-from-decimal parse.

#![cfg(all(unix, feature = "differential-arb"))]

use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};
use std::sync::Mutex;

use rug::float::Special;
use rug::Float;

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
    fn enclose(&self, f: FnId, input: u32, working_prec: u32) -> Enclosure {
        let (fn_id, order) = fnid_to_worker_args(f);
        let request = format!("{fn_id} {order} {input:08x} {working_prec}");
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
fn parse_response(line: &str, working_prec: u32) -> Result<Enclosure, String> {
    let mut parts = line.split_whitespace();
    let tag = parts.next().ok_or_else(|| "empty response".to_string())?;
    if tag == "ERR" {
        let msg: String = parts.collect::<Vec<_>>().join(" ");
        return Err(format!("worker reported: {msg}"));
    }
    if tag != "OK" {
        return Err(format!("unexpected response tag: {tag}"));
    }
    let lo_str = parts
        .next()
        .ok_or_else(|| "missing lo endpoint".to_string())?;
    let hi_str = parts
        .next()
        .ok_or_else(|| "missing hi endpoint".to_string())?;
    let lo = parse_float_endpoint(lo_str, working_prec)?;
    let hi = parse_float_endpoint(hi_str, working_prec)?;
    Ok(Enclosure { lo, hi })
}

/// Parse one decimal-or-special endpoint into a rug `Float` at
/// `working_prec` bits.
fn parse_float_endpoint(s: &str, working_prec: u32) -> Result<Float, String> {
    match s {
        "nan" => Ok(Float::with_val(working_prec, Special::Nan)),
        "inf" => Ok(Float::with_val(working_prec, Special::Infinity)),
        "-inf" => Ok(Float::with_val(working_prec, Special::NegInfinity)),
        _ => Float::parse(s)
            .map(|inc| Float::with_val(working_prec, inc))
            .map_err(|e| format!("parse `{s}`: {e}")),
    }
}
