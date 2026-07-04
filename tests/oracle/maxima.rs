//! Maxima oracle backend per ADR-0035 Tier 6: the third independent
//! oracle, used as a sampling layer (hand-derived corpus +
//! tie-breakers + N-sample per release) rather than a full f32
//! sweep.
//!
//! Maxima's Macsyma lineage from 1968 is completely independent of
//! FLINT/Arb and mpmath; three-way agreement between Arb, mpmath,
//! and Maxima on a sampled input is the strongest evidence we can
//! get short of formal proof that the certified `f32` is correct.
//!
//! Function coverage caveats per slice p1.10's probe:
//!
//! - `bessel_i(n, x)` at very small subnormals can trigger
//!   Maxima's "Exceeded maximum allowed fpprec"; the worker
//!   reports those as `INC` rather than a certified `f32`.
//!   Arb + mpmath two-oracle agreement covers that input class.
//! - `li` is composed via `ei(log(x))` (Maxima has no direct
//!   logarithmic integral primitive); the composition's accuracy
//!   inherits Maxima's `expintegral_ei` precision.
//!
//! The Rust side mirrors the [`super::arb::ArbOracle`] +
//! [`super::mpmath::MpmathOracle`] shape; the wire protocol is
//! identical so the response parser is reused.
//!
//! Cost: ~500ms-1s per request (Maxima startup per call); not
//! suitable for full-sweep use. The sampling layer keeps total
//! cost bounded (~50-100 requests per release).

#![cfg(all(unix, feature = "differential-arb"))]

use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};
use std::sync::Mutex;

use pfloat::RoundingMode;

use super::arb::ArbError;
use super::types::{Enclosed, FnId, OracleBackend};

/// Path to the in-tree Maxima worker launcher script. The
/// launcher is a nix-shell wrapper so the test harness can invoke
/// Maxima without requiring a system install; the script's
/// shebang pulls Maxima and Python from nixpkgs.
fn worker_script_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scripts/maxima_oracle_worker.sh")
}

struct MaximaWorker {
    child: Child,
    stdin: BufWriter<ChildStdin>,
    stdout: BufReader<ChildStdout>,
}

impl MaximaWorker {
    fn spawn(script: &PathBuf) -> Result<Self, ArbError> {
        let mut child = Command::new(script)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .map_err(ArbError::Spawn)?;
        let stdin = BufWriter::new(child.stdin.take().expect("stdin pipe requested above"));
        let stdout = BufReader::new(child.stdout.take().expect("stdout pipe requested above"));
        let mut worker = Self {
            child,
            stdin,
            stdout,
        };
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

impl Drop for MaximaWorker {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

/// The Maxima oracle backend. Owns one nix-shell + Maxima
/// subprocess; restarts it once on a failed read.
pub struct MaximaOracle {
    worker: Mutex<MaximaWorker>,
    script: PathBuf,
}

impl MaximaOracle {
    pub fn new() -> Result<Self, ArbError> {
        let script = worker_script_path();
        if !script.exists() {
            return Err(ArbError::WorkerScriptNotFound(script));
        }
        let worker = MaximaWorker::spawn(&script)?;
        Ok(Self {
            worker: Mutex::new(worker),
            script,
        })
    }

    fn request(&self, line: &str) -> Result<String, String> {
        let mut worker = self.worker.lock().expect("MaximaOracle mutex poisoned");
        match worker.request_raw(line) {
            Ok(r) => Ok(r),
            Err(first_err) => match MaximaWorker::spawn(&self.script) {
                Ok(new_worker) => {
                    *worker = new_worker;
                    worker
                        .request_raw(line)
                        .map_err(|e| format!("restart succeeded but request failed: {e}"))
                }
                Err(restart_err) => Err(format!(
                    "request failed ({first_err}); restart also failed ({restart_err})"
                )),
            },
        }
    }
}

impl OracleBackend for MaximaOracle {
    fn enclose(&self, f: FnId, input: u32, mode: RoundingMode, working_prec: u32) -> Enclosed {
        let (fn_id, order) = super::arb::fnid_to_worker_args(f);
        let mode_str = mode_to_str(mode);
        let request = format!("{fn_id} {order} {input:08x} {mode_str}");
        let response = self
            .request(&request)
            .unwrap_or_else(|e| panic!("Maxima oracle request `{request}` failed: {e}"));
        super::arb::parse_response_external(&response, working_prec).unwrap_or_else(|e| {
            panic!(
                "Maxima oracle response parse failed: request=`{request}` \
                 response=`{response}` error={e}"
            )
        })
    }

    fn name(&self) -> &'static str {
        "Maxima"
    }

    fn is_authoritative(&self) -> bool {
        true
    }
}

fn mode_to_str(mode: RoundingMode) -> &'static str {
    match mode {
        RoundingMode::NearestEven => "NE",
        RoundingMode::NearestAway => "RNA",
        RoundingMode::TowardZero => "RZ",
        RoundingMode::TowardPositive => "RP",
        RoundingMode::TowardNegative => "RM",
    }
}
