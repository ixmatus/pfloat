//! mpmath oracle backend per ADR-0035 Tier 2: a second independent
//! oracle that runs the same functions as the Arb backend through
//! mpmath (BSD-licensed pure-Python multi-precision library, no
//! shared code lineage with FLINT/Arb).
//!
//! Structurally identical to [`super::arb::ArbOracle`]: subprocess
//! ownership, worker protocol, single-point [`Enclosure`] return.
//! The two oracles must agree per input on the certified `f32`;
//! any divergence is a strong signal that one of the two has a
//! silent defect, since correlated bugs between two independent
//! libraries are vanishingly unlikely.
//!
//! Currently gated under the existing `differential-arb` feature
//! since the venv setup includes mpmath alongside python-flint. A
//! future slice may split `differential-mpmath` out if there is
//! demand for the mpmath backend alone.
//!
//! See `scripts/mpmath_oracle_worker.py` for the worker
//! implementation and protocol notes.

#![cfg(all(unix, feature = "differential-arb"))]

use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};
use std::sync::Mutex;

use pfloat::RoundingMode;

use super::arb::{
    ArbError,
    // We reuse ArbError to avoid duplicating the error enum; the
    // failure modes are identical.
};
use super::types::{Enclosure, FnId, OracleBackend};

/// Resolve the venv path: env var override, otherwise
/// `${HOME}/.cache/pfloat-arb-oracle/venv` (same venv as the Arb
/// worker; mpmath sits alongside python-flint).
fn default_venv_path() -> PathBuf {
    if let Ok(p) = std::env::var("PFLOAT_ARB_ORACLE_VENV") {
        return PathBuf::from(p);
    }
    let home = std::env::var("HOME").expect("HOME must be set for the default oracle venv path");
    PathBuf::from(home).join(".cache/pfloat-arb-oracle/venv")
}

fn worker_script_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("scripts/mpmath_oracle_worker.py")
}

struct MpmathWorker {
    child: Child,
    stdin: BufWriter<ChildStdin>,
    stdout: BufReader<ChildStdout>,
}

impl MpmathWorker {
    fn spawn(venv_python: &PathBuf, script: &PathBuf) -> Result<Self, ArbError> {
        let mut child = Command::new(venv_python)
            .arg(script)
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

impl Drop for MpmathWorker {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

/// The mpmath oracle backend. Owns one worker subprocess.
pub struct MpmathOracle {
    worker: Mutex<MpmathWorker>,
    venv_python: PathBuf,
    script: PathBuf,
}

impl MpmathOracle {
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
        let worker = MpmathWorker::spawn(&venv_python, &script)?;
        Ok(Self {
            worker: Mutex::new(worker),
            venv_python,
            script,
        })
    }

    fn request(&self, line: &str) -> Result<String, String> {
        let mut worker = self.worker.lock().expect("MpmathOracle mutex poisoned");
        match worker.request_raw(line) {
            Ok(r) => Ok(r),
            Err(first_err) => match MpmathWorker::spawn(&self.venv_python, &self.script) {
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

impl OracleBackend for MpmathOracle {
    fn enclose(&self, f: FnId, input: u32, mode: RoundingMode, working_prec: u32) -> Enclosure {
        let (fn_id, order) = super::arb::fnid_to_worker_args(f);
        let mode_str = mode_to_str(mode);
        let request = format!("{fn_id} {order} {input:08x} {mode_str}");
        let response = self
            .request(&request)
            .unwrap_or_else(|e| panic!("mpmath oracle request `{request}` failed: {e}"));
        super::arb::parse_response_external(&response, working_prec).unwrap_or_else(|e| {
            panic!(
                "mpmath oracle response parse failed: request=`{request}` \
                 response=`{response}` error={e}"
            )
        })
    }

    fn name(&self) -> &'static str {
        "mpmath"
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
