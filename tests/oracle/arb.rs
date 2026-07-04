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
//! [`Enclosed::Bracket`] at the certified `f32`, or as
//! [`Enclosed::Inconclusive`] for `INC`, which the verifier reports
//! as `OracleInconclusive`. A certified NaN (`OK 7fc00000`, a
//! genuinely-undefined true value) is a [`Enclosed::Bracket`] with
//! NaN endpoints and stays distinct from `INC`: conflating the two
//! (both were once NaN-endpoint enclosures) let an inconclusive
//! Arb verdict silently pass as agreement whenever the kernel also
//! returned NaN (pf-41ou). [`ArbOracle::is_authoritative`] returns
//! `true` so the verifier short-circuits its outer Ziv loop and
//! accepts the worker's single answer.

#![cfg(all(unix, feature = "differential-arb"))]

use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};
use std::sync::Mutex;

use rug::float::Special;
use rug::{Complete, Float};

use pfloat::RoundingMode;

use super::types::{Enclosed, Enclosure, FnId, OracleBackend};

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

impl ArbOracle {
    /// Compute the rigorous-enclosure midpoint of the function at
    /// `input` (a binary32 bit pattern), at oracle precision
    /// `oracle_prec >= working_prec + 64`. Returns a [`rug::Float`]
    /// at `oracle_prec` precision.
    ///
    /// The returned midpoint is the centre of the Arb ball at
    /// `oracle_prec`; the ball's radius bounds how far the midpoint
    /// can be from the true value, and at the recommended
    /// `oracle_prec >= working_prec + 64` the radius is well within
    /// the pf-tqzz cross-check tolerance
    /// `2^(error_guard - working_prec) * |midpoint|`.
    ///
    /// The mode parameter is omitted from the request because the
    /// midpoint is mode-independent: every IEEE rounding mode of
    /// the same input value produces the same midpoint candidate
    /// (the rounding happens downstream when comparing against the
    /// kernel's eval(w) intermediate).
    ///
    /// pf-tqzz (slice p1g.3, ADR-0039). Panics on a non-Arb-primary
    /// `FnId`; the caller (cross-check harness) is responsible for
    /// routing other `FnId`s through MPFR (see [`super::cross_check`]).
    pub fn midpoint(&self, f: FnId, input: u32, oracle_prec: u32) -> Result<Float, MidpointError> {
        let (fn_id, order) = fnid_to_worker_args(f);
        let request = format!("MIDPOINT {fn_id} {order} {input:08x} {oracle_prec}");
        let response = self.request(&request).map_err(MidpointError::Worker)?;
        parse_midpoint_response(&response, oracle_prec)
    }
}

/// Errors returned by [`ArbOracle::midpoint`].
#[derive(Debug)]
pub enum MidpointError {
    /// The worker reported `INC` (NaN or unbounded ball); the
    /// midpoint has no finite representation.
    Inconclusive,
    /// The worker reported `ERR <msg>`; the message is preserved
    /// verbatim.
    WorkerError(String),
    /// The Arb subprocess request itself failed (broken pipe,
    /// failed restart, etc.).
    Worker(String),
    /// The response did not parse as a MIDPOINT triple.
    Malformed(String),
}

impl std::fmt::Display for MidpointError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Inconclusive => write!(f, "midpoint inconclusive (ball NaN or unbounded)"),
            Self::WorkerError(msg) => write!(f, "worker error: {msg}"),
            Self::Worker(msg) => write!(f, "worker request failed: {msg}"),
            Self::Malformed(msg) => write!(f, "malformed midpoint response: {msg}"),
        }
    }
}

impl std::error::Error for MidpointError {}

/// Parse a worker MIDPOINT response into a `rug::Float` at
/// `oracle_prec`.
///
/// Response shape per `scripts/arb_oracle_worker.py::handle_midpoint`:
///
/// - `OK <sign> <mantissa_hex> <exponent>` — `sign ∈ {+, -}`;
///   `mantissa_hex` is the absolute integer mantissa as lowercase
///   hex (no `0x` prefix); `exponent` is signed decimal. The value
///   is `sign * mantissa * 2^exponent`. The zero case emits
///   `OK + 0 0` (mantissa `0`, exponent `0`).
/// - `INC` — the ball is NaN or unbounded.
/// - `ERR <message>` — worker-side error.
fn parse_midpoint_response(line: &str, oracle_prec: u32) -> Result<Float, MidpointError> {
    let mut parts = line.split_whitespace();
    let tag = parts
        .next()
        .ok_or_else(|| MidpointError::Malformed("empty response".to_string()))?;
    if tag == "INC" {
        return Err(MidpointError::Inconclusive);
    }
    if tag == "ERR" {
        let msg: String = parts.collect::<Vec<_>>().join(" ");
        return Err(MidpointError::WorkerError(msg));
    }
    if tag != "OK" {
        return Err(MidpointError::Malformed(format!(
            "unexpected response tag: {tag}"
        )));
    }
    let sign_str = parts
        .next()
        .ok_or_else(|| MidpointError::Malformed("missing sign".to_string()))?;
    let mant_hex = parts
        .next()
        .ok_or_else(|| MidpointError::Malformed("missing mantissa".to_string()))?;
    let exp_str = parts
        .next()
        .ok_or_else(|| MidpointError::Malformed("missing exponent".to_string()))?;
    if parts.next().is_some() {
        return Err(MidpointError::Malformed(format!(
            "trailing data in midpoint response: `{line}`"
        )));
    }
    let exp: i64 = exp_str
        .parse()
        .map_err(|e| MidpointError::Malformed(format!("exponent `{exp_str}`: {e}")))?;

    // Zero is encoded as mantissa `0`, exponent `0`; preserve the
    // unsigned zero so the magnitude comparison downstream stays
    // total.
    if mant_hex == "0" {
        return Ok(Float::with_val(oracle_prec, 0));
    }

    // Parse the absolute mantissa as a big integer, apply sign,
    // lift to Float at oracle_prec, then scale by 2^exponent via
    // a binary shift (rug::Float's Shl/Shr scale by powers of two
    // exactly).
    let abs_int = rug::Integer::parse_radix(mant_hex, 16)
        .map_err(|e| MidpointError::Malformed(format!("mantissa hex `{mant_hex}`: {e}")))?
        .complete();
    let signed = if sign_str == "-" {
        -abs_int
    } else if sign_str == "+" {
        abs_int
    } else {
        return Err(MidpointError::Malformed(format!(
            "sign `{sign_str}` not in {{+, -}}"
        )));
    };
    let base = Float::with_val(oracle_prec, &signed);
    let scaled = if exp >= 0 {
        let shift =
            u32::try_from(exp).map_err(|e| MidpointError::Malformed(format!("exp shl: {e}")))?;
        base << shift
    } else {
        let shift =
            u32::try_from(-exp).map_err(|e| MidpointError::Malformed(format!("exp shr: {e}")))?;
        base >> shift
    };
    Ok(scaled)
}

impl OracleBackend for ArbOracle {
    fn enclose(&self, f: FnId, input: u32, mode: RoundingMode, working_prec: u32) -> Enclosed {
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
///
/// Exposed under `pub(super)` so the mpmath oracle (which speaks
/// the same worker protocol modulo the script path) can reuse the
/// FnId-to-string mapping.
pub(super) fn fnid_to_worker_args(f: FnId) -> (&'static str, String) {
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
        // Reciprocal trig (pfloat 1.1, ADR-0056): Arb-primary, native.
        FnId::Cot => ("cot", "-".to_string()),
        FnId::Sec => ("sec", "-".to_string()),
        FnId::Csc => ("csc", "-".to_string()),
        _ => panic!("ArbOracle::enclose called with non-Arb-primary FnId: {f:?}"),
    }
}

/// Parse one worker response line into an [`Enclosed`].
///
/// This is the pure classifier at the heart of pf-41ou: the worker's
/// four response classes map to distinct Rust outcomes, and `INC` is
/// never conflated with a certified NaN. Under the ADR-0035 protocol
/// the response is one of:
///
/// - `OK <f32_bits_hex>`: the certified `f32` bit pattern. Returned
///   as an [`Enclosed::Bracket`] single-point bracket at that `f32`
///   value (both endpoints equal). The verifier's
///   `certified_round_f32` on the single point under any rounding
///   mode returns the same `f32`. When the certified value is NaN
///   (`OK 7fc00000`, a genuinely-undefined true value) both endpoints
///   are NaN and the bracket certifies `Some(f32::NAN)`.
/// - `INC`: the worker's Ziv loop could not certify a unique `f32`
///   at its max precision. Returned as [`Enclosed::Inconclusive`], a
///   value the verifier maps to [`super::types::Verdict::OracleInconclusive`].
///   It is emphatically NOT a certified NaN: an inconclusive Arb
///   verdict must never count as agreement (the earlier code returned
///   a NaN-endpoint enclosure here, which the verifier read as a
///   certified NaN and silently passed as `Ok` whenever the kernel
///   also returned NaN — pf-41ou).
/// - `ERR <message>`: an error processing the request. Propagated
///   as a `Result::Err` for the caller to panic.
///
/// `working_prec` is threaded through only for the `OK` finite-value
/// path; the `INC` and `ERR` classifications are precision-free.
fn parse_response(line: &str, working_prec: u32) -> Result<Enclosed, String> {
    let mut parts = line.split_whitespace();
    let tag = parts.next().ok_or_else(|| "empty response".to_string())?;
    if tag == "ERR" {
        let msg: String = parts.collect::<Vec<_>>().join(" ");
        return Err(format!("worker reported: {msg}"));
    }
    if tag == "INC" {
        // The worker abstained. This is NOT a certified NaN: return
        // the dedicated inconclusive outcome so the verifier reports
        // OracleInconclusive rather than folding it into an Ok/NaN
        // agreement (pf-41ou).
        return Ok(Enclosed::Inconclusive);
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
    // A certified NaN (`OK 7fc00000`) becomes a bracket with NaN
    // endpoints here; that is a genuinely-undefined true value the
    // verifier accepts against a NaN kernel output, and it stays
    // distinct from the `INC` outcome above.
    let value = f32_to_float_endpoint(bits, working_prec);
    Ok(Enclosed::Bracket(Enclosure {
        lo: value.clone(),
        hi: value,
    }))
}

/// Wrapper that exposes [`parse_response`] across modules; the
/// mpmath oracle uses identical response parsing because both
/// workers speak the same wire protocol.
pub(super) fn parse_response_external(line: &str, working_prec: u32) -> Result<Enclosed, String> {
    parse_response(line, working_prec)
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

#[cfg(test)]
mod tests {
    use super::*;

    // pf-41ou. The worker's four response classes must map to four
    // distinct Rust outcomes; the honesty bug was `INC` collapsing
    // into the same NaN-endpoint enclosure a *certified* NaN uses,
    // so an inconclusive Arb verdict silently passed as agreement.
    // These exercise the pure line -> Enclosed classifier directly,
    // so they run without a live Arb worker (no python-flint needed).

    #[test]
    fn inc_classifies_as_inconclusive_not_a_bracket() {
        match parse_response("INC", 64) {
            Ok(Enclosed::Inconclusive) => {}
            other => panic!("`INC` must map to Enclosed::Inconclusive, got {other:?}"),
        }
    }

    #[test]
    fn certified_nan_is_a_bracket_distinct_from_inc() {
        // `OK 7fc00000` is a *certified* NaN true value: a bracket
        // with NaN endpoints, emphatically NOT the inconclusive
        // outcome. This is the distinction pf-41ou preserves.
        match parse_response("OK 7fc00000", 64) {
            Ok(Enclosed::Bracket(enc)) => {
                assert!(
                    enc.lo.is_nan() && enc.hi.is_nan(),
                    "certified NaN must have NaN endpoints"
                );
            }
            other => panic!("certified NaN must be Enclosed::Bracket, got {other:?}"),
        }
    }

    #[test]
    fn certified_finite_is_a_single_point_bracket() {
        // `OK 3f800000` = 1.0_f32.
        match parse_response("OK 3f800000", 64) {
            Ok(Enclosed::Bracket(enc)) => {
                assert_eq!(enc.lo.to_f32_round(rug::float::Round::Nearest), 1.0);
                assert_eq!(enc.hi.to_f32_round(rug::float::Round::Nearest), 1.0);
            }
            other => panic!("certified finite must be Enclosed::Bracket, got {other:?}"),
        }
    }

    #[test]
    fn err_propagates_as_error() {
        assert!(
            parse_response("ERR something broke", 64).is_err(),
            "`ERR` must propagate as Result::Err"
        );
    }
}
