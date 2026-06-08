//! Subprocess driver and exact dyadic codec for the complex-Arb (`acb`) worker
//! (`scripts/acb_complex_worker.py`), the independent componentwise
//! certified-rounding backstop (ADR-0092, C5).
//!
//! Reaches Arb purely out of process through the python-flint worker; nothing
//! here links FLINT/Arb. The codec is exact in both directions: each operand
//! component is sent as `sign * mantissa * 2^exp` and the worker's reply (the
//! rigorous per-component enclosure, also dyadic) is lifted back to a `BigFloat`
//! losslessly, so no decimal crosses the boundary. This mirrors
//! `pfloat-ball/tests/common/arb_bracket.rs`, extended to the two components of
//! a complex result.

use core::cmp::Ordering;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};

use pfloat::{BigFloat, Parts, RoundingMode, Sign};

const NE: RoundingMode = RoundingMode::NearestEven;

/// One component (real or imaginary) of an `acb` result.
#[derive(Debug, Clone)]
pub enum Comp {
    /// `lo <= value <= hi`, both exact dyadics.
    Finite {
        lo: BigFloat,
        hi: BigFloat,
    },
    /// The component is NaN.
    Nan,
    /// The component is entirely +inf / -inf.
    PosInf,
    NegInf,
    /// Sign indeterminate at this precision.
    Inconclusive,
}

/// A complex result's two component brackets.
#[derive(Debug, Clone)]
pub struct ComplexBracket {
    pub re: Comp,
    pub im: Comp,
}

/// A live handle to the python-flint `acb` worker, speaking `CBRACKET`.
pub struct AcbComplexWorker {
    child: Child,
    stdin: BufWriter<ChildStdin>,
    stdout: BufReader<ChildStdout>,
}

impl AcbComplexWorker {
    /// Spawn the worker, resolving the venv (`PFLOAT_ARB_ORACLE_VENV` or
    /// `~/.cache/pfloat-arb-oracle/venv`, the same venv the ball lane uses) and
    /// the in-tree worker script. Panics with a clear message if the venv or
    /// script is missing; the lane is env-gated so a developer without the venv
    /// simply does not run it.
    pub fn spawn() -> Self {
        let venv = std::env::var("PFLOAT_ARB_ORACLE_VENV")
            .map(PathBuf::from)
            .unwrap_or_else(|_| {
                PathBuf::from(std::env::var("HOME").expect("HOME set for the Arb venv path"))
                    .join(".cache/pfloat-arb-oracle/venv")
            });
        let python = venv.join("bin/python3");
        assert!(
            python.exists(),
            "Arb venv python not found at {python:?}; run scripts/setup_arb_oracle.sh"
        );
        let script =
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../scripts/acb_complex_worker.py");
        assert!(script.exists(), "worker script not found at {script:?}");

        let mut child = Command::new(&python)
            .arg(&script)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .expect("spawn acb worker");
        let stdin = BufWriter::new(child.stdin.take().expect("stdin pipe"));
        let stdout = BufReader::new(child.stdout.take().expect("stdout pipe"));
        let mut w = Self {
            child,
            stdin,
            stdout,
        };
        let resp = w.request("ready?");
        assert_eq!(resp, "OK ready", "worker handshake failed: got `{resp}`");
        w
    }

    fn request(&mut self, line: &str) -> String {
        writeln!(self.stdin, "{line}").expect("write worker stdin");
        self.stdin.flush().expect("flush worker stdin");
        let mut resp = String::new();
        let n = self
            .stdout
            .read_line(&mut resp)
            .expect("read worker stdout");
        assert!(n > 0, "worker closed stdout");
        resp.trim_end().to_string()
    }

    /// The rigorous componentwise enclosure of `fn_id(z [, w])` at `oracle_prec`
    /// Arb bits. `w` is `Some` for the binary ops (cadd/csub/cmul/cdiv).
    pub fn cbracket(
        &mut self,
        fn_id: &str,
        oracle_prec: u32,
        z: (&BigFloat, &BigFloat),
        w: Option<(&BigFloat, &BigFloat)>,
    ) -> ComplexBracket {
        let mut line = format!("CBRACKET {fn_id} {oracle_prec} {} {}", tri(z.0), tri(z.1));
        if let Some((wr, wi)) = w {
            line.push_str(&format!(" {} {}", tri(wr), tri(wi)));
        }
        let resp = self.request(&line);
        parse_cbracket(fn_id, &resp)
    }
}

impl Drop for AcbComplexWorker {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

/// `<sign> <mantissa_hex> <exp>` of an exact finite `BigFloat`.
fn tri(bf: &BigFloat) -> String {
    let (s, m, e) = bigfloat_to_dyadic(bf).expect("finite operand");
    format!("{s} {m} {e}")
}

/// Parse a `CBRACKET` reply (`OK <re_component> <im_component>`) into a
/// [`ComplexBracket`]; panic on `ERR` / malformed.
fn parse_cbracket(fn_id: &str, resp: &str) -> ComplexBracket {
    let toks: Vec<&str> = resp.split_whitespace().collect();
    assert!(
        !toks.is_empty() && toks[0] == "OK",
        "cbracket {fn_id}: unexpected response `{resp}`"
    );
    let mut idx = 1;
    let re = parse_component(&toks, &mut idx, fn_id, resp);
    let im = parse_component(&toks, &mut idx, fn_id, resp);
    assert_eq!(
        idx,
        toks.len(),
        "cbracket {fn_id}: trailing tokens in `{resp}`"
    );
    ComplexBracket { re, im }
}

/// Parse one component starting at `*idx`, advancing it past the tokens
/// consumed. `F` consumes 7 tokens (`F` + two dyadic triples); the non-finite
/// markers consume 1.
fn parse_component(toks: &[&str], idx: &mut usize, fn_id: &str, resp: &str) -> Comp {
    let tag = toks.get(*idx).unwrap_or_else(|| {
        panic!("cbracket {fn_id}: missing component at {idx} in `{resp}`");
    });
    match *tag {
        "F" => {
            assert!(
                *idx + 6 < toks.len(),
                "cbracket {fn_id}: short F in `{resp}`"
            );
            let lo = dyadic_to_bigfloat(toks[*idx + 1], toks[*idx + 2], toks[*idx + 3]);
            let hi = dyadic_to_bigfloat(toks[*idx + 4], toks[*idx + 5], toks[*idx + 6]);
            *idx += 7;
            Comp::Finite { lo, hi }
        }
        "N" => {
            *idx += 1;
            Comp::Nan
        }
        "P" => {
            *idx += 1;
            Comp::PosInf
        }
        "M" => {
            *idx += 1;
            Comp::NegInf
        }
        "Q" => {
            *idx += 1;
            Comp::Inconclusive
        }
        other => panic!("cbracket {fn_id}: unknown component tag `{other}` in `{resp}`"),
    }
}

/// Whether the worker venv is present (the lane skips when it is not).
pub fn venv_available() -> bool {
    let venv = std::env::var("PFLOAT_ARB_ORACLE_VENV")
        .map(PathBuf::from)
        .unwrap_or_else(|_| {
            PathBuf::from(std::env::var("HOME").unwrap_or_default())
                .join(".cache/pfloat-arb-oracle/venv")
        });
    venv.join("bin/python3").exists()
}

/// Gate an acb-driven test: `true` to proceed, `false` to skip (venv absent).
/// When `PFLOAT_ARB_REQUIRED` is set -- the per-release gate -- a missing venv
/// is a HARD failure instead of a silent skip, so the independent backstop
/// cannot quietly no-op and still report green.
pub fn acb_lane_available(test_name: &str) -> bool {
    if venv_available() {
        return true;
    }
    assert!(
        std::env::var("PFLOAT_ARB_REQUIRED").is_err(),
        "PFLOAT_ARB_REQUIRED is set but the Arb venv is absent: {test_name} cannot run the \
         independent containment check (run scripts/setup_arb_oracle.sh)"
    );
    eprintln!("skip: Arb venv absent ({test_name}); run scripts/setup_arb_oracle.sh");
    false
}

/// `(sign_str, abs_mantissa_hex, exp)` of an exact finite `BigFloat` value
/// `sign * mantissa * 2^exp`; `None` for non-finite. The integer mantissa is
/// the left-aligned limbs read big-endian (the same shift the ball codec uses).
pub fn bigfloat_to_dyadic(bf: &BigFloat) -> Option<(&'static str, String, i64)> {
    match bf.parts() {
        Parts::Zero { sign } => Some((sign_str(sign), "0".to_string(), 0)),
        Parts::Normal {
            sign,
            exponent,
            mantissa,
            ..
        } => {
            let stored_bits = mantissa.len() as i64 * 64;
            let exp = exponent + 1 - stored_bits;
            let hex: String = mantissa.iter().rev().map(|l| format!("{l:016x}")).collect();
            Some((sign_str(sign), hex, exp))
        }
        _ => None,
    }
}

fn sign_str(s: Sign) -> &'static str {
    if matches!(s, Sign::Negative) {
        "-"
    } else {
        "+"
    }
}

/// Exact `BigFloat` for `sign * int(man_hex) * 2^exp`, built at a working
/// precision wide enough to hold the mantissa losslessly: Horner over the hex
/// digits (each `* 16 + d` exact while the accumulator stays under `2^work`),
/// then an exact power-of-two scale.
pub fn dyadic_to_bigfloat(sign: &str, man_hex: &str, exp: &str) -> BigFloat {
    let work = (man_hex.len() as u32) * 4 + 8;
    let sixteen = BigFloat::try_from_i64_exact(16, work).unwrap();
    let mut acc = BigFloat::try_from_i64_exact(0, work).unwrap();
    for ch in man_hex.chars() {
        let d = i64::from(ch.to_digit(16).expect("hex digit"));
        acc = acc
            .mul(&sixteen, NE)
            .0
            .add(&BigFloat::try_from_i64_exact(d, work).unwrap(), NE)
            .0;
    }
    if sign == "-" {
        acc = acc.negated();
    }
    let e: i64 = exp.parse().expect("exp integer");
    acc.scale_by_pow2(e).0
}

/// Test helper: round-trip a `BigFloat` through the dyadic wire form. The codec
/// is lossless, so the result must equal the input exactly.
pub fn encode_decode(bf: &BigFloat) -> BigFloat {
    let (s, m, e) = bigfloat_to_dyadic(bf).expect("finite");
    dyadic_to_bigfloat(s, &m, &e.to_string())
}

/// Whether `a == b` in value AND sign (so `+0` and `-0` differ); both NaN also
/// counts as equal. The componentwise certified-rounding comparison.
pub fn exact_signed_eq(a: &BigFloat, b: &BigFloat) -> bool {
    if a.is_nan() || b.is_nan() {
        return a.is_nan() && b.is_nan();
    }
    matches!(a.partial_cmp(b).0, Some(Ordering::Equal))
        && a.is_sign_negative() == b.is_sign_negative()
}
