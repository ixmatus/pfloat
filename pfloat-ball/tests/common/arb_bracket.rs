//! Subprocess driver and exact dyadic codec for the Arb worker's
//! `BRACKET` verb (pf-fe5f, ADR-0078 follow-up).
//!
//! Reaches Arb purely out of process through pfloat's python-flint
//! worker; nothing here links FLINT/Arb. The codec is exact in both
//! directions: a `BigFloat` is sent as `sign * mantissa * 2^exp` and the
//! worker's reply (the rigorous enclosure, also dyadic) is lifted back to
//! a `BigFloat` losslessly, so no decimal crosses the boundary.

use core::cmp::Ordering;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};

use pfloat::{BigFloat, Parts, RoundingMode, Sign};
use pfloat_ball::Ball;

const NE: RoundingMode = RoundingMode::NearestEven;

/// One `BRACKET` response: the rigorous enclosure, or a non-finite verdict.
#[derive(Debug)]
pub enum Bracket {
    /// `lo <= f(x) <= hi`, both exact.
    Finite {
        lo: BigFloat,
        hi: BigFloat,
    },
    /// Arb's ball is NaN (a pole / out-of-domain point).
    Nan,
    /// The result lies entirely at +inf / -inf.
    PosInf,
    NegInf,
    /// Sign indeterminate at this precision.
    Inconclusive,
}

/// A live handle to the python-flint Arb worker, speaking `BRACKET`.
pub struct ArbBracketWorker {
    child: Child,
    stdin: BufWriter<ChildStdin>,
    stdout: BufReader<ChildStdout>,
}

impl ArbBracketWorker {
    /// Spawn the worker, resolving the venv (`PFLOAT_ARB_ORACLE_VENV` or
    /// `~/.cache/pfloat-arb-oracle/venv`) and the in-tree worker script
    /// (`../scripts/arb_oracle_worker.py`, relative to this crate). Panics
    /// with a clear message if the venv or script is missing; the lane is
    /// env-gated so a developer without the venv simply does not run it.
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
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../scripts/arb_oracle_worker.py");
        assert!(script.exists(), "worker script not found at {script:?}");

        let mut child = Command::new(&python)
            .arg(&script)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .expect("spawn Arb worker");
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

    /// The rigorous enclosure of `fn_id` over input(s) at `oracle_prec`
    /// Arb bits. `y` is `Some` for binary functions (add/sub/mul/div/
    /// atan2/hypot).
    pub fn bracket(
        &mut self,
        fn_id: &str,
        oracle_prec: u32,
        x: &BigFloat,
        y: Option<&BigFloat>,
    ) -> Bracket {
        let (xs, xm, xe) = bigfloat_to_dyadic(x).expect("finite x operand");
        let mut line = format!("BRACKET {fn_id} {oracle_prec} {xs} {xm} {xe}");
        if let Some(y) = y {
            let (ys, ym, ye) = bigfloat_to_dyadic(y).expect("finite y operand");
            line.push_str(&format!(" {ys} {ym} {ye}"));
        }
        let resp = self.request(&line);
        let parts: Vec<&str> = resp.split_whitespace().collect();
        match parts.as_slice() {
            ["OK", ls, lm, le, hs, hm, he] => Bracket::Finite {
                lo: dyadic_to_bigfloat(ls, lm, le),
                hi: dyadic_to_bigfloat(hs, hm, he),
            },
            ["NAN"] => Bracket::Nan,
            ["POS_INF"] => Bracket::PosInf,
            ["NEG_INF"] => Bracket::NegInf,
            ["INC"] => Bracket::Inconclusive,
            _ => panic!("BRACKET {fn_id}: unexpected response `{resp}`"),
        }
    }
}

impl Drop for ArbBracketWorker {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

/// Whether the Arb worker venv is present. The lane skips when it is not,
/// so a `--features differential-arb` build on a machine without the venv
/// does not fail; this is a per-release lane needing
/// `scripts/setup_arb_oracle.sh`.
pub fn venv_available() -> bool {
    let venv = std::env::var("PFLOAT_ARB_ORACLE_VENV")
        .map(PathBuf::from)
        .unwrap_or_else(|_| {
            PathBuf::from(std::env::var("HOME").unwrap_or_default())
                .join(".cache/pfloat-arb-oracle/venv")
        });
    venv.join("bin/python3").exists()
}

/// SOUND containment: the ball's denoted interval contains Arb's rigorous
/// enclosure `[lo, hi]` of f(x). Then `f(x) in [lo, hi] subset of ball`,
/// so `f(x) in ball` follows whatever f(x) is inside `[lo, hi]`. The
/// reverse check (`ball subset of [lo, hi]`) would PASS exactly when the
/// ball is too narrow to be sound, so it is never the assertion.
pub fn ball_contains(ball: &Ball<BigFloat>, lo: &BigFloat, hi: &BigFloat) -> bool {
    ball.lower().partial_cmp(lo).0 != Some(Ordering::Greater)
        && ball.upper().partial_cmp(hi).0 != Some(Ordering::Less)
}

/// `(sign_str, abs_mantissa_hex, exp)` of an exact finite `BigFloat`
/// value `sign * mantissa * 2^exp`. `None` for non-finite. Matches the
/// worker's dyadic input format; the construction is exact (the integer
/// mantissa is the left-aligned limbs read big-endian, the same shift
/// `tests/differential/mod.rs::bigfloat_to_rug` uses).
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
/// precision wide enough to hold the mantissa losslessly: Horner over the
/// hex digits (each `* 16 + d` is exact while the accumulator stays under
/// `2^work`), then an exact power-of-two scale.
pub fn dyadic_to_bigfloat(sign: &str, man_hex: &str, exp: &str) -> BigFloat {
    let work = (man_hex.len() as u32) * 4 + 8;
    let sixteen = BigFloat::try_from_i64_exact(16, work).unwrap();
    let mut acc = BigFloat::try_from_i64_exact(0, work).unwrap();
    for c in man_hex.chars() {
        let d = i64::from(c.to_digit(16).expect("hex digit"));
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

/// Test helper: round-trip a `BigFloat` through the dyadic wire form. The
/// codec is lossless, so the result must equal the input exactly.
pub fn encode_decode(bf: &BigFloat) -> BigFloat {
    let (s, m, e) = bigfloat_to_dyadic(bf).expect("finite");
    dyadic_to_bigfloat(s, &m, &e.to_string())
}
