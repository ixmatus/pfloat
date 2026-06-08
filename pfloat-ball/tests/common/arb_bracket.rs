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
        parse_bracket(fn_id, &resp)
    }

    /// The rigorous enclosure of `fn_id` over the input INTERVAL(s)
    /// `[mid - rad, mid + rad]`, via the worker's `BRACKETI` verb. `y` is
    /// `Some((y_mid, y_rad))` for binary functions (add/sub/mul/div/atan2/
    /// hypot).
    ///
    /// Where [`bracket`](Self::bracket) evaluates `f` at a single point, this
    /// brackets `f` over the whole input ball, so it can witness a result ball
    /// that fails to enclose an interior extremum of a non-monotonic function
    /// (the range-soundness check five point witnesses are structurally blind
    /// to). The radii are passed as `BigFloat`s (a `Mag` lifts via
    /// `to_bigfloat`); a zero radius makes `BRACKETI` reduce to `bracket`
    /// bit-for-bit.
    pub fn bracket_interval(
        &mut self,
        fn_id: &str,
        oracle_prec: u32,
        x_mid: &BigFloat,
        x_rad: &BigFloat,
        y: Option<(&BigFloat, &BigFloat)>,
    ) -> Bracket {
        let (xms, xmm, xme) = bigfloat_to_dyadic(x_mid).expect("finite x midpoint");
        let (xrs, xrm, xre) = bigfloat_to_dyadic(x_rad).expect("finite x radius");
        let mut line =
            format!("BRACKETI {fn_id} {oracle_prec} {xms} {xmm} {xme} {xrs} {xrm} {xre}");
        if let Some((y_mid, y_rad)) = y {
            let (yms, ymm, yme) = bigfloat_to_dyadic(y_mid).expect("finite y midpoint");
            let (yrs, yrm, yre) = bigfloat_to_dyadic(y_rad).expect("finite y radius");
            line.push_str(&format!(" {yms} {ymm} {yme} {yrs} {yrm} {yre}"));
        }
        let resp = self.request(&line);
        parse_bracket(fn_id, &resp)
    }
}

/// Parse a worker BRACKET / BRACKETI reply into a [`Bracket`]. The point and
/// interval verbs share the reply grammar (the worker encodes both through one
/// helper), so they share one parser.
fn parse_bracket(fn_id: &str, resp: &str) -> Bracket {
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
        _ => panic!("bracket {fn_id}: unexpected response `{resp}`"),
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

/// Gate an Arb-driven test: `true` to proceed, `false` to skip (venv
/// absent). When `PFLOAT_ARB_REQUIRED` is set -- the per-release gate -- a
/// missing venv is a HARD failure instead of a silent skip, so the
/// independent backstop cannot quietly no-op and still report green. The
/// release runs `PFLOAT_ARB_REQUIRED=1 cargo test -p pfloat-ball --features
/// differential-arb` and learns immediately if the worker never ran.
pub fn arb_lane_available(test_name: &str) -> bool {
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

/// The SOUND independent check: the ball's denoted interval must not be
/// provably DISJOINT from Arb's rigorous enclosure `[lo, hi]` of
/// `f(witness)`. Since `f(witness) in [lo, hi]` (Arb is rigorous) and FTIA
/// requires `f(witness) in ball`, a ball whose interval lies entirely
/// outside `[lo, hi]` excludes `f(witness)` and is unsound. Returns false
/// exactly in that provable-exclusion case:
///
/// ```text
/// ball.lower() <= hi   AND   ball.upper() >= lo
/// ```
///
/// This is the robust form (it never false-fails a sound ball: a tighter
/// `ball superset of [lo, hi]` would, when `f(witness)` sits within Arb's
/// sub-ULP half-width of a ball edge), and it is the same predicate the
/// self-consistency lane uses (`property_ftia::contains_bracket`), so the
/// two lanes test the identical FTIA claim with different oracles. The
/// REVERSE direction (`ball subset of [lo, hi]`) is the false backstop: it
/// would pass exactly when the ball is too narrow to be sound, so it is
/// never the assertion.
pub fn contains_bracket(ball: &Ball<BigFloat>, lo: &BigFloat, hi: &BigFloat) -> bool {
    ball.lower().partial_cmp(hi).0 != Some(Ordering::Greater)
        && ball.upper().partial_cmp(lo).0 != Some(Ordering::Less)
}

/// The SOUND independent check for INTERVAL input (`BRACKETI`): the result
/// ball must be a SUPERSET of Arb's rigorous enclosure `[lo, hi]` of `f` over
/// the whole input interval. Since `f(input_interval) ⊆ [lo, hi]` (Arb is
/// rigorous) and range soundness requires `f(input_interval) ⊆ ball`, a ball
/// that contains `[lo, hi]` certainly contains the image:
///
/// ```text
/// ball.lower() <= lo   AND   ball.upper() >= hi
/// ```
///
/// This is the OPPOSITE direction from [`contains_bracket`], and deliberately
/// so. For POINT input the oracle brackets `f` at one value and the sound
/// check is overlap; a superset check there would false-fail when the true
/// value sits within Arb's sub-ULP half-width of a ball edge. For INTERVAL
/// input the oracle brackets the whole IMAGE, and overlap becomes UNSOUND: a
/// ball that misses an interior extremum still overlaps the image enclosure,
/// so only the superset direction proves range soundness.
///
/// Caveat (the reason this is asserted narrowly): Arb's interval image carries
/// an outward overshoot, because the input ball's radius is an inflated ~30-bit
/// `mag` and a steep `f` propagates that overshoot to the output. A correct
/// result ball can therefore be TIGHTER than `[lo, hi]` away from the extrema
/// (measured: at p=113 the great majority of general balls have `ball ⊉
/// [lo, hi]` for exactly this reason). So the hard superset assertion is used
/// only at an extremum straddle, where `|f'| -> 0` makes `[lo, hi]` tight while
/// the ball stays Lipschitz-wide; the general width relationship is MEASURED by
/// the tightness lane, not asserted.
pub fn contains_interval(ball: &Ball<BigFloat>, lo: &BigFloat, hi: &BigFloat) -> bool {
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
