//! Independent Arb containment lane for pfloat-ball (ADR-0078 follow-up,
//! pf-fe5f). Unlike `property_ftia.rs` (self-consistency: the same kernel
//! computes the midpoint and the witness oracle), this lane brackets the
//! witness with an INDEPENDENT Arb rigorous interval and asserts the ball
//! contains it. Reaches Arb out of process via the python-flint worker's
//! `BRACKET` verb; no FLINT/Arb in the link graph.
//!
//! Per-release / env-gated: it needs the worker venv
//! (`scripts/setup_arb_oracle.sh`). The codec round-trip test runs
//! without the venv; the worker tests skip when the venv is absent.
//!
//! This slice (pf-fe5f.3 prerequisite work) lands the subprocess driver +
//! dyadic codec and proves the pipe; the containment + tightness lanes
//! build on it.

#![cfg(feature = "differential-arb")]

mod common;

use common::arb_bracket::{
    bigfloat_to_dyadic, encode_decode, venv_available, ArbBracketWorker, Bracket,
};
use common::bf;
use pfloat::{BigFloat, RoundingMode};

const NE: RoundingMode = RoundingMode::NearestEven;

/// `n / d` at precision `p`, exactly the dyadic the test means to send.
fn ratio(n: i64, d: i64, p: u32) -> BigFloat {
    bf(n, p).div(&bf(d, p), NE).0
}

#[test]
fn dyadic_codec_round_trips_exactly() {
    // A spread of finite values: integers, dyadic fractions, both signs,
    // large and small binary exponents, at several precisions. The codec
    // is lossless, so encode-then-decode must equal the input bit-exactly.
    let cases = [
        bf(0, 53),
        bf(1, 53),
        bf(-1, 113),
        bf(1_000_003, 53),
        ratio(3, 2, 53),
        ratio(-7, 4, 113),
        ratio(1, 3, 200), // 1/3 is irrational in binary -> full p-bit mantissa
        bf(5, 24).scale_by_pow2(100).0,
        bf(5, 24).scale_by_pow2(-100).0,
        ratio(355, 113, 256),
    ];
    for x in cases {
        let rt = encode_decode(&x);
        assert_eq!(
            x.partial_cmp(&rt).0,
            Some(core::cmp::Ordering::Equal),
            "codec round-trip changed the value: {x} != {rt}"
        );
        // The dyadic encoding itself parses as expected (sane tokens).
        let (s, m, e) = bigfloat_to_dyadic(&x).expect("finite");
        assert!(s == "+" || s == "-");
        assert!(!m.is_empty());
        let _ = e;
    }
}

#[test]
fn bracket_pipe_brackets_the_reference_value() {
    if !venv_available() {
        eprintln!("skip: Arb venv absent (run scripts/setup_arb_oracle.sh)");
        return;
    }
    let mut w = ArbBracketWorker::spawn();
    // oracle_prec well above the 400-bit reference's own rounding so the
    // bracket comfortably contains the pfloat reference; this slice only
    // proves the pipe (rigor + independence are S1's and S3's tests).
    let prec = 128;
    let refp = 400;

    // (fn_id, input at p=53, reference computed by pfloat at p=400).
    let unary: [(&str, BigFloat, BigFloat); 3] = [
        ("exp", ratio(3, 2, 53), ratio(3, 2, refp).exp(NE).0),
        ("sqrt", bf(2, 53), bf(2, refp).sqrt(NE).0),
        ("ln", bf(2, 53), bf(2, refp).ln(NE).0),
    ];
    for (fn_id, x, reference) in unary {
        match w.bracket(fn_id, prec, &x, None) {
            Bracket::Finite { lo, hi } => {
                assert!(
                    lo.partial_cmp(&reference).0 != Some(core::cmp::Ordering::Greater)
                        && hi.partial_cmp(&reference).0 != Some(core::cmp::Ordering::Less),
                    "{fn_id}: reference {reference} not in bracket [{lo}, {hi}]"
                );
            }
            other => panic!("{fn_id}: expected Finite bracket, got {other:?}"),
        }
    }

    // A binary op through the pipe: add(1.5, 0.25) = 1.75.
    let a = ratio(3, 2, 53);
    let b = ratio(1, 4, 53);
    let reference = a.add(&b, NE).0.round_to_precision(refp, NE).unwrap().0;
    match w.bracket("add", prec, &a, Some(&b)) {
        Bracket::Finite { lo, hi } => assert!(
            lo.partial_cmp(&reference).0 != Some(core::cmp::Ordering::Greater)
                && hi.partial_cmp(&reference).0 != Some(core::cmp::Ordering::Less),
            "add: reference {reference} not in bracket [{lo}, {hi}]"
        ),
        other => panic!("add: expected Finite bracket, got {other:?}"),
    }
}

// ---------- S3: the unary containment lane (pf-fe5f.4) ----------
//
// The independent soundness backstop: for each ball op and each witness
// inside the input ball, the result ball must not be provably disjoint
// from Arb's rigorous bracket of f(witness). This is what lifts the crate
// from self-consistency to independently-verified soundness.
//
// Tightness (ADR-0078's secondary "logs tightness per bucket") is NOT in
// this slice: a meaningful slack needs Arb's enclosure of f over the whole
// input INTERVAL (the witness-image span is rounding-dominated for
// near-exact inputs and so measures nothing), which is a BRACKET
// interval-input extension. Filed as a follow-up; the soundness backstop,
// the load-bearing claim, lands here.

use common::arb_bracket::contains_bracket;
use common::{random_ball, witnesses, Rng};
use pfloat::Parts;
use pfloat_ball::{Ball, Mag};

/// Always-defined or simple-domain unary functions. The inverse-trig
/// domain edges (asin/acos/acosh/atanh near their boundaries) are the S5
/// edge lane's job; here a witness that lands out of domain yields a `Nan`
/// bracket and is skipped.
const UNARY: &[&str] = &[
    "exp", "expm1", "exp2", "exp10", "sinh", "cosh", "tanh", "atan", "asinh", "cbrt", "sin", "cos",
    "tan", "sqrt", "ln", "log2", "log10", "log1p",
];

/// Dispatch a ball unary op by name.
fn ball_unary(a: &Ball<BigFloat>, fn_id: &str) -> Ball<BigFloat> {
    match fn_id {
        "exp" => a.exp().0,
        "expm1" => a.expm1().0,
        "exp2" => a.exp2().0,
        "exp10" => a.exp10().0,
        "sinh" => a.sinh().0,
        "cosh" => a.cosh().0,
        "tanh" => a.tanh().0,
        "atan" => a.atan().0,
        "asinh" => a.asinh().0,
        "cbrt" => a.cbrt().0,
        "sin" => a.sin().0,
        "cos" => a.cos().0,
        "tan" => a.tan().0,
        "sqrt" => a.sqrt().0,
        "ln" => a.ln().0,
        "log2" => a.log2().0,
        "log10" => a.log10().0,
        "log1p" => a.log1p().0,
        other => panic!("unknown ball unary {other}"),
    }
}

/// Binary exponent (floor(log2|v|)) of a finite non-zero `BigFloat`,
/// `i64::MIN` for zero / non-finite.
fn exponent_of(v: &BigFloat) -> i64 {
    match v.parts() {
        Parts::Normal { exponent, .. } => exponent,
        _ => i64::MIN,
    }
}

/// `lo <= x <= hi`.
fn between(lo: &BigFloat, x: &BigFloat, hi: &BigFloat) -> bool {
    lo.partial_cmp(x).0 != Some(core::cmp::Ordering::Greater)
        && hi.partial_cmp(x).0 != Some(core::cmp::Ordering::Less)
}

#[test]
fn arb_containment_unary() {
    if !venv_available() {
        eprintln!("skip: Arb venv absent (run scripts/setup_arb_oracle.sh)");
        return;
    }
    let mut w = ArbBracketWorker::spawn();
    let deep = std::env::var("PFLOAT_DEEP").is_ok();
    let balls_per_fn = if deep { 400 } else { 40 };

    for &fn_id in UNARY {
        let seed = fn_id.bytes().fold(0xA5B6_C7D8u64, |a, b| {
            a.wrapping_mul(131).wrapping_add(b as u64)
        });
        let mut rng = Rng(seed);
        let mut usable = 0u32;
        for _ in 0..balls_per_fn {
            let p = [24u32, 53, 113][(rng.next() % 3) as usize];
            let a = random_ball(&mut rng, p);
            let result = ball_unary(&a, fn_id);
            if result.is_entire() {
                continue; // pole-straddling / overflow: encloses everything
            }
            let prec = p + 128;
            let mut any = false;
            for wit in witnesses(&a, 400) {
                if let Bracket::Finite { lo, hi } = w.bracket(fn_id, prec, &wit, None) {
                    // SOUNDNESS: Arb rigorously brackets f(wit) by [lo, hi];
                    // FTIA requires f(wit) in result; so result must not lie
                    // entirely outside [lo, hi].
                    assert!(
                        contains_bracket(&result, &lo, &hi),
                        "{fn_id} p={p}: ball [{}, {}] is disjoint from Arb's [{}, {}] for f(witness) -- UNSOUND",
                        result.lower(),
                        result.upper(),
                        lo,
                        hi
                    );
                    any = true;
                }
            }
            if any {
                usable += 1;
            }
        }
        assert!(
            usable >= 5,
            "{fn_id}: only {usable} in-domain samples (coverage gap)"
        );
    }
}

#[test]
fn witness_inside_invariant() {
    // The exact witness reconstruction must keep every witness inside the
    // input ball's denoted interval; otherwise a real soundness bug could
    // hide behind a witness that was never in the ball, turning a violation
    // into a vacuous pass.
    let mut rng = Rng(0x1357_9bdf_2468_ace0);
    for _ in 0..2000 {
        let p = [24u32, 53, 113][(rng.next() % 3) as usize];
        let a = random_ball(&mut rng, p);
        if a.is_entire() {
            continue;
        }
        let lo = a.lower();
        let hi = a.upper();
        for wit in witnesses(&a, 400) {
            assert!(
                between(&lo, &wit, &hi),
                "witness {wit} outside ball [{lo}, {hi}]"
            );
        }
    }
}

#[test]
fn negative_control_too_narrow_ball_is_caught() {
    // A deliberately too-narrow ball MUST be detected as unsound, proving
    // the sound direction has teeth (otherwise the backstop is vacuous).
    // Build the correct result, quarter its radius, and assert at least one
    // witness's Arb bracket is provably outside the narrowed ball.
    if !venv_available() {
        eprintln!("skip: Arb venv absent");
        return;
    }
    let mut w = ArbBracketWorker::spawn();
    let mut rng = Rng(0x0fed_cba9_8765_4321);
    let mut proved = 0u32;
    for _ in 0..200 {
        let p = [53u32, 113][(rng.next() % 2) as usize];
        let a = random_ball(&mut rng, p);
        let result = a.exp().0;
        let rad = match result.radius() {
            Mag::Finite { .. } => result.radius(),
            _ => continue, // exact result: nothing to narrow
        };
        // Quarter the radius -> the narrowed ball excludes the image near
        // the original edges.
        let narrow_rad = Mag::from_pow2(exponent_of(&rad.to_bigfloat()).saturating_sub(2));
        let bad = match Ball::new(result.midpoint().clone(), narrow_rad) {
            Ok(b) if !b.is_entire() => b,
            _ => continue,
        };
        for wit in witnesses(&a, 400) {
            if let Bracket::Finite { lo, hi } = w.bracket("exp", p + 128, &wit, None) {
                if !contains_bracket(&bad, &lo, &hi) {
                    proved += 1;
                    break;
                }
            }
        }
    }
    assert!(
        proved >= 10,
        "negative control proved unsoundness on only {proved} narrowed balls; the check may lack teeth"
    );
}

// ---------- S4: the binary containment lane (pf-fe5f.5) ----------

const BINARY: &[&str] = &["add", "sub", "mul", "div", "hypot", "atan2"];

fn ball_binary(a: &Ball<BigFloat>, b: &Ball<BigFloat>, fn_id: &str) -> Ball<BigFloat> {
    match fn_id {
        "add" => a.add(b).0,
        "sub" => a.sub(b).0,
        "mul" => a.mul(b).0,
        "div" => a.div(b).0,
        "hypot" => a.hypot(b).0,
        "atan2" => a.atan2(b).0,
        other => panic!("unknown ball binary {other}"),
    }
}

#[test]
fn arb_containment_binary() {
    if !venv_available() {
        eprintln!("skip: Arb venv absent (run scripts/setup_arb_oracle.sh)");
        return;
    }
    let mut w = ArbBracketWorker::spawn();
    let deep = std::env::var("PFLOAT_DEEP").is_ok();
    let pairs_per_fn = if deep { 200 } else { 25 };

    for &fn_id in BINARY {
        let seed = fn_id.bytes().fold(0xBE11_2233u64, |a, b| {
            a.wrapping_mul(137).wrapping_add(b as u64)
        });
        let mut rng = Rng(seed);
        let mut usable = 0u32;
        for _ in 0..pairs_per_fn {
            let p = [24u32, 53, 113][(rng.next() % 3) as usize];
            let a = random_ball(&mut rng, p);
            let b = random_ball(&mut rng, p);
            let result = ball_binary(&a, &b, fn_id);
            if result.is_entire() {
                continue;
            }
            let prec = p + 128;
            let mut any = false;
            for wx in witnesses(&a, 400) {
                for wy in witnesses(&b, 400) {
                    if fn_id == "div" && wy.is_zero() {
                        continue; // x/0 is the entire/flag case, not a bracket
                    }
                    // atan2(0,0), div edges, etc. give an Arb NaN / inf
                    // bracket (no finite interval to enclose) and are skipped.
                    if let Bracket::Finite { lo, hi } = w.bracket(fn_id, prec, &wx, Some(&wy)) {
                        assert!(
                            contains_bracket(&result, &lo, &hi),
                            "{fn_id} p={p}: ball [{}, {}] disjoint from Arb's [{}, {}] for f(wx,wy) -- UNSOUND",
                            result.lower(), result.upper(), lo, hi
                        );
                        any = true;
                    }
                }
            }
            if any {
                usable += 1;
            }
        }
        assert!(
            usable >= 5,
            "{fn_id}: only {usable} usable pairs (coverage gap)"
        );
    }
}

#[test]
fn div_by_zero_straddling_ball_is_entire() {
    // A divisor ball straddling zero makes the quotient unbounded; the ball
    // div must report `entire` (it encloses everything), so the containment
    // lane's `result.is_entire()` skip is exercised, not a spurious bracket.
    let p = 53;
    let num = random_ball(&mut Rng(7), p);
    // A ball centered at 0 with positive radius straddles zero.
    let denom = Ball::new(bf(0, p), Mag::from_pow2(-3)).unwrap();
    let (q, _) = num.div(&denom);
    assert!(
        q.is_entire(),
        "div by a zero-straddling ball must be entire"
    );
}

// ---------- S5: the domain-edge lane (pf-fe5f.6) ----------
//
// The inverse-trig functions S3 deferred (asin/acos near ±1, acosh near 1,
// atanh near ±1), plus the reconciliation the BRACKET verb surfaced: Arb
// returns a NaN ball at a point outside the domain, while pfloat-ball uses
// a Status::INVALID + (entire | sound-over-the-in-domain-part) convention.
// The lane checks both: over the in-domain witnesses the ball still
// contains Arb's bracket (soundness holds on the defined part), and at an
// out-of-domain point Arb's NaN lines up with the ball's INVALID flag.

use pfloat::Status;

fn ball_unary_status(a: &Ball<BigFloat>, fn_id: &str) -> (Ball<BigFloat>, Status) {
    match fn_id {
        "asin" => a.asin(),
        "acos" => a.acos(),
        "acosh" => a.acosh(),
        "atanh" => a.atanh(),
        other => panic!("unknown edge fn {other}"),
    }
}

/// A ball near the domain boundary of `fn_id`: some entirely in domain,
/// some straddling the edge (so the INVALID / entire convention is
/// exercised). Midpoints are exact dyadics `m * 2^-20`.
fn edge_ball(rng: &mut Rng, p: u32, fn_id: &str) -> Ball<BigFloat> {
    let radexp = -(10 + (rng.next() % 24) as i64);
    let rad = Mag::from_pow2(radexp);
    // Half the balls are constructed to STRADDLE a boundary so the INVALID
    // / clamp reconciliation path is actually exercised; the rest spread
    // across the in-domain and out-of-domain sides. (A uniform spread on a
    // 2^-20 grid with a <= 2^-10 radius straddles a boundary with
    // probability ~1/1000, so at the default sample count it never did --
    // the pre-merge review's D5 finding.)
    let mid = if rng.next() % 2 == 0 {
        // Snap to a boundary plus a sub-radius offset (|offset| <= 0.75*rad
        // < rad), so the boundary lies strictly inside [mid - rad, mid + rad]
        // and the ball provably straddles it.
        let b = match fn_id {
            "acosh" => bf(1, p), // domain [1, inf): 1 is the only boundary
            _ if rng.next() % 2 == 0 => bf(1, p), // asin / acos / atanh: +-1
            _ => bf(-1, p),
        };
        let (offset, _) = bf(rng.int(3), p).scale_by_pow2(radexp - 2);
        b.add(&offset, NE).0
    } else {
        match fn_id {
            "acosh" => {
                // [0.8, 4): below 1 is out of domain, at/above 1 in domain.
                let lo = (8i64 << 20) / 10;
                let span = (4u64 << 20) - lo as u64;
                bf(lo + (rng.next() % span) as i64, p).scale_by_pow2(-20).0
            }
            _ => {
                // asin / acos / atanh spread over ~[-1.2, 1.2].
                let span = (12i64 << 20) / 10;
                bf(rng.int(span), p).scale_by_pow2(-20).0
            }
        }
    };
    Ball::new(mid, rad).unwrap()
}

/// Whether the ball's denoted interval has `b` strictly interior -- a
/// domain-boundary straddle (one side in domain, the other out), the case
/// that exercises pfloat-ball's clamp + INVALID reconciliation.
fn straddles(a: &Ball<BigFloat>, b: &BigFloat) -> bool {
    a.lower().partial_cmp(b).0 == Some(core::cmp::Ordering::Less)
        && a.upper().partial_cmp(b).0 == Some(core::cmp::Ordering::Greater)
}

#[test]
fn arb_containment_edge() {
    if !venv_available() {
        eprintln!("skip: Arb venv absent (run scripts/setup_arb_oracle.sh)");
        return;
    }
    let mut w = ArbBracketWorker::spawn();
    let deep = std::env::var("PFLOAT_DEEP").is_ok();
    let n = if deep { 200 } else { 30 };

    for &fn_id in &["asin", "acos", "acosh", "atanh"] {
        let seed = fn_id.bytes().fold(0xED6E_0000u64, |a, b| {
            a.wrapping_mul(139).wrapping_add(b as u64)
        });
        let mut rng = Rng(seed);
        let mut usable = 0u32;
        let mut invalid = 0u32;
        let mut straddled = 0u32;
        let mut clamp_verified = 0u32;
        let one = bf(1, 53);
        let neg_one = bf(-1, 53);
        let bounds: &[&BigFloat] = if fn_id == "acosh" {
            &[&one]
        } else {
            &[&one, &neg_one]
        };
        for _ in 0..n {
            let p = [24u32, 53, 113][(rng.next() % 3) as usize];
            let a = edge_ball(&mut rng, p, fn_id);
            // A ball whose interval contains a boundary exercises the
            // clamp + INVALID reconciliation; count it as the real measure
            // (INVALID alone also fires for fully-out-of-domain balls).
            let is_straddle = bounds.iter().any(|b| straddles(&a, b));
            if is_straddle {
                straddled += 1;
            }
            let (result, status) = ball_unary_status(&a, fn_id);
            if status.invalid() {
                invalid += 1;
            }
            if result.is_entire() {
                continue;
            }
            let mut any = false;
            for wit in witnesses(&a, 400) {
                // An out-of-domain witness yields a NaN bracket (Arb agrees
                // the function is undefined there) and is skipped; over the
                // in-domain witnesses the ball must still be sound.
                if let Bracket::Finite { lo, hi } = w.bracket(fn_id, p + 128, &wit, None) {
                    assert!(
                        contains_bracket(&result, &lo, &hi),
                        "{fn_id} edge p={p}: ball [{}, {}] disjoint from Arb's [{}, {}] over the in-domain part -- UNSOUND",
                        result.lower(), result.upper(), lo, hi
                    );
                    any = true;
                }
            }
            if any {
                usable += 1;
                // A straddler with a finite (clamped) result and an
                // in-domain witness independently confirms the clamp is
                // SOUND, not merely that INVALID was flagged.
                if is_straddle {
                    clamp_verified += 1;
                }
            }
        }
        assert!(
            usable >= 3,
            "{fn_id}: only {usable} edge balls had in-domain witnesses"
        );
        assert!(
            straddled >= 3,
            "{fn_id}: only {straddled} boundary-straddling balls; the clamp / INVALID reconciliation path is under-exercised"
        );
        assert!(
            invalid >= 1,
            "{fn_id}: never raised INVALID on a boundary ball"
        );
        // asin / acos / acosh clamp the out-of-domain part to the boundary
        // value and return a finite ball, so the clamp's soundness is
        // independently checkable; atanh straddles a pole and correctly
        // goes entire, so it has no finite clamp to verify.
        if fn_id != "atanh" {
            assert!(
                clamp_verified >= 1,
                "{fn_id}: clamp-to-finite soundness never independently verified on a straddling ball"
            );
        }
    }
}

#[test]
fn domain_edge_invalid_matches_arb_nan() {
    // Reconciliation: at a point outside the domain, pfloat-ball raises
    // INVALID and Arb returns a NaN ball. Pin a few.
    let p = 53;
    let cases: [(&str, BigFloat); 3] = [
        ("asin", ratio(3, 2, p)),  // 1.5 > 1
        ("atanh", ratio(3, 2, p)), // 1.5 > 1
        ("acosh", ratio(1, 2, p)), // 0.5 < 1
    ];
    let mut worker = if venv_available() {
        Some(ArbBracketWorker::spawn())
    } else {
        None
    };
    for (fn_id, x) in cases {
        let a = Ball::point(x.clone()).unwrap();
        let (_, status) = ball_unary_status(&a, fn_id);
        assert!(status.invalid(), "{fn_id} out of domain must raise INVALID");
        if let Some(w) = worker.as_mut() {
            assert!(
                matches!(w.bracket(fn_id, p + 128, &x, None), Bracket::Nan),
                "{fn_id} out of domain: Arb should return a NaN ball"
            );
        }
    }
}
