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
