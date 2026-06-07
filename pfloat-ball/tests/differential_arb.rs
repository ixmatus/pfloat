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
