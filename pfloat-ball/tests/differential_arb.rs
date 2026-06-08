//! Independent Arb containment and range-soundness lanes for pfloat-ball
//! (ADR-0078, extended by ADR-0082; pf-fe5f and pf-fe5f.7). Where
//! `property_ftia.rs` is self-consistency (the same kernel computes the
//! midpoint and the witness oracle), the lanes here check each ball against
//! an INDEPENDENT Arb rigorous enclosure reached out of process through the
//! python-flint worker, so no FLINT/Arb enters the link graph.
//!
//! Two predicates run, against two worker verbs:
//!
//! - Point containment (`BRACKET`): bracket the true value at sampled
//!   witnesses inside the input ball and assert the ball contains the
//!   bracket. This is overlap: the ball must admit Arb's enclosure of the
//!   true value.
//! - Interval range soundness (`BRACKETI`): bracket `f` over the WHOLE
//!   input interval and assert the result ball is a SUPERSET of that image.
//!   Asserted only at an extremum straddle, where it is provably clean; off
//!   the extrema a correct ball can be tighter than Arb's inflated interval
//!   image, so the general width relationship is measured, not asserted, by
//!   the per-bucket tightness regression table.
//!
//! Per-release / env-gated: the worker lanes need the worker venv
//! (`scripts/setup_arb_oracle.sh`), and the tightness table needs the pinned
//! Arb build, so it is `#[ignore]`. The codec round-trip test runs without
//! the venv; the worker tests skip when the venv is absent.

#![cfg(feature = "differential-arb")]

mod common;

use common::arb_bracket::{
    arb_lane_available, bigfloat_to_dyadic, encode_decode, venv_available, ArbBracketWorker,
    Bracket,
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
    if !arb_lane_available("bracket_pipe_brackets_the_reference_value") {
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

/// Two brackets denote the same verdict: identical finite endpoints (exact
/// dyadics, so bit-for-bit), or the same non-finite class.
fn brackets_equal(a: &Bracket, b: &Bracket) -> bool {
    use core::cmp::Ordering::Equal;
    match (a, b) {
        (Bracket::Finite { lo: la, hi: ha }, Bracket::Finite { lo: lb, hi: hb }) => {
            la.partial_cmp(lb).0 == Some(Equal) && ha.partial_cmp(hb).0 == Some(Equal)
        }
        (Bracket::Nan, Bracket::Nan)
        | (Bracket::PosInf, Bracket::PosInf)
        | (Bracket::NegInf, Bracket::NegInf)
        | (Bracket::Inconclusive, Bracket::Inconclusive) => true,
        _ => false,
    }
}

#[test]
fn degenerate_interval_bracket_equals_point_bracket() {
    // S1 (pf-fe5f.7): a BRACKETI call with a zero radius denotes the exact
    // point [mid, mid], so it must reduce to the point BRACKET bit-for-bit.
    // This pins the interval verb against the established point verb and is
    // the base case the range-soundness lane builds on -- a degenerate ball
    // is a point, and its image is f(mid). `union(p, p) = p` for an exact `p`,
    // so the two requests reach `dispatch_elementary` with the same operand.
    if !arb_lane_available("degenerate_interval_bracket_equals_point_bracket") {
        return;
    }
    let mut w = ArbBracketWorker::spawn();
    let prec = 256;
    let zero = bf(0, 53);

    // Unary spread: includes cbrt's negative half (the real-odd path, where
    // Arb's principal root and pfloat-ball's real root diverge) and the
    // base-changed / hyperbolic surface.
    let unary: [(&str, BigFloat); 6] = [
        ("exp", ratio(3, 2, 53)),
        ("ln", bf(2, 53)),
        ("sin", ratio(11, 8, 53)),
        ("cbrt", bf(-8, 53)),
        ("tanh", ratio(-5, 4, 53)),
        ("log2", ratio(7, 2, 53)),
    ];
    for (fn_id, x) in unary {
        let point = w.bracket(fn_id, prec, &x, None);
        let interval = w.bracket_interval(fn_id, prec, &x, &zero, None);
        assert!(
            brackets_equal(&point, &interval),
            "{fn_id}: degenerate interval {interval:?} != point {point:?}"
        );
    }

    // Binary: degenerate intervals on both operands.
    let a = ratio(3, 2, 53);
    let b = ratio(1, 4, 53);
    for fn_id in ["add", "sub", "mul", "div", "hypot", "atan2"] {
        let point = w.bracket(fn_id, prec, &a, Some(&b));
        let interval = w.bracket_interval(fn_id, prec, &a, &zero, Some((&b, &zero)));
        assert!(
            brackets_equal(&point, &interval),
            "{fn_id}: degenerate interval {interval:?} != point {point:?}"
        );
    }
}

// ---------- S3: the unary containment lane (pf-fe5f.4) ----------
//
// The independent soundness backstop: for each ball op and each witness
// inside the input ball, the result ball must not be provably disjoint
// from Arb's rigorous bracket of f(witness). This is what lifts the crate
// from self-consistency to independently-verified soundness.
//
// Two things need Arb's enclosure of f over the whole input INTERVAL (a
// BRACKET interval-input extension, deferred to pf-fe5f.7), not just at
// sampled witnesses. (1) Tightness (ADR-0078's secondary "logs tightness
// per bucket"): a meaningful slack is rounding-dominated at the witnesses
// for near-exact inputs and so measures nothing. (2) Range soundness for
// non-monotonic functions: five witnesses cannot see a ball that excludes
// an interior extremum of sin/cos/tan, so this lane -- like the
// self-consistency lane -- is structurally blind to a missed-extremum bug;
// Phase 1's exhaustive scalar sweep is the current guard for the kernel.
// The point-sampling soundness backstop -- the load-bearing every-witness
// containment claim -- lands here.

use common::arb_bracket::{contains_bracket, contains_interval};
use common::{bin_exponent, random_ball, witnesses, Rng};
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
    if !arb_lane_available("arb_containment_unary") {
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
        let mut checked = 0u32;
        let mut skipped = 0u32;
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
                match w.bracket(fn_id, prec, &wit, None) {
                    Bracket::Finite { lo, hi } => {
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
                        checked += 1;
                        any = true;
                    }
                    // NaN (out-of-domain witness) / INF / INC: no finite
                    // interval to enclose. Counted, not silently dropped, so a
                    // function whose witnesses mostly skip (the cbrt-negative
                    // class) is visible rather than passing on thin coverage.
                    _ => skipped += 1,
                }
            }
            if any {
                usable += 1;
            }
        }
        eprintln!(
            "arb_containment_unary {fn_id}: {usable} usable balls, {checked} witness brackets checked, {skipped} skipped (NaN/INF/INC)"
        );
        assert!(
            usable >= 5,
            "{fn_id}: only {usable} in-domain samples (coverage gap)"
        );
    }
}

/// The unary functions seeded from their own hard-to-round corpus. `cbrt` /
/// `sqrt` have no corpus, and the inverse-trig edge functions are the edge
/// lane's job, so neither is seeded here.
const HTR_UNARY: &[&str] = &[
    "exp", "expm1", "exp2", "exp10", "sinh", "cosh", "tanh", "atan", "asinh", "sin", "cos", "tan",
    "ln", "log2", "log10", "log1p",
];

#[test]
fn arb_containment_unary_hard_to_round() {
    // The independent soundness lane, seeded with the Lefevre-Muller /
    // CORE-MATH hard-to-round corpus: f(mid) is boundary-close, so the ball's
    // directed endpoints are maximally stressed against Arb's rigorous
    // bracket. Same point-witness containment claim as arb_containment_unary,
    // hardest available inputs. pf-vcqh, ADR-0078 deferral.
    if !arb_lane_available("arb_containment_unary_hard_to_round") {
        return;
    }
    let mut w = ArbBracketWorker::spawn();
    let mut rng = Rng(0x4c4d_5345_4544_4152);
    // seeded_ball's precision is 53 or 113; oracle prec above the larger.
    let prec = 113 + 128;
    for &fn_id in HTR_UNARY {
        let cases = common::lm_cases_for(fn_id).expect("seeded fn has a corpus");
        let mut usable = 0u32;
        let mut checked = 0u32;
        for &(xbits, _) in cases {
            if !common::is_finite_nonzero_f64(xbits) {
                continue;
            }
            let a = common::seeded_ball(&mut rng, xbits);
            let result = ball_unary(&a, fn_id);
            if result.is_entire() {
                continue;
            }
            let mut any = false;
            for wit in witnesses(&a, 400) {
                if let Bracket::Finite { lo, hi } = w.bracket(fn_id, prec, &wit, None) {
                    assert!(
                        contains_bracket(&result, &lo, &hi),
                        "{fn_id} HtR: ball [{}, {}] disjoint from Arb's [{}, {}] for f(witness) -- UNSOUND",
                        result.lower(),
                        result.upper(),
                        lo,
                        hi
                    );
                    checked += 1;
                    any = true;
                }
            }
            if any {
                usable += 1;
            }
        }
        eprintln!(
            "arb_containment_unary_hard_to_round {fn_id}: {usable} usable seeded balls, {checked} witness brackets checked"
        );
        assert!(
            usable >= 5,
            "{fn_id} HtR: only {usable} usable seeded balls (corpus coverage gap)"
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
    if !arb_lane_available("negative_control_too_narrow_ball_is_caught") {
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
    if !arb_lane_available("arb_containment_binary") {
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
        let mut checked = 0u32;
        let mut skipped = 0u32;
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
                    // bracket (no finite interval to enclose); counted, not
                    // silently dropped.
                    match w.bracket(fn_id, prec, &wx, Some(&wy)) {
                        Bracket::Finite { lo, hi } => {
                            assert!(
                                contains_bracket(&result, &lo, &hi),
                                "{fn_id} p={p}: ball [{}, {}] disjoint from Arb's [{}, {}] for f(wx,wy) -- UNSOUND",
                                result.lower(), result.upper(), lo, hi
                            );
                            checked += 1;
                            any = true;
                        }
                        _ => skipped += 1,
                    }
                }
            }
            if any {
                usable += 1;
            }
        }
        eprintln!(
            "arb_containment_binary {fn_id}: {usable} usable pairs, {checked} brackets checked, {skipped} skipped (NaN/INF/INC)"
        );
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
    let mid = if rng.next().is_multiple_of(2) {
        // Snap to a boundary plus a sub-radius offset (|offset| <= 0.75*rad
        // < rad), so the boundary lies strictly inside [mid - rad, mid + rad]
        // and the ball provably straddles it.
        let b = match fn_id {
            "acosh" => bf(1, p), // domain [1, inf): 1 is the only boundary
            _ if rng.next().is_multiple_of(2) => bf(1, p), // asin / acos / atanh: +-1
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
    if !arb_lane_available("arb_containment_edge") {
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

// ---------- pf-fe5f.7 S3: interval range-soundness lane ----------
//
// The point lanes above sample witnesses INSIDE the input ball, so they are
// structurally blind to a result ball that fails to enclose an interior
// extremum of a non-monotonic function: no witness reaches sin's peak at
// pi/2, but the true image does. The interval BRACKETI verb closes this: Arb's
// image of the WHOLE input interval sees the extremum, and a sound result ball
// must be a SUPERSET of that image (`contains_interval`, the opposite of the
// point predicate; see its doc comment for why overlap is unsound here).
//
// The hard superset assertion is used only at an extremum straddle. There
// |f'| -> 0 makes Arb's image tight while pfloat-ball's Lipschitz radius stays
// wide, so a correct ball is robustly a superset (measured: 0 false-fails at
// p in {24,53,113}). Off the extrema a correct ball can be tighter than Arb's
// interval image (its input radius is an inflated ~30-bit mag), so the general
// width relationship is the tightness lane's job (S4), not a pass/fail here.

/// The interior extremum of `f` (sin/cos) at index `n`, computed at precision
/// `work`: sin extrema at `pi/2 + n*pi`, cos extrema at `n*pi`; the value is
/// `(-1)^n` (`+1` a maximum, `-1` a minimum).
fn extremum_loc(f: &str, n: i64, work: u32) -> BigFloat {
    let pi = pfloat::constants::pi(work, NE).0;
    let npi = pi.mul(&bf(n, work), NE).0;
    match f {
        "sin" => pfloat::constants::pi_over_2(work, NE).0.add(&npi, NE).0,
        "cos" => npi,
        other => panic!("extremum_loc only for sin/cos, got {other}"),
    }
}

/// A ball straddling an interior extremum of `f`: the midpoint is the extremum
/// location plus a sub-radius offset, so the extremum lies strictly inside
/// `[mid - rad, mid + rad]` yet is reached by no endpoint witness. Returns the
/// ball and the high-precision extremum location (for the straddle check).
fn extremum_ball(rng: &mut Rng, f: &str, p: u32) -> (Ball<BigFloat>, BigFloat) {
    let work = p + 80;
    let n = rng.int(3); // a spread of maxima and minima at modest |location|
    let loc = extremum_loc(f, n, work);
    let rad_exp = -6 - (rng.next() % 10) as i64; // radius 2^-6 .. 2^-15
    let off_sign = if rng.next().is_multiple_of(2) { 1 } else { -1 };
    // |offset| = 2^(rad_exp - 2) is a quarter of the radius, so the extremum
    // sits strictly interior but away from the midpoint.
    let (offset, _) = bf(off_sign, p).scale_by_pow2(rad_exp - 2);
    let mid = loc.round_to_precision(p, NE).unwrap().0.add(&offset, NE).0;
    (Ball::new(mid, Mag::from_pow2(rad_exp)).unwrap(), loc)
}

/// Whether the high-precision `loc` is strictly interior to the ball.
fn straddles_loc(a: &Ball<BigFloat>, loc: &BigFloat) -> bool {
    a.lower().partial_cmp(loc).0 == Some(core::cmp::Ordering::Less)
        && a.upper().partial_cmp(loc).0 == Some(core::cmp::Ordering::Greater)
}

#[test]
fn arb_range_soundness_extrema() {
    if !arb_lane_available("arb_range_soundness_extrema") {
        return;
    }
    let mut w = ArbBracketWorker::spawn();
    for &f in &["sin", "cos"] {
        let seed = f.bytes().fold(0x5241_4e47_0000_0001u64, |a, b| {
            a.wrapping_mul(131).wrapping_add(b as u64)
        });
        for &p in &[24u32, 53, 113] {
            let mut rng = Rng(seed ^ p as u64);
            let mut checked = 0u32;
            let mut straddled = 0u32;
            for _ in 0..40 {
                let (a, loc) = extremum_ball(&mut rng, f, p);
                if !straddles_loc(&a, &loc) {
                    continue; // not an actual straddle (rare at small radius)
                }
                straddled += 1;
                let result = ball_unary(&a, f);
                if result.is_entire() {
                    continue;
                }
                let rad_bf = a.radius().to_bigfloat();
                if let Bracket::Finite { lo, hi } =
                    w.bracket_interval(f, p + 128, a.midpoint(), &rad_bf, None)
                {
                    // SOUNDNESS: Arb's image [lo, hi] encloses f over the whole
                    // straddle (extremum included); a sound ball must be a
                    // superset of it. A ball that stops short -- a missed
                    // interior extremum -- is caught here.
                    assert!(
                        contains_interval(&result, &lo, &hi),
                        "{f} p={p}: ball [{}, {}] is NOT a superset of Arb's interval image [{}, {}] over an extremum straddle -- a missed interior extremum",
                        result.lower(),
                        result.upper(),
                        lo,
                        hi
                    );
                    checked += 1;
                }
            }
            eprintln!(
                "arb_range_soundness_extrema {f} p={p}: {checked} interval brackets checked, {straddled} straddles"
            );
            assert!(
                straddled >= 10,
                "{f} p={p}: only {straddled} extremum straddles generated"
            );
            assert!(
                checked >= 10,
                "{f} p={p}: only {checked} interval brackets checked"
            );
        }
    }
}

#[test]
fn range_soundness_catches_missed_extremum() {
    // Negative control with teeth: the naive endpoint-only enclosure (treating
    // sin / cos as if monotonic) misses the interior extremum, since both
    // endpoints sit strictly past the peak. The superset check MUST reject it,
    // or the lane is vacuous.
    if !arb_lane_available("range_soundness_catches_missed_extremum") {
        return;
    }
    let mut w = ArbBracketWorker::spawn();
    let mut caught = 0u32;
    for &f in &["sin", "cos"] {
        let seed = f.bytes().fold(0x4d49_5353_0000_0001u64, |a, b| {
            a.wrapping_mul(137).wrapping_add(b as u64)
        });
        for &p in &[24u32, 53, 113] {
            let mut rng = Rng(seed ^ p as u64);
            for _ in 0..40 {
                let (a, loc) = extremum_ball(&mut rng, f, p);
                if !straddles_loc(&a, &loc) {
                    continue;
                }
                // The endpoint-only enclosure: [min(f(alo), f(ahi)), max(...)],
                // outward-rounded. For a straddle both endpoints lie past the
                // extremum, so this ball never reaches the peak value.
                let alo = a.lower();
                let ahi = a.upper();
                let slo = scalar_sin_cos(f, &alo);
                let shi = scalar_sin_cos(f, &ahi);
                let (lo_end, hi_end) =
                    if slo.partial_cmp(&shi).0 == Some(core::cmp::Ordering::Greater) {
                        (shi, slo)
                    } else {
                        (slo, shi)
                    };
                let naive = Ball::from_interval(&lo_end, &hi_end).unwrap();
                let rad_bf = a.radius().to_bigfloat();
                if let Bracket::Finite { lo, hi } =
                    w.bracket_interval(f, p + 128, a.midpoint(), &rad_bf, None)
                {
                    if !contains_interval(&naive, &lo, &hi) {
                        caught += 1;
                    }
                }
            }
        }
    }
    assert!(
        caught >= 20,
        "missed-extremum negative control caught only {caught}; the range check lacks teeth"
    );
}

/// The scalar `sin` / `cos` kernel at NearestEven, for the endpoint-only
/// enclosure the negative control builds.
fn scalar_sin_cos(f: &str, x: &BigFloat) -> BigFloat {
    match f {
        "sin" => x.sin(NE).0,
        "cos" => x.cos(NE).0,
        other => panic!("scalar_sin_cos only for sin/cos, got {other}"),
    }
}

// ---------- pf-fe5f.7 S4: interval tightness regression table ----------
//
// The slack `log2(width(ball) / width(arb_interval_image))` measures how much
// looser pfloat-ball's enclosure is than Arb's rigorous image of the input
// interval. The witness-image span measures nothing here (it is rounding-
// dominated for near-exact inputs), which is exactly why the interval bracket
// is needed. The slack is a per-(function, precision, magnitude) DISTRIBUTION,
// recorded as an expected COUNT PER BUCKET in a checked-in table -- never an
// aggregate floor, which would let one bucket's regression hide behind
// another's improvement. A change to the enclosure algorithm shifts the
// counts; the table is then regenerated deliberately.
//
// The widths come from Arb (the image width is the worker's), so the table is
// python-flint / Arb version-sensitive: regenerate under the pinned venv with
//   PFLOAT_ARB_TIGHTNESS_REGEN=1 cargo test -p pfloat-ball \
//     --features differential-arb arb_tightness_regression -- --ignored
// (the lane is `#[ignore]` so the per-push run, which has no pinned oracle,
// does not depend on the exact Arb build; the per-release run includes it).

/// The functions whose tightness is tracked: monotonic (`exp`, `ln`, `atan`),
/// Lipschitz (`sin`, `cos`), and composed (`tan`) cover the three enclosure
/// shapes.
const TIGHTNESS_FNS: &[&str] = &["exp", "ln", "atan", "sin", "cos", "tan"];

/// A "fat" ball: the radius is `2^-4 .. 2^-12` relative to the midpoint, well
/// above the `p`-bit rounding floor, so the widths reflect `f`'s variation
/// rather than the kernel's rounding (the near-exact regime the bead excludes).
fn tightness_ball(rng: &mut Rng, p: u32) -> Ball<BigFloat> {
    let m = match rng.int(1 << 20) {
        0 => 1,
        v => v,
    };
    let scale = rng.int(12);
    let (mid, _) = bf(m, p).scale_by_pow2(scale);
    let e = bin_exponent(&mid);
    let radexp = e - (4 + (rng.next() % 9) as i64); // e-4 .. e-12 (relative 2^-4 .. 2^-12)
    Ball::new(mid, Mag::from_pow2(radexp)).unwrap()
}

/// Three coarse magnitude bands of the midpoint, by binary exponent.
fn mag_bucket(mid: &BigFloat) -> i64 {
    match bin_exponent(mid) {
        i64::MIN => 0,
        e if e < -8 => -1,
        e if e >= 8 => 1,
        _ => 0,
    }
}

/// `log2(width)` as a binary exponent, or `None` for a zero / non-finite width
/// (a degenerate sample the tightness metric cannot speak to).
fn width_log2(lo: &BigFloat, hi: &BigFloat) -> Option<i64> {
    let w = hi.sub(lo, NE).0;
    match bin_exponent(&w) {
        i64::MIN => None,
        e => Some(e),
    }
}

/// Build the live tightness histogram: `(fn, precision, mag_bucket, slack_band)
/// -> count`. The slack band is `clamp(log2(width_R) - log2(width_J), -2, 8)`,
/// a robust integer proxy for `floor(log2(width_R / width_J))`.
fn tightness_histogram(
    w: &mut ArbBracketWorker,
) -> std::collections::BTreeMap<(String, u32, i64, i64), u32> {
    let mut hist = std::collections::BTreeMap::new();
    for &fn_id in TIGHTNESS_FNS {
        let seed = fn_id.bytes().fold(0x7137_0000_0000_0001u64, |a, b| {
            a.wrapping_mul(149).wrapping_add(b as u64)
        });
        for &p in &[24u32, 53, 113] {
            let mut rng = Rng(seed ^ p as u64);
            for _ in 0..150 {
                let a = tightness_ball(&mut rng, p);
                let result = ball_unary(&a, fn_id);
                if result.is_entire() {
                    continue;
                }
                let rad_bf = a.radius().to_bigfloat();
                let Bracket::Finite { lo, hi } =
                    w.bracket_interval(fn_id, p + 128, a.midpoint(), &rad_bf, None)
                else {
                    continue;
                };
                let (Some(er), Some(ej)) = (
                    width_log2(&result.lower(), &result.upper()),
                    width_log2(&lo, &hi),
                ) else {
                    continue; // a degenerate width: the metric says nothing
                };
                let band = (er - ej).clamp(-2, 8);
                let mag = mag_bucket(a.midpoint());
                *hist.entry((fn_id.to_string(), p, mag, band)).or_insert(0) += 1;
            }
        }
    }
    hist
}

fn tightness_table_path() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/arb_tightness_expected.txt")
}

fn format_tightness(hist: &std::collections::BTreeMap<(String, u32, i64, i64), u32>) -> String {
    let mut out = String::from(
        "# pfloat-ball interval tightness baseline (pf-fe5f.7 S4). DO NOT hand-edit.\n\
         # slack band = clamp(log2(width_ball) - log2(width_arb_image), -2, 8).\n\
         # columns: fn precision mag_bucket slack_band count\n\
         # Arb/python-flint version-sensitive; regenerate under the pinned venv with\n\
         #   PFLOAT_ARB_TIGHTNESS_REGEN=1 cargo test -p pfloat-ball \\\n\
         #     --features differential-arb arb_tightness_regression -- --ignored --nocapture\n",
    );
    for ((fn_id, p, mag, band), count) in hist {
        out.push_str(&format!("{fn_id} {p} {mag} {band} {count}\n"));
    }
    out
}

fn parse_tightness(text: &str) -> std::collections::BTreeMap<(String, u32, i64, i64), u32> {
    let mut map = std::collections::BTreeMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        assert_eq!(f.len(), 5, "malformed tightness row: `{line}`");
        let key = (
            f[0].to_string(),
            f[1].parse().expect("precision"),
            f[2].parse().expect("mag_bucket"),
            f[3].parse().expect("slack_band"),
        );
        map.insert(key, f[4].parse().expect("count"));
    }
    map
}

#[test]
#[ignore = "per-release: needs the pinned Arb venv; the image width is version-sensitive"]
fn arb_tightness_regression() {
    if !arb_lane_available("arb_tightness_regression") {
        return;
    }
    let mut w = ArbBracketWorker::spawn();
    let live = tightness_histogram(&mut w);
    let path = tightness_table_path();

    if std::env::var("PFLOAT_ARB_TIGHTNESS_REGEN").is_ok() {
        std::fs::write(&path, format_tightness(&live)).expect("write tightness table");
        eprintln!("regenerated {} ({} buckets)", path.display(), live.len());
        return;
    }

    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!(
            "tightness table {} missing ({e}); regenerate with PFLOAT_ARB_TIGHTNESS_REGEN=1",
            path.display()
        )
    });
    let expected = parse_tightness(&text);

    // Per-bucket comparison: report every cell that differs, rather than an
    // aggregate count, so a regression in one bucket cannot hide behind a
    // compensating shift in another.
    let mut diffs = Vec::new();
    let keys: std::collections::BTreeSet<_> = live.keys().chain(expected.keys()).collect();
    for k in keys {
        let (l, e) = (
            live.get(k).copied().unwrap_or(0),
            expected.get(k).copied().unwrap_or(0),
        );
        if l != e {
            diffs.push(format!(
                "  {} p={} mag={} band={}: expected {e}, got {l}",
                k.0, k.1, k.2, k.3
            ));
        }
    }
    assert!(
        diffs.is_empty(),
        "tightness histogram drifted from the committed baseline ({} buckets differ):\n{}\nregenerate with PFLOAT_ARB_TIGHTNESS_REGEN=1 if this is an intended enclosure change",
        diffs.len(),
        diffs.join("\n")
    );
}
