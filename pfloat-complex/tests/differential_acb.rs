//! Independent complex-Arb (`acb`) componentwise certified-rounding lane
//! (ADR-0092, the C5 soundness backstop).
//!
//! The enumerated tables (`annex_g_special_values.rs`) and the identities
//! (`identities.rs`) are written by reading the standard and the kernels, so a
//! shared transcription error could make a wrong test agree with a wrong
//! kernel. This lane breaks that circularity: it computes each operation's true
//! value in an INDEPENDENT engine (python-flint's rigorous `acb` ball
//! arithmetic, reached out of process so no FLINT/Arb enters the link graph)
//! and checks pfloat-complex's output is the correct rounding, bit-for-bit, of
//! each component.
//!
//! The check is the certified-rounding test, run per component: Arb rigorously
//! encloses the true value `v in [lo, hi]`; if `round(lo, p, mode) ==
//! round(hi, p, mode)` (value AND sign) then that is the unique correct
//! rounding `cr` (rounding is monotone, so `v` rounds there too), and
//! pfloat-complex's component must equal `cr`. When the enclosure still
//! straddles a `p`-bit boundary, the oracle precision is increased; if it never
//! certifies it is a hard-to-round Ziv-cap residual (counted, not failed).
//!
//! Arb has NO signed zero, so this lane certifies only the NUMERIC value of
//! finite, nonzero components; the signed-zero branch rows and the inf/NaN
//! special values are pinned by the enumerated lane instead.
//!
//! Per-release / env-gated: needs the worker venv (`scripts/setup_arb_oracle.sh`).
//! `PFLOAT_ARB_REQUIRED=1` turns a missing venv into a hard failure so the
//! backstop cannot silently no-op. `PFLOAT_DEEP=1` widens the input grids.

#![cfg(feature = "differential-acb")]

mod common;

use common::acb_bracket::{
    acb_lane_available, encode_decode, exact_signed_eq, AcbComplexWorker, Comp,
};
use common::{bf, bf_of_f64_bits, is_finite_nonzero_f64, lm_inputs_for, Rng, ALL_MODES};
use pfloat::{BigFloat, RoundingMode};
use pfloat_complex::Complex;

const NE: RoundingMode = RoundingMode::NearestEven;
const PRECS: [u32; 3] = [24, 53, 113];

#[derive(Clone, Copy)]
enum Which {
    Re,
    Im,
}

/// Compute `fn_id(z [, w])` through the public pfloat-complex API at the input
/// components' precision.
fn compute(
    fn_id: &str,
    z: (&BigFloat, &BigFloat),
    w: Option<(&BigFloat, &BigFloat)>,
    mode: RoundingMode,
) -> Complex<BigFloat> {
    let zc = Complex::new(z.0.clone(), z.1.clone());
    match fn_id {
        "csqrt" => zc.sqrt(mode).0,
        "cexp" => zc.exp(mode).0,
        "clog" => zc.log(mode).0,
        "cmul" | "cdiv" => {
            let (wr, wi) = w.expect("binary op needs w");
            let wc = Complex::new(wr.clone(), wi.clone());
            if fn_id == "cmul" {
                zc.mul(&wc, mode).0
            } else {
                zc.div(&wc, mode).0
            }
        }
        other => panic!("unknown fn_id {other}"),
    }
}

/// The certified rounding of one component to precision `p` under `mode`, via
/// the acb oracle, tightening the oracle precision until it certifies a unique
/// value. `None` if the bracket still straddles a `p`-bit boundary at the
/// deepest oracle precision (a Ziv-cap residual) or the component is non-finite
/// (a special the enumerated lane owns).
fn certify(
    w: &mut AcbComplexWorker,
    fn_id: &str,
    z: (&BigFloat, &BigFloat),
    w_op: Option<(&BigFloat, &BigFloat)>,
    which: Which,
    p: u32,
    mode: RoundingMode,
) -> Option<BigFloat> {
    for &extra in &[160u32, 512, 2048] {
        let cb = w.cbracket(fn_id, p + extra, z, w_op);
        let comp = match which {
            Which::Re => cb.re,
            Which::Im => cb.im,
        };
        let Comp::Finite { lo, hi } = comp else {
            return None; // non-finite acb component: the specials lane owns it
        };
        let (lr, _) = lo.round_to_precision(p, mode).expect("p >= 1");
        let (hr, _) = hi.round_to_precision(p, mode).expect("p >= 1");
        if exact_signed_eq(&lr, &hr) {
            return Some(lr);
        }
        // Bracket straddles a p-bit boundary at this oracle precision: tighten.
    }
    None
}

#[derive(Default)]
struct Stats {
    checked: u32,
    skipped: u32,
    residual: u32,
}

/// Assert componentwise correct rounding of `fn_id(z [, w])` under `mode`
/// against the acb oracle, accumulating coverage stats. A finite nonzero
/// component must equal the oracle's certified rounding bit-for-bit; zero /
/// non-finite components and Ziv-cap residuals are counted, not asserted.
fn check_cr(
    w: &mut AcbComplexWorker,
    fn_id: &str,
    z: (&BigFloat, &BigFloat),
    w_op: Option<(&BigFloat, &BigFloat)>,
    p: u32,
    mode: RoundingMode,
    stats: &mut Stats,
) {
    let result = compute(fn_id, z, w_op, mode);
    for (which, pf) in [(Which::Re, &result.re), (Which::Im, &result.im)] {
        if pf.is_nan() || pf.is_infinite() {
            stats.skipped += 1;
            continue;
        }
        match certify(w, fn_id, z, w_op, which, p, mode) {
            Some(cr) => {
                if cr.is_zero() {
                    stats.skipped += 1; // sign-of-zero is the enumerated lane's job
                    continue;
                }
                let part = match which {
                    Which::Re => "re",
                    Which::Im => "im",
                };
                assert!(
                    exact_signed_eq(pf, &cr),
                    "{fn_id}.{part} CR MISMATCH at z=({}, {}) w={:?} p={p} {mode:?}:\n  pfloat = {pf}\n  acb cr = {cr}",
                    z.0,
                    z.1,
                    w_op.map(|(a, b)| format!("({a}, {b})"))
                );
                stats.checked += 1;
            }
            None => stats.residual += 1,
        }
    }
}

/// A deterministic spread of finite complex inputs at precision `p`, both
/// components nonzero (off the negative-real cut) and of modest magnitude (so
/// `cexp` does not overflow): `re, im = n/d` reduced at `p`.
fn input_grid(rng: &mut Rng, p: u32, n: usize) -> Vec<(BigFloat, BigFloat)> {
    let mut g = Vec::with_capacity(n);
    while g.len() < n {
        let rn = rng.int(20);
        let rd = rng.int(8).abs() + 1;
        let inn = rng.int(20);
        let ind = rng.int(8).abs() + 1;
        if inn == 0 {
            continue; // keep im != 0 (off the cut, no zero-output ambiguity)
        }
        let re = bf(rn, p).div(&bf(rd, p), NE).0;
        let im = bf(inn, p).div(&bf(ind, p), NE).0;
        g.push((re, im));
    }
    g
}

fn deep() -> bool {
    std::env::var("PFLOAT_DEEP").is_ok()
}

#[test]
fn dyadic_codec_round_trips_exactly() {
    // Runs without the venv: the wire codec is lossless, so encode-then-decode
    // must equal the input bit-for-bit across a spread of finite values.
    let cases = [
        bf(0, 53),
        bf(1, 53),
        bf(-1, 113),
        bf(1_000_003, 53),
        bf(3, 53).div(&bf(2, 53), NE).0,
        bf(-7, 113).div(&bf(4, 113), NE).0,
        bf(1, 200).div(&bf(3, 200), NE).0, // 1/3, full p-bit mantissa
        bf(5, 24).scale_by_pow2(100).0,
        bf(5, 24).scale_by_pow2(-100).0,
    ];
    for x in cases {
        let rt = encode_decode(&x);
        assert!(
            exact_signed_eq(&x, &rt),
            "codec round-trip changed the value: {x} != {rt}"
        );
    }
}

#[test]
fn acb_csqrt_componentwise_cr() {
    if !acb_lane_available("acb_csqrt_componentwise_cr") {
        return;
    }
    let mut w = AcbComplexWorker::spawn();
    let n = if deep() { 60 } else { 16 };
    let mut stats = Stats::default();
    for &p in &PRECS {
        let mut rng = Rng(0x5172_0001 ^ p as u64);
        for (re, im) in input_grid(&mut rng, p, n) {
            for &mode in &ALL_MODES {
                check_cr(&mut w, "csqrt", (&re, &im), None, p, mode, &mut stats);
            }
        }
    }
    eprintln!(
        "acb csqrt: {} checked, {} skipped, {} ziv-residual",
        stats.checked, stats.skipped, stats.residual
    );
    assert!(
        stats.checked >= 100,
        "csqrt coverage too thin: {}",
        stats.checked
    );
}

#[test]
fn acb_cexp_componentwise_cr() {
    if !acb_lane_available("acb_cexp_componentwise_cr") {
        return;
    }
    let mut w = AcbComplexWorker::spawn();
    let n = if deep() { 60 } else { 16 };
    let mut stats = Stats::default();
    for &p in &PRECS {
        let mut rng = Rng(0xC0FF_EE02 ^ p as u64);
        for (re, im) in input_grid(&mut rng, p, n) {
            for &mode in &ALL_MODES {
                check_cr(&mut w, "cexp", (&re, &im), None, p, mode, &mut stats);
            }
        }
    }
    eprintln!(
        "acb cexp: {} checked, {} skipped, {} ziv-residual",
        stats.checked, stats.skipped, stats.residual
    );
    assert!(
        stats.checked >= 100,
        "cexp coverage too thin: {}",
        stats.checked
    );
}

#[test]
fn acb_clog_componentwise_cr() {
    if !acb_lane_available("acb_clog_componentwise_cr") {
        return;
    }
    let mut w = AcbComplexWorker::spawn();
    let n = if deep() { 60 } else { 16 };
    let mut stats = Stats::default();
    for &p in &PRECS {
        let mut rng = Rng(0x10A6_0003 ^ p as u64);
        for (re, im) in input_grid(&mut rng, p, n) {
            for &mode in &ALL_MODES {
                check_cr(&mut w, "clog", (&re, &im), None, p, mode, &mut stats);
            }
        }
    }
    eprintln!(
        "acb clog: {} checked, {} skipped, {} ziv-residual",
        stats.checked, stats.skipped, stats.residual
    );
    assert!(
        stats.checked >= 100,
        "clog coverage too thin: {}",
        stats.checked
    );
}

#[test]
fn acb_cmul_cdiv_componentwise_cr() {
    if !acb_lane_available("acb_cmul_cdiv_componentwise_cr") {
        return;
    }
    let mut w = AcbComplexWorker::spawn();
    let n = if deep() { 50 } else { 14 };
    let mut stats = Stats::default();
    for &fn_id in &["cmul", "cdiv"] {
        for &p in &PRECS {
            let mut rng = Rng(0xD1F0_0004 ^ (p as u64) ^ fn_id.len() as u64);
            let zs = input_grid(&mut rng, p, n);
            let ws = input_grid(&mut rng, p, n);
            for (z, wv) in zs.iter().zip(ws.iter()) {
                for &mode in &ALL_MODES {
                    check_cr(
                        &mut w,
                        fn_id,
                        (&z.0, &z.1),
                        Some((&wv.0, &wv.1)),
                        p,
                        mode,
                        &mut stats,
                    );
                }
            }
        }
    }
    eprintln!(
        "acb cmul/cdiv: {} checked, {} skipped, {} ziv-residual",
        stats.checked, stats.skipped, stats.residual
    );
    assert!(
        stats.checked >= 150,
        "cmul/cdiv coverage too thin: {}",
        stats.checked
    );
}

#[test]
fn acb_lefevre_muller_seeded_cr() {
    // Seed the components with hard-to-round binary64 inputs for the scalar
    // sub-kernels, so each composed component sits boundary-close and the
    // certified rounding is maximally stressed.
    if !acb_lane_available("acb_lefevre_muller_seeded_cr") {
        return;
    }
    let mut w = AcbComplexWorker::spawn();
    let exp_in = lm_inputs_for("exp").expect("exp corpus");
    let ln_in = lm_inputs_for("ln").expect("ln corpus");
    let mut stats = Stats::default();
    let p = 113;
    let mut used = 0u32;
    for (k, &(xb, _)) in exp_in.iter().enumerate() {
        if !is_finite_nonzero_f64(xb) {
            continue;
        }
        let x = bf_of_f64_bits(xb, p);
        if matches!(
            x.abs().partial_cmp(&bf(8, p)).0,
            Some(core::cmp::Ordering::Greater)
        ) {
            continue; // keep e^x from overflowing the oracle bound
        }
        let (yb, _) = ln_in[k % ln_in.len()];
        if !is_finite_nonzero_f64(yb) {
            continue;
        }
        let mut y = bf_of_f64_bits(yb, p);
        while matches!(
            y.abs().partial_cmp(&bf(2, p)).0,
            Some(core::cmp::Ordering::Greater)
        ) {
            y = y.scale_by_pow2(-1).0;
        }
        for &fn_id in &["csqrt", "cexp", "clog"] {
            check_cr(&mut w, fn_id, (&x, &y), None, p, NE, &mut stats);
        }
        used += 1;
        if used >= 30 {
            break;
        }
    }
    eprintln!(
        "acb L-M seeded: {} checked, {} skipped, {} ziv-residual ({used} seeds)",
        stats.checked, stats.skipped, stats.residual
    );
    assert!(
        stats.checked >= 30,
        "L-M coverage too thin: {}",
        stats.checked
    );
}

#[test]
fn negative_control_one_ulp_off_is_rejected() {
    // The lane has teeth: a value one ULP off the oracle's certified rounding
    // MUST be rejected by the certified-rounding comparison. Build the correct
    // cr, nudge it by next_up, and confirm the comparison fails.
    if !acb_lane_available("negative_control_one_ulp_off_is_rejected") {
        return;
    }
    let mut w = AcbComplexWorker::spawn();
    let mut proved = 0u32;
    for &p in &PRECS {
        let mut rng = Rng(0x0BAD_0005 ^ p as u64);
        for (re, im) in input_grid(&mut rng, p, 12) {
            for &fn_id in &["csqrt", "cexp", "clog"] {
                if let Some(cr) = certify(&mut w, fn_id, (&re, &im), None, Which::Re, p, NE) {
                    if cr.is_zero() {
                        continue;
                    }
                    let bad = cr.next_up().0;
                    assert!(
                        !exact_signed_eq(&bad, &cr),
                        "next_up did not move the value off cr"
                    );
                    proved += 1;
                }
            }
        }
    }
    assert!(
        proved >= 20,
        "negative control proved teeth on only {proved} cases"
    );
}

#[test]
fn ziv_cap_residual_probe() {
    // Probe the documented hard-to-round regimes the kernels name as their
    // failure mode: cexp near y = k*pi/2 (a trig factor near 0) and clog near
    // |z| = 1 (ln(hypot) near 0). Where pfloat's component is finite nonzero AND
    // the oracle certifies a unique rounding, it MUST match; where the oracle
    // itself cannot certify even at +2048 bits, it is logged as a residual, not
    // a failure (pfloat's cap-5 enclosure and the oracle are both bounded).
    if !acb_lane_available("ziv_cap_residual_probe") {
        return;
    }
    let mut w = AcbComplexWorker::spawn();
    let p = 113;
    let mut stats = Stats::default();

    // cexp near y = pi/2 + small: build pi/2 at high precision, perturb.
    let half_pi = bf(1, p + 64).atan2(&bf(0, p + 64), NE).0;
    for k in 0..24i64 {
        let off = bf(if k % 2 == 0 { 1 } else { -1 }, p + 64)
            .scale_by_pow2(-(40 + k))
            .0;
        let y = half_pi.add(&off, NE).0.round_to_precision(p, NE).unwrap().0;
        let x = bf(1, p).div(&bf(4, p), NE).0; // re = 1/4
        for &mode in &ALL_MODES {
            check_cr(&mut w, "cexp", (&x, &y), None, p, mode, &mut stats);
        }
    }

    // clog near |z| = 1: z = (1 + eps) + delta*i with |z| ~ 1.
    for k in 0..24i64 {
        let eps = bf(1, p).scale_by_pow2(-(30 + k)).0;
        let re = bf(1, p).add(&eps, NE).0;
        let im = bf(1, p).scale_by_pow2(-(30 + k)).0;
        for &mode in &ALL_MODES {
            check_cr(&mut w, "clog", (&re, &im), None, p, mode, &mut stats);
        }
    }

    eprintln!(
        "ziv-cap residual probe: {} checked, {} skipped, {} ziv-residual",
        stats.checked, stats.skipped, stats.residual
    );
    // The probe is a measurement; it asserts the checked cases matched (the
    // assertions inside check_cr) and reports the residual count. A nonzero
    // residual is the documented measure-zero hard-to-round caveat, not a bug.
    assert!(
        stats.checked >= 20,
        "residual probe certified too few cases to be meaningful: {}",
        stats.checked
    );
}
