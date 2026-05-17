//! MPFR differential for `BigFloat::ai` (Airy `Ai`).
//!
//! `rug` 1.27 exposes MPFR's `mpfr_ai` (Airy `Ai` of the first kind)
//! for all real arguments — but no `Bi`, `Ai′`, `Bi′`. So `Ai` gets
//! a true MPFR oracle here; `Bi`/`Ai′`/`Bi′` are covered by the
//! checked-in authoritative table, self-consistency, and the
//! Wronskian in `differential_bi.rs`.
//!
//! The oracle is iterated over a bounded representative point table
//! (the `differential_si` shape: a fixed set spanning the Maclaurin
//! regime, plus a few large-`|x|` points at low precision for the
//! DLMF 9.7 asymptotic), not a random sweep — each `Ai` evaluation
//! calls the gamma kernel twice, so an unbounded high-precision
//! sweep is impractically slow.
//!
//! Tolerance is `p − 2` bits, not bit-exact: pfloat does not yet do
//! Ziv correct-rounding (roadmap slice 7c) and the small-`|x|` `Ai`
//! Maclaurin path carries a `c1·f − c2·g` cancellation (the
//! `differential_si` posture).

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{bigfloat_to_rug, mpfr_round_of, rug_from_i64, TRANSCENDENTAL_PRECISIONS};
use pfloat::{BigFloat, RoundingMode};

/// `Ai` rounded at `p` vs MPFR `mpfr_ai`, agreeing to `p − slack`
/// bits. Comparison is done in `rug` space via the bit-exact
/// `BigFloat → rug` converter.
fn assert_ai_close(num: i64, den: i64, p: u32, slack: u32) {
    let x = {
        let n = BigFloat::try_from_i64_exact(num, p).unwrap();
        if den == 1 {
            n
        } else {
            n.div(
                &BigFloat::try_from_i64_exact(den, p).unwrap(),
                RoundingMode::NearestEven,
            )
            .0
        }
    };
    let bf = {
        let (r, _s) = x.ai(RoundingMode::NearestEven);
        bigfloat_to_rug(&r)
    };
    let rg = {
        let n = rug_from_i64(num, p);
        let xr = if den == 1 {
            n
        } else {
            rug::Float::with_val(p, n / rug_from_i64(den, p))
        };
        rug::Float::with_val_round(p, xr.ai_ref(), mpfr_round_of(RoundingMode::NearestEven)).0
    };
    let diff = rug::Float::with_val(p + 16, &bf - &rg).abs();
    let mut tol = rug::Float::with_val(p + 16, rg.clone().abs());
    if tol == 0 {
        tol = rug::Float::with_val(p + 16, 1);
    }
    tol >>= p - slack;
    assert!(
        diff <= tol,
        "Ai({num}/{den}) at p={p}: bf={bf}, rg={rg}, diff={diff}"
    );
}

/// Representative points spanning the convergent Maclaurin regime
/// (both signs; Airy has no parity), checked at every transcendental
/// precision against the MPFR oracle.
const AI_POINTS: &[(i64, i64)] = &[
    (1, 2),
    (1, 1),
    (2, 1),
    (7, 2),
    (5, 1),
    (-1, 1),
    (-2, 1),
    (-5, 1),
    (-7, 2),
];

#[test]
fn ai_matches_mpfr_maclaurin_regime() {
    // p ≤ 256 only: each Ai call invokes the gamma kernel twice, and
    // MPFR's own mpfr_ai is itself costly, so p=1024 over the full
    // table is impractical here. The p=1024 Maclaurin path is pinned
    // independently by the airy-module unit tests (boundary constants
    // and the Wronskian to ~176 bits at p=200).
    for &p in TRANSCENDENTAL_PRECISIONS.iter().filter(|&&p| p <= 256) {
        for &(num, den) in AI_POINTS {
            assert_ai_close(num, den, p, 2);
        }
    }
}

#[test]
fn ai_matches_mpfr_asymptotic_regime() {
    // Just past the p=53 asymptotic threshold (|x| ≳ 128), where
    // pfloat takes the DLMF 9.7.5 / 9.7.9 path while MPFR's mpfr_ai
    // is still tractable. Larger |x| is covered analytically by the
    // airy-module asymptotic unit tests (x=±300, p=150).
    for &n in &[150i64, -150, 180, -180] {
        assert_ai_close(n, 1, 53, 4);
    }
}
