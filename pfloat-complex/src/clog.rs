//! Componentwise correctly-rounded principal complex logarithm `clog` with
//! C99 Annex G §G.6.3.2 branch cuts (ADR-0091).
//!
//! `clog(z) = ln|z| + i·arg(z) = ln(hypot(x, y)) + i·atan2(y, x)`. The
//! imaginary part is exactly `atan2(y, x)`, so it delegates entirely to the
//! scalar `atan2` kernel, which already carries the Annex G branch cut on the
//! negative real axis and the signed-zero discrimination (`arg(−1 + 0i) = +π`,
//! `arg(−1 − 0i) = −π`) in all five rounding modes. The four poles fall out of
//! this composition: `clog(±0 ± 0i) = −∞ + i·atan2(±0, ±0)`, with the `−∞` and
//! its `DIV_BY_ZERO` from `ln(+0)`.
//!
//! The real part `ln(hypot(x, y))` is enclosed by a directed pair: `hypot`
//! bracketed, then `ln` (monotone increasing) applied to each end. The named
//! failure mode is catastrophic cancellation when `|z|` is near 1, where
//! `ln(hypot) ≈ 0` loses leading bits; the enclosure observes the bracket
//! straddling 0 and grows the working precision until it resolves, bounded by
//! the same measure-zero cap caveat the divide carries. `clog(1 + 0i) = +0`
//! falls out exactly, because `hypot(1, 0) = 1` and the scalar `ln(1) = +0`.

use pfloat::{BigFloat, RoundingMode, Sign, Status};

use crate::enclosure::{resolve_bracket, GUARDS};

/// `clog(a + bi)` at output precision `p`, returning `(re, im, status)`. Runs
/// in `BigFloat`; the generic complex `log` bridges through
/// `RealScalar::to_big`/`from_big`.
pub(crate) fn clog_big(
    a: &BigFloat,
    b: &BigFloat,
    p: u32,
    mode: RoundingMode,
) -> (BigFloat, BigFloat, Status) {
    let inv = if a.is_signaling_nan() || b.is_signaling_nan() {
        Status::INVALID
    } else {
        Status::OK
    };

    // The imaginary part is arg(z) = atan2(y, x): one correctly-rounded scalar
    // kernel carrying the whole signed-zero/∞ branch table.
    let (im, s_im) = b.atan2_round(a, p, mode).expect("p >= 1");

    // The real part ln(hypot(x, y)). The hypot class is fixed by the inputs, so
    // its special values are dispatched directly; only a finite nonzero |z|
    // needs the enclosure.
    let (re, s_re) = if a.is_infinite() || b.is_infinite() {
        // hypot with an infinite part is +∞; ln(+∞) = +∞ (∞ dominates a NaN
        // part, so clog(NaN + ∞i) has real part +∞, not NaN).
        (pos_inf(p), Status::OK)
    } else if a.is_nan() || b.is_nan() {
        // No infinity: hypot is NaN, ln(NaN) = NaN.
        (nan(p), Status::OK)
    } else if a.is_zero() && b.is_zero() {
        // hypot(±0, ±0) = +0; ln(+0) = −∞ + DIV_BY_ZERO (the poles).
        (neg_inf(p), Status::DIV_BY_ZERO)
    } else {
        ln_hypot_enclosure(a, b, p, mode)
    };

    (re, im, s_re | s_im | inv)
}

/// The directed-pair enclosure of `ln(hypot(x, y))` for finite `x, y` not both
/// zero. `hypot` is bracketed `[h_lo, h_hi]`, then `ln` (increasing) applied to
/// each end gives `[ln(h_lo), ln(h_hi)] ⊇ ln(hypot)`. Near `|z| = 1` the
/// bracket straddles `ln(1) = 0`; the loop grows the guard until it resolves.
fn ln_hypot_enclosure(
    a: &BigFloat,
    b: &BigFloat,
    p: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    use RoundingMode::{TowardNegative as TN, TowardPositive as TP};
    let mut last = None;
    for &guard in &GUARDS {
        let w = p.saturating_add(guard);
        let h_lo = a.hypot_round(b, w, TN).expect("w >= 1").0;
        let h_hi = a.hypot_round(b, w, TP).expect("w >= 1").0;
        let re_lo = h_lo.ln_round(w, TN).expect("w >= 1").0;
        let re_hi = h_hi.ln_round(w, TP).expect("w >= 1").0;
        let r = resolve_bracket(&re_lo, &re_hi, p, mode);
        if r.converged {
            return (r.value, r.status);
        }
        last = Some(r);
    }
    let r = last.expect("GUARDS is non-empty");
    (r.value, r.status)
}

fn pos_inf(p: u32) -> BigFloat {
    BigFloat::try_new_infinity(Sign::Positive, p).expect("p >= 1")
}

fn neg_inf(p: u32) -> BigFloat {
    BigFloat::try_new_infinity(Sign::Negative, p).expect("p >= 1")
}

fn nan(p: u32) -> BigFloat {
    BigFloat::try_new_quiet_nan(Sign::Positive, p, &[]).expect("p >= 1")
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    fn bf(n: i64) -> BigFloat {
        BigFloat::try_from_i64_exact(n, 64).unwrap()
    }
    fn bfp(n: i64, p: u32) -> BigFloat {
        BigFloat::try_from_i64_exact(n, p).unwrap()
    }
    fn pz() -> BigFloat {
        BigFloat::try_new_zero(Sign::Positive, 64).unwrap()
    }
    fn nz() -> BigFloat {
        BigFloat::try_new_zero(Sign::Negative, 64).unwrap()
    }
    fn pinf() -> BigFloat {
        BigFloat::try_new_infinity(Sign::Positive, 64).unwrap()
    }
    fn ninf() -> BigFloat {
        BigFloat::try_new_infinity(Sign::Negative, 64).unwrap()
    }
    fn qnan() -> BigFloat {
        BigFloat::try_new_quiet_nan(Sign::Positive, 64, &[]).unwrap()
    }
    const NE: RoundingMode = RoundingMode::NearestEven;

    fn near(x: &BigFloat, target: &BigFloat, p: u32) -> bool {
        let d = x.sub(target, NE).0.abs();
        matches!(
            d.partial_cmp(&bfp(1, p).scale_by_pow2(-(p as i64) + 20).0)
                .0,
            Some(Ordering::Less)
        )
    }

    #[test]
    fn clog_one_is_exact_zero() {
        // clog(1 + 0i) = +0 + 0i, exact (hypot = 1, ln(1) = +0; atan2(0,1) = +0).
        let (re, im, s) = clog_big(&bf(1), &pz(), 64, NE);
        assert!(re.is_zero() && re.is_sign_positive());
        assert!(im.is_zero() && im.is_sign_positive());
        assert!(!s.inexact(), "clog(1) is exact");
    }

    #[test]
    fn clog_real_axis_matches_scalar_ln() {
        // clog(e + 0i) real part = ln(hypot(e, 0)) = ln(e) = 1, imaginary +0.
        let p = 200;
        let e = bfp(1, p).exp_round(p, NE).0;
        let (re, im, _) = clog_big(&e, &bfp(0, p), p, NE);
        assert!(near(&re, &bfp(1, p), p), "ln(e) ≈ 1, got {re}");
        assert!(im.is_zero() && im.is_sign_positive());
    }

    #[test]
    fn clog_modulus_and_phase() {
        // clog(0 + 2i): |z| = 2, arg = +π/2. re = ln 2, im = π/2.
        let p = 200;
        let (re, im, _) = clog_big(&bfp(0, p), &bfp(2, p), p, NE);
        let ln2 = bfp(2, p).ln(NE).0;
        let half_pi = bfp(1, p).atan2(&bfp(0, p), NE).0;
        assert!(near(&re, &ln2, p), "re = ln 2");
        assert!(near(&im, &half_pi, p), "im = π/2");
    }

    #[test]
    fn clog_branch_cut_signed_zero() {
        // The negative real axis cut: clog(−2 + 0i) = ln 2 + iπ,
        // clog(−2 − 0i) = ln 2 − iπ. Real parts equal, imaginary parts opposite.
        let p = 200;
        let ln2 = bfp(2, p).ln(NE).0;
        let pi = bfp(1, p).atan2(&bfp(0, p), NE).0.scale_by_pow2(1).0; // 2·(π/2) = π
        let (re_u, im_u, _) = clog_big(&bfp(-2, p), &bfp(0, p).copysign(&pz()), p, NE);
        let (re_l, im_l, _) = clog_big(&bfp(-2, p), &bfp(0, p).copysign(&nz()), p, NE);
        assert!(
            near(&re_u, &ln2, p) && near(&re_l, &ln2, p),
            "both re = ln 2"
        );
        assert!(
            im_u.is_sign_positive() && near(&im_u, &pi, p),
            "upper cut → +π"
        );
        assert!(im_l.is_sign_negative(), "lower cut → −π");
    }

    #[test]
    fn clog_exp_inverse() {
        // clog(cexp-like input): clog(z) for z = 3 + 4i should give
        // ln 5 + i·atan2(4,3) (|z| = 5). Check re = ln 5.
        let p = 200;
        let (re, im, _) = clog_big(&bfp(3, p), &bfp(4, p), p, NE);
        let ln5 = bfp(5, p).ln(NE).0;
        let phase = bfp(4, p).atan2(&bfp(3, p), NE).0;
        assert!(near(&re, &ln5, p), "re = ln 5");
        assert!(near(&im, &phase, p), "im = atan2(4, 3)");
    }

    #[test]
    fn clog_poles() {
        // The four poles: clog(±0 ± 0i) = −∞ + i·atan2(±0, ±0), DIV_BY_ZERO.
        // x=+0,y=+0 → −∞ + i0; x=−0,y=+0 → −∞ + iπ; x=+0,y=−0 → −∞ − i0;
        // x=−0,y=−0 → −∞ − iπ.
        let (re, im, s) = clog_big(&pz(), &pz(), 64, NE);
        assert!(re.is_infinite() && re.is_sign_negative());
        assert!(im.is_zero() && im.is_sign_positive());
        assert!(s.div_by_zero());
        let (_, im2, _) = clog_big(&nz(), &pz(), 64, NE);
        assert!(
            im2.is_sign_positive() && !im2.is_zero(),
            "clog(−0+0i) im = +π"
        );
        let (_, im3, _) = clog_big(&pz(), &nz(), 64, NE);
        assert!(
            im3.is_zero() && im3.is_sign_negative(),
            "clog(+0−0i) im = −0"
        );
        let (_, im4, _) = clog_big(&nz(), &nz(), 64, NE);
        assert!(
            im4.is_sign_negative() && !im4.is_zero(),
            "clog(−0−0i) im = −π"
        );
    }

    #[test]
    fn clog_infinite_real_part_dominates_nan() {
        // clog(NaN + ∞i): hypot(NaN, ∞) = +∞, so re = +∞ (NOT NaN); im = NaN.
        let (re, im, _) = clog_big(&qnan(), &pinf(), 64, NE);
        assert!(re.is_infinite() && re.is_sign_positive());
        assert!(im.is_nan());
        // clog(−∞ + 3i): re = +∞, im = atan2(3, −∞) = +π.
        let (re2, im2, _) = clog_big(&ninf(), &bf(3), 64, NE);
        assert!(re2.is_infinite() && re2.is_sign_positive());
        assert!(im2.is_sign_positive() && !im2.is_nan());
    }

    #[test]
    fn clog_nan_without_infinity_is_nan() {
        let (re, im, s) = clog_big(&qnan(), &bf(2), 64, NE);
        assert!(re.is_nan() && im.is_nan());
        assert!(!s.invalid(), "quiet NaN does not signal");
    }

    #[test]
    fn clog_signaling_nan_raises_invalid() {
        let snan = BigFloat::try_new_signaling_nan(Sign::Positive, 64, &[]).unwrap();
        let (_, _, s) = clog_big(&snan, &bf(2), 64, NE);
        assert!(s.invalid());
    }

    #[test]
    fn clog_below_unit_circle_is_negative_real_part() {
        // |z| < 1 → ln(hypot) < 0. clog(0.5 + 0i) re = ln(0.5) = −ln 2 < 0.
        let p = 200;
        let half = bfp(1, p).scale_by_pow2(-1).0; // 0.5
        let (re, _, _) = clog_big(&half, &bfp(0, p), p, NE);
        assert!(re.is_sign_negative(), "ln(0.5) < 0");
        let neg_ln2 = bfp(2, p).ln(NE).0.negated();
        assert!(near(&re, &neg_ln2, p), "re = −ln 2");
    }

    #[test]
    fn clog_directed_modes_bracket() {
        // clog(3 + 4i) real part ln 5 is irrational: TN ≤ NE ≤ TP, TN < TP.
        let (re_ne, _, _) = clog_big(&bf(3), &bf(4), 53, NE);
        let (re_tn, _, _) = clog_big(&bf(3), &bf(4), 53, RoundingMode::TowardNegative);
        let (re_tp, _, _) = clog_big(&bf(3), &bf(4), 53, RoundingMode::TowardPositive);
        assert_eq!(re_tn.partial_cmp(&re_tp).0, Some(Ordering::Less));
        assert!(matches!(
            re_ne.partial_cmp(&re_tn).0,
            Some(Ordering::Greater | Ordering::Equal)
        ));
        assert!(matches!(
            re_ne.partial_cmp(&re_tp).0,
            Some(Ordering::Less | Ordering::Equal)
        ));
    }
}
