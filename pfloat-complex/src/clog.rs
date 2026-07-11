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

use core::cmp::Ordering;

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
    // Near |z| = 1 the bracket straddles ln(1) = 0 to a depth set by
    // the INPUT structure, not the output precision: for z = 1 + εi
    // the real part is ~ε²/2, so a static guard schedule topping out
    // at p + 1024 cannot resolve ε deeper than ~512 bits and the
    // exhausted fall-through silently returned a collapsed +0 with
    // INEXACT (pf-qm8a, ADR-0100: clog(1 + 2^-545 i) at p64; the
    // constructible positive-measure band past ADR-0091's
    // measure-zero caveat). The depth is bounded by the components'
    // own grids: a² + b² − 1, if nonzero, is a dyadic no smaller
    // than 2^bot with bot = min(2(e_a − p_a + 1), 2(e_b − p_b + 1)),
    // so growing the guard to −bot + 64 always either resolves the
    // bracket or hits the exactly-Pythagorean case (a² + b² = 1),
    // where the directed hypot pair collapses onto exactly 1 and
    // ln(1) = 0 converges exactly. The schedule is unchanged for
    // everything the old GUARDS already handled.
    let sq_bot = |v: &BigFloat| match v.parts() {
        pfloat::Parts::Normal { exponent, .. } => exponent
            .saturating_sub(i64::from(v.precision()))
            .saturating_add(1)
            .saturating_mul(2),
        _ => 0i64,
    };
    let bot = sq_bot(a).min(sq_bot(b));
    // bot > 0 (both components' grids sit above 2^0) means |z| is far
    // from 1 and the first guard resolves; clamp instead of feeding a
    // negative through the conversion (the unclamped form silently
    // produced a u32::MAX cap — caught by the slice's adversarial
    // verification).
    let guard_cap = u32::try_from(bot.saturating_neg().saturating_add(64).max(0))
        .unwrap_or(u32::MAX)
        .max(GUARDS[GUARDS.len() - 1]);

    let mut last = None;
    let mut guard = GUARDS[0];
    loop {
        let w = p.saturating_add(guard);
        let h_lo = a.hypot_round(b, w, TN).expect("w >= 1").0;
        let h_hi = a.hypot_round(b, w, TP).expect("w >= 1").0;

        // A bracket whose BOTH ends are exactly 1 with both components
        // nonzero is always a lie: no nontrivial dyadic point lies on
        // the unit circle (x² + y² = 2^2k has only the trivial
        // Gaussian-integer solutions), so the truth is strictly off 1
        // and the scalar hypot has merely exhausted its own Ziv cap,
        // returning a falsely-exact 1 (its target was too small to see
        // the 2^-2d offset). Accepting it would let ln(1) = 0 converge
        // "exactly" on the first iteration and bypass the depth-scaled
        // growth below entirely — the slice's adversarial verification
        // reproduced exactly that at depth ≥ 576 (component status OK,
        // worse than the original defect). Treat it as unresolved and
        // keep growing; by depth ≤ −bot/2 the widened target lets
        // hypot resolve genuinely.
        let one = BigFloat::try_from_i64_exact(1, 1).expect("precision >= 1");
        let both_ends_one = matches!(h_lo.partial_cmp(&one).0, Some(Ordering::Equal))
            && matches!(h_hi.partial_cmp(&one).0, Some(Ordering::Equal));
        let collapsed_lie = both_ends_one && !a.is_zero() && !b.is_zero();

        if !collapsed_lie {
            let re_lo = h_lo.ln_round(w, TN).expect("w >= 1").0;
            let re_hi = h_hi.ln_round(w, TP).expect("w >= 1").0;
            let r = resolve_bracket(&re_lo, &re_hi, p, mode);
            if r.converged {
                return (r.value, r.status);
            }
            last = Some(r);
        }
        if guard >= guard_cap {
            break;
        }
        guard = guard.saturating_mul(2).min(guard_cap);
    }
    // Past the structure-derived cap: unreachable for finite-precision
    // inputs by the bound above; the residual is the documented
    // measure-zero caveat (ADR-0091 posture), reported INEXACT —
    // never OK. `last` is None only if every iteration was the
    // collapsed-lie bracket (equally unreachable: the cap exceeds the
    // depth the lie requires); return an explicit INEXACT zero rather
    // than panic or claim exactness.
    match last {
        Some(r) => (r.value, r.status),
        None => (
            BigFloat::try_new_zero(Sign::Positive, p).expect("p >= 1"),
            Status::INEXACT,
        ),
    }
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
        let e = bfp(1, p).exp_round(p, NE).unwrap().0;
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
