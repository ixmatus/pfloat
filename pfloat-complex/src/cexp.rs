//! Componentwise correctly-rounded complex exponential `cexp` with C99 Annex G
//! §G.6.3.1 special values (ADR-0091).
//!
//! `cexp(x + iy) = e^x (cos y + i sin y)`. `cexp` is entire, so there is no
//! branch cut; the subtleties are the special-value table, which is classified
//! on the INPUT real-part class (a finite huge `x` whose `e^x` overflows is
//! still the finite path, with `OVERFLOW` correct, never merged with the true
//! `x = +∞` rows), and the signed zeros. The load-bearing asymmetry: with
//! `x = +∞` an indeterminate angle gives `(+∞, NaN) + INVALID`, but with
//! `x = −∞` the `e^{−∞} = +0` factor dominates the undefined `cos/sin(∞)`, so
//! the result is `(+0, +0)` with NO `INVALID`.
//!
//! For finite `x` and finite nonzero `y` the two components `e^x cos y` and
//! `e^x sin y` are products of two transcendentals. Because `e^x > 0`, the
//! product's sign equals the trig factor's sign, so the directed enclosure
//! picks the interval endpoints by the trig sign alone (simpler than the
//! divide's denominator sign logic). A trig bracket that straddles zero (`y`
//! near `kπ/2`) has not resolved the component's sign; the guard grows. By
//! Niven, `cos y` and `sin y` are never exactly zero for a finite nonzero
//! dyadic `y`, so the straddle is transient.

use pfloat::{BigFloat, RoundingMode, Sign, Status};

use crate::enclosure::{resolve_bracket, Resolved, GUARDS};

/// `cexp(a + bi)` at output precision `p`, returning `(re, im, status)`. Runs
/// in `BigFloat`; the generic complex `exp` bridges through
/// `RealScalar::to_big`/`from_big`.
pub(crate) fn cexp_big(
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

    // Classify on the real part x first.
    if a.is_nan() {
        // E21-E25: Re = NaN. Im = copysign(0, y) for y = ±0, else NaN. A quiet
        // NaN does not signal.
        let im = if b.is_zero() {
            signed_zero(b, p)
        } else {
            nan(p)
        };
        return (nan(p), im, inv);
    }

    if a.is_infinite() {
        if a.is_sign_positive() {
            // x = +∞.
            if b.is_zero() {
                return (pos_inf(p), signed_zero(b, p), inv); // E9/E10
            }
            if b.is_infinite() {
                return (pos_inf(p), nan(p), Status::INVALID); // E12/E13 [REP +∞]
            }
            if b.is_nan() {
                return (pos_inf(p), nan(p), inv); // E14 [REP +∞]
            }
            // E11: finite nonzero y. (sign(cos y)·∞, sign(sin y)·∞).
            return (
                signed_inf(trig_sign(b, Trig::Cos), p),
                signed_inf(trig_sign(b, Trig::Sin), p),
                inv,
            );
        }
        // x = −∞: e^{−∞} = +0 dominates; no INVALID on the indeterminate-angle
        // rows (the asymmetry with x = +∞).
        if b.is_zero() {
            return (pos_zero(p), signed_zero(b, p), inv); // E15/E16
        }
        if b.is_infinite() || b.is_nan() {
            return (pos_zero(p), pos_zero(p), inv); // E18-E20 [REP +0]
        }
        // E17: finite nonzero y. (sign(cos y)·0, sign(sin y)·0).
        return (
            signed_zero_of(trig_sign(b, Trig::Cos), p),
            signed_zero_of(trig_sign(b, Trig::Sin), p),
            inv,
        );
    }

    // x finite from here.
    if b.is_nan() {
        return (nan(p), nan(p), inv); // E8 (quiet NaN does not signal)
    }
    if b.is_infinite() {
        return (nan(p), nan(p), Status::INVALID); // E6/E7: cos/sin of ∞ undefined
    }
    if b.is_zero() {
        // E1-E4: (e^x, copysign(0, y)). e^0 = 1 is exact via the scalar exp;
        // the imaginary zero is stamped, not formed by e^x·0 (which would be
        // NaN when e^x overflowed).
        let (re, s) = a
            .exp_round(p, mode)
            .expect("p >= 1 (cexp target precision)");
        return (re, signed_zero(b, p), s);
    }

    // E5: x finite, y finite nonzero. The sign-aware product enclosure.
    let mut last: Option<(Resolved, Resolved, Status)> = None;
    for &guard in &GUARDS {
        let w = p.saturating_add(guard);
        let (ex_lo, s_lo) = a
            .exp_round(w, RoundingMode::TowardNegative)
            .expect("w >= 1 (p + guard)");
        let (ex_hi, s_hi) = a
            .exp_round(w, RoundingMode::TowardPositive)
            .expect("w >= 1 (p + guard)");
        // The exp directed pair carries INEXACT (always, e^x is transcendental
        // for x ≠ 0) and OVERFLOW if it saturated; the result is transcendental
        // so it is INEXACT regardless of how the bracket collapses.
        let exp_status = s_lo | s_hi | Status::INEXACT;

        let cos_lo = cos_round(b, w, RoundingMode::TowardNegative);
        let cos_hi = cos_round(b, w, RoundingMode::TowardPositive);
        let sin_lo = sin_round(b, w, RoundingMode::TowardNegative);
        let sin_hi = sin_round(b, w, RoundingMode::TowardPositive);

        let (re_lo, re_hi) = product_bracket(&ex_lo, &ex_hi, &cos_lo, &cos_hi);
        let (im_lo, im_hi) = product_bracket(&ex_lo, &ex_hi, &sin_lo, &sin_hi);
        let re = resolve_bracket(&re_lo, &re_hi, p, mode);
        let im = resolve_bracket(&im_lo, &im_hi, p, mode);
        if re.converged && im.converged {
            return (re.value, im.value, exp_status);
        }
        last = Some((re, im, exp_status));
    }
    let (re, im, exp_status) = last.expect("GUARDS is non-empty");
    (re.value, im.value, exp_status)
}

/// The directed bracket of `e^x · t` where `e^x ∈ [ex_lo, ex_hi]` (both > 0)
/// and `t ∈ [t_lo, t_hi]`. Because `e^x > 0` the product sign follows `t`, so
/// the endpoints are chosen by the sign of the `t` bracket. A `t` bracket that
/// straddles zero yields a zero-spanning product bracket, which
/// [`resolve_bracket`] reports as not converged so the caller grows the guard.
fn product_bracket(
    ex_lo: &BigFloat,
    ex_hi: &BigFloat,
    t_lo: &BigFloat,
    t_hi: &BigFloat,
) -> (BigFloat, BigFloat) {
    use RoundingMode::{TowardNegative as TN, TowardPositive as TP};
    let t_lo_neg = t_lo.is_sign_negative() && !t_lo.is_zero();
    let t_hi_pos = t_hi.is_sign_positive() && !t_hi.is_zero();
    if !t_lo_neg {
        // t ≥ 0: product ∈ [ex_lo·t_lo, ex_hi·t_hi].
        (ex_lo.mul(t_lo, TN).0, ex_hi.mul(t_hi, TP).0)
    } else if !t_hi_pos {
        // t ≤ 0: product ∈ [ex_hi·t_lo, ex_lo·t_hi].
        (ex_hi.mul(t_lo, TN).0, ex_lo.mul(t_hi, TP).0)
    } else {
        // t straddles 0: widest span uses the largest e^x for both ends.
        (ex_hi.mul(t_lo, TN).0, ex_hi.mul(t_hi, TP).0)
    }
}

#[derive(Clone, Copy)]
enum Trig {
    Cos,
    Sin,
}

/// The sign of `cos y` or `sin y` for a finite nonzero `y`, resolved by
/// evaluating at growing precision until the value is nonzero. By Niven both
/// are nonzero for a finite nonzero dyadic `y`, so this terminates; the
/// fallback sign is reached only on a measure-zero unresolved input.
fn trig_sign(y: &BigFloat, which: Trig) -> Sign {
    for &guard in &GUARDS {
        let w = 64u32.saturating_add(guard);
        let v = match which {
            Trig::Cos => cos_round(y, w, RoundingMode::NearestEven),
            Trig::Sin => sin_round(y, w, RoundingMode::NearestEven),
        };
        if !v.is_zero() && !v.is_nan() {
            return v.sign();
        }
    }
    Sign::Positive
}

fn cos_round(y: &BigFloat, w: u32, mode: RoundingMode) -> BigFloat {
    y.cos_round(w, mode).expect("w >= 1").0
}

fn sin_round(y: &BigFloat, w: u32, mode: RoundingMode) -> BigFloat {
    y.sin_round(w, mode).expect("w >= 1").0
}

fn pos_inf(p: u32) -> BigFloat {
    BigFloat::try_new_infinity(Sign::Positive, p).expect("p >= 1")
}

fn pos_zero(p: u32) -> BigFloat {
    BigFloat::try_new_zero(Sign::Positive, p).expect("p >= 1")
}

fn nan(p: u32) -> BigFloat {
    BigFloat::try_new_quiet_nan(Sign::Positive, p, &[]).expect("p >= 1")
}

fn signed_inf(sign: Sign, p: u32) -> BigFloat {
    BigFloat::try_new_infinity(sign, p).expect("p >= 1")
}

fn signed_zero_of(sign: Sign, p: u32) -> BigFloat {
    BigFloat::try_new_zero(sign, p).expect("p >= 1")
}

/// `copysign(+0, src)`: a signed zero carrying `src`'s sign.
fn signed_zero(src: &BigFloat, p: u32) -> BigFloat {
    BigFloat::try_new_zero(Sign::Positive, p)
        .expect("p >= 1")
        .copysign(src)
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

    // |a + bi| within 2^-(p-20) of the integer n (a coarse closeness check).
    fn near_int(x: &BigFloat, n: i64, p: u32) -> bool {
        let d = x.sub(&bfp(n, p), NE).0.abs();
        matches!(
            d.partial_cmp(&bfp(1, p).scale_by_pow2(-(p as i64) + 20).0)
                .0,
            Some(Ordering::Less)
        )
    }

    #[test]
    fn cexp_zero_is_one_exact() {
        // cexp(0 + 0i) = 1 + 0i, exact (e^0 = 1, sin(0) = 0).
        let (re, im, s) = cexp_big(&pz(), &pz(), 64, NE);
        assert!(matches!(re.partial_cmp(&bf(1)).0, Some(Ordering::Equal)));
        assert!(im.is_zero() && im.is_sign_positive());
        assert!(!s.inexact(), "cexp(0) is exact");
    }

    #[test]
    fn cexp_real_axis_matches_scalar_exp() {
        // cexp(x + 0i) = (e^x, +0); the real part is the scalar exp bit-for-bit.
        let (re, im, s) = cexp_big(&bf(3), &pz(), 113, NE);
        let scalar = bfp(3, 113).exp_round(113, NE).unwrap().0;
        assert_eq!(re.partial_cmp(&scalar).0, Some(Ordering::Equal));
        assert!(im.is_zero() && im.is_sign_positive());
        assert!(s.inexact(), "e^3 is irrational");
    }

    #[test]
    fn cexp_real_axis_signed_zero_imaginary() {
        // cexp(2 − 0i) = (e^2, −0): the imaginary zero follows y's sign.
        let (_, im, _) = cexp_big(&bf(2), &nz(), 64, NE);
        assert!(im.is_zero() && im.is_sign_negative());
    }

    #[test]
    fn cexp_modulus_is_exp_real_part() {
        // |cexp(1 + 2i)| = e^1. Check hypot(re, im) ≈ e.
        let (re, im, _) = cexp_big(&bf(1), &bf(2), 200, NE);
        let modulus = re.hypot(&im, NE).0;
        let e = bfp(1, 200).exp_round(200, NE).unwrap().0;
        let d = modulus.sub(&e, NE).0.abs();
        assert!(
            matches!(
                d.partial_cmp(&bfp(1, 200).scale_by_pow2(-170).0).0,
                Some(Ordering::Less)
            ),
            "|cexp(1+2i)| ≈ e, got {modulus}"
        );
    }

    #[test]
    fn cexp_additive_law() {
        // cexp(z + w) = cexp(z)·cexp(w). z = 1 + 1i, w = 1 + 1i, z + w = 2 + 2i.
        let p = 200;
        let (lr, li, _) = cexp_big(&bfp(2, p), &bfp(2, p), p, NE);
        let (zr, zi, _) = cexp_big(&bfp(1, p), &bfp(1, p), p, NE);
        // (zr + zi·i)² = (zr² − zi²) + (2·zr·zi)i.
        let rr = zr.mul_sub_mul(&zr, &zi, &zi, NE).0;
        let ri = zr.mul(&zi, NE).0.scale_by_pow2(1).0;
        let dr = lr.sub(&rr, NE).0.abs();
        let di = li.sub(&ri, NE).0.abs();
        let tol = bfp(1, p).scale_by_pow2(-170).0;
        assert!(
            matches!(dr.partial_cmp(&tol).0, Some(Ordering::Less)),
            "re mismatch {dr}"
        );
        assert!(
            matches!(di.partial_cmp(&tol).0, Some(Ordering::Less)),
            "im mismatch {di}"
        );
    }

    #[test]
    fn cexp_quarter_turn() {
        // cexp(i·π/2) = i: re ≈ 0, im ≈ 1. Build π/2 at high precision.
        let p = 256;
        let half_pi = {
            // π/2 = 2·atan(1). atan2(1, 0) = π/2.
            bfp(1, p).atan2(&bfp(0, p), NE).0
        };
        let (re, im, _) = cexp_big(&bfp(0, p), &half_pi, p, NE);
        // |re| is tiny, im ≈ 1.
        assert!(
            matches!(
                re.abs().partial_cmp(&bfp(1, p).scale_by_pow2(-200).0).0,
                Some(Ordering::Less)
            ),
            "re(cexp(iπ/2)) ≈ 0, got {re}"
        );
        assert!(near_int(&im, 1, p), "im(cexp(iπ/2)) ≈ 1, got {im}");
    }

    #[test]
    fn cexp_finite_plus_infinite_imaginary_is_invalid() {
        // E6/E7: cexp(1 + ∞i) = NaN + NaN·i + INVALID.
        let (re, im, s) = cexp_big(&bf(1), &pinf(), 64, NE);
        assert!(re.is_nan() && im.is_nan());
        assert!(s.invalid());
    }

    #[test]
    fn cexp_neg_infinity_dominates_without_invalid() {
        // E18: cexp(−∞ + ∞i) = +0 + 0i, NO INVALID (the asymmetry with +∞).
        let (re, im, s) = cexp_big(&ninf(), &pinf(), 64, NE);
        assert!(re.is_zero() && re.is_sign_positive());
        assert!(im.is_zero());
        assert!(!s.invalid(), "−∞ real part suppresses INVALID");
    }

    #[test]
    fn cexp_pos_infinity_plus_infinite_imaginary_is_invalid() {
        // E12: cexp(+∞ + ∞i) = +∞ + NaN·i + INVALID.
        let (re, im, s) = cexp_big(&pinf(), &pinf(), 64, NE);
        assert!(re.is_infinite() && re.is_sign_positive());
        assert!(im.is_nan());
        assert!(s.invalid());
    }

    #[test]
    fn cexp_nan_minus_zero_keeps_zero_sign() {
        // E22: cexp(NaN − 0i) = NaN − 0i (the imaginary-zero sign is pinned by
        // the conjugation symmetry cexp(conj z) = conj(cexp z)).
        let (re, im, s) = cexp_big(&qnan(), &nz(), 64, NE);
        assert!(re.is_nan());
        assert!(im.is_zero() && im.is_sign_negative());
        assert!(!s.invalid());
    }

    #[test]
    fn cexp_pos_infinity_real_axis() {
        // E9: cexp(+∞ + 0i) = +∞ + 0i.
        let (re, im, _) = cexp_big(&pinf(), &pz(), 64, NE);
        assert!(re.is_infinite() && re.is_sign_positive());
        assert!(im.is_zero() && im.is_sign_positive());
    }

    #[test]
    fn cexp_neg_infinity_real_axis() {
        // E15: cexp(−∞ + 0i) = +0 + 0i.
        let (re, im, _) = cexp_big(&ninf(), &pz(), 64, NE);
        assert!(re.is_zero() && re.is_sign_positive());
        assert!(im.is_zero() && im.is_sign_positive());
    }

    #[test]
    fn cexp_pos_infinity_finite_y_signs_from_trig() {
        // E11: cexp(+∞ + 1i) = (sign(cos 1)·∞, sign(sin 1)·∞). cos 1 > 0,
        // sin 1 > 0, so both parts are +∞.
        let (re, im, _) = cexp_big(&pinf(), &bf(1), 64, NE);
        assert!(re.is_infinite() && re.is_sign_positive());
        assert!(im.is_infinite() && im.is_sign_positive());
        // cexp(+∞ + 3i): cos 3 < 0, sin 3 > 0 → (−∞, +∞).
        let (re3, im3, _) = cexp_big(&pinf(), &bf(3), 64, NE);
        assert!(re3.is_infinite() && re3.is_sign_negative());
        assert!(im3.is_infinite() && im3.is_sign_positive());
    }

    #[test]
    fn cexp_directed_modes_bracket() {
        // cexp(1 + 1i) real part e·cos(1) is irrational; directed rounding
        // gives TN ≤ NE ≤ TP, TN < TP.
        let (re_ne, _, _) = cexp_big(&bf(1), &bf(1), 53, NE);
        let (re_tn, _, _) = cexp_big(&bf(1), &bf(1), 53, RoundingMode::TowardNegative);
        let (re_tp, _, _) = cexp_big(&bf(1), &bf(1), 53, RoundingMode::TowardPositive);
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
