//! Componentwise correctly-rounded principal complex square root `csqrt`
//! with C99 Annex G §G.6.4.2 branch cuts (ADR-0091).
//!
//! The principal branch has `Re(csqrt) ≥ 0` and a branch cut along the
//! negative real axis, continuous from above: `csqrt(−4 + 0i) = +2i` while
//! `csqrt(−4 − 0i) = −2i`. The cut and the signed-zero discrimination are a
//! semantic convention layered above rounding, fixed by Annex G; the named
//! failure mode is a wrong-branch result when a caller feeds an unsigned zero
//! where the sign of zero was the only distinguishing information.
//!
//! The interior formula is Kahan's cancellation-robust reformulation. The
//! naive `u = sqrt((|z| + x)/2)`, `v = sign(y)·sqrt((|z| − x)/2)` loses every
//! leading bit of `|z| − x` when `x > 0` and `|y|` is small. Computing instead
//! `w = sqrt((|x| + hypot(x, y))/2)`, which only ever ADDS, and deriving the
//! smaller component by division has no cancellation:
//!
//! - `x ≥ 0`: `u = w`, `v = y / (2w)`.
//! - `x < 0`: `v = copysign(w, y)`, `u = |y| / (2w)`.
//!
//! `w` is irrational, so each component is enclosed by a directed pair at a
//! growing working precision (the [`crate::enclosure`] machinery) and rounded
//! once. The real-axis cases (`y = ±0`) are dispatched exactly ahead of the
//! loop: their zero component is stamped with `copysign(0, y)` directly, never
//! routed through a division (whose `±0` sign would come from the rounding
//! mode, not the input).

use pfloat::{BigFloat, RoundingMode, Sign, Status};

use crate::enclosure::{resolve_bracket, Resolved, GUARDS};

/// Principal `csqrt(a + bi)` at output precision `p`, returning
/// `(re, im, status)`. Runs in `BigFloat`; the generic complex `sqrt` bridges
/// through `RealScalar::to_big`/`from_big`, because the enclosure's working
/// precision exceeds any `FixedFloat<PREC>`.
pub(crate) fn csqrt_big(
    a: &BigFloat,
    b: &BigFloat,
    p: u32,
    mode: RoundingMode,
) -> (BigFloat, BigFloat, Status) {
    // A signaling NaN raises INVALID even where an infinity overrides the
    // value (the §9.2.1 hypot rule, applied componentwise).
    let inv = if a.is_signaling_nan() || b.is_signaling_nan() {
        Status::INVALID
    } else {
        Status::OK
    };

    // Q1/Q2: an infinite imaginary part dominates everything, even x = ±∞ or
    // x = NaN. csqrt(any ± ∞i) = +∞ ± ∞i.
    if b.is_infinite() {
        return (pos_inf(p), inf_signed(b, p), inv);
    }
    // Q3/Q4: x = +∞. Re = +∞; Im = NaN if y is NaN, else copysign(+0, y).
    if a.is_infinite() && a.is_sign_positive() {
        let im = if b.is_nan() {
            nan(p)
        } else {
            signed_zero(b, p)
        };
        return (pos_inf(p), im, inv);
    }
    // Q5/Q6: x = −∞. Im = copysign(+∞, y) (the sign of a NaN y is the [REP]
    // free choice); Re = NaN if y is NaN, else +0.
    if a.is_infinite() && a.is_sign_negative() {
        let re = if b.is_nan() {
            nan(p)
        } else {
            BigFloat::try_new_zero(Sign::Positive, p).expect("p >= 1")
        };
        return (re, inf_signed(b, p), inv);
    }
    // Q7/Q8/Q9: a NaN operand with no infinity. (NaN, NaN); INVALID only for a
    // signaling NaN (a quiet NaN propagates without signaling).
    if a.is_nan() || b.is_nan() {
        return (nan(p), nan(p), inv);
    }

    // Both parts finite from here.

    // Real axis y = ±0: dispatch exactly, stamping the zero component's sign
    // with copysign(0, y) rather than deriving it from a division.
    if b.is_zero() {
        if a.is_zero() {
            // Q10/Q11: origin. (+0, copysign(0, y)).
            return (
                BigFloat::try_new_zero(Sign::Positive, p).expect("p >= 1"),
                signed_zero(b, p),
                Status::OK,
            );
        }
        if a.is_sign_positive() {
            // Q14/Q15: x > 0. (sqrt(x), copysign(0, y)).
            let (re, s) = a.sqrt_round(p, mode).expect("p >= 1");
            return (re, signed_zero(b, p), s);
        }
        // Q12/Q13: x < 0. (+0, copysign(sqrt(|x|), y)).
        let (root, s) = a.abs().sqrt_round(p, mode).expect("p >= 1");
        return (
            BigFloat::try_new_zero(Sign::Positive, p).expect("p >= 1"),
            root.copysign(b),
            s,
        );
    }

    // Q16: x finite, y finite nonzero. The Kahan enclosure.
    let mut last: Option<(Resolved, Resolved)> = None;
    for &guard in &GUARDS {
        let w = p.saturating_add(guard);
        let (re, im) = kahan_brackets(a, b, w, p, mode);
        if re.converged && im.converged {
            return (re.value, im.value, re.status | im.status);
        }
        last = Some((re, im));
    }
    let (re, im) = last.expect("GUARDS is non-empty");
    (re.value, im.value, re.status | im.status)
}

/// One enclosure iteration: bracket `w = sqrt((|x| + hypot(x,y))/2)` with a
/// directed pair at working precision `w`, derive the `(u, v)` component
/// brackets by Kahan's branch, and resolve each to precision `p`.
fn kahan_brackets(
    a: &BigFloat,
    b: &BigFloat,
    w: u32,
    p: u32,
    mode: RoundingMode,
) -> (Resolved, Resolved) {
    use RoundingMode::{TowardNegative as TN, TowardPositive as TP};

    // Bracket hypot(x, y), then the radicand (|x| + hypot)/2 (the /2 is exact),
    // then w = sqrt(radicand). Each lower bound rounds toward −∞ end to end and
    // each upper bound toward +∞, so [s_lo, s_hi] encloses the true w.
    let h_lo = a.hypot_round(b, w, TN).expect("w >= 1").0;
    let h_hi = a.hypot_round(b, w, TP).expect("w >= 1").0;
    let absa = a.abs();
    let r_lo = absa.add(&h_lo, TN).0.scale_by_pow2(-1).0;
    let r_hi = absa.add(&h_hi, TP).0.scale_by_pow2(-1).0;
    let s_lo = r_lo.sqrt_round(w, TN).expect("w >= 1").0;
    let s_hi = r_hi.sqrt_round(w, TP).expect("w >= 1").0;
    let two_s_lo = s_lo.scale_by_pow2(1).0;
    let two_s_hi = s_hi.scale_by_pow2(1).0;

    let x_neg = a.is_sign_negative() && !a.is_zero();
    let (u_lo, u_hi, v_lo, v_hi) = if !x_neg {
        // x ≥ 0: u = w; v = y/(2w), sign(v) = sign(y), |v| decreasing in w.
        let (v_lo, v_hi) = if b.is_sign_negative() {
            (
                b.div_round(&two_s_lo, w, TN).expect("w >= 1").0,
                b.div_round(&two_s_hi, w, TP).expect("w >= 1").0,
            )
        } else {
            (
                b.div_round(&two_s_hi, w, TN).expect("w >= 1").0,
                b.div_round(&two_s_lo, w, TP).expect("w >= 1").0,
            )
        };
        (s_lo, s_hi, v_lo, v_hi)
    } else {
        // x < 0: |v| = w, v = copysign(w, y); u = |y|/(2w), positive,
        // decreasing in w.
        let absb = b.abs();
        let u_lo = absb.div_round(&two_s_hi, w, TN).expect("w >= 1").0;
        let u_hi = absb.div_round(&two_s_lo, w, TP).expect("w >= 1").0;
        let (v_lo, v_hi) = if b.is_sign_negative() {
            (s_hi.negated(), s_lo.negated())
        } else {
            (s_lo.clone(), s_hi.clone())
        };
        (u_lo, u_hi, v_lo, v_hi)
    };

    (
        resolve_bracket(&u_lo, &u_hi, p, mode),
        resolve_bracket(&v_lo, &v_hi, p, mode),
    )
}

fn pos_inf(p: u32) -> BigFloat {
    BigFloat::try_new_infinity(Sign::Positive, p).expect("p >= 1")
}

fn nan(p: u32) -> BigFloat {
    BigFloat::try_new_quiet_nan(Sign::Positive, p, &[]).expect("p >= 1")
}

/// `copysign(+∞, src)`: a signed infinity carrying `src`'s sign.
fn inf_signed(src: &BigFloat, p: u32) -> BigFloat {
    BigFloat::try_new_infinity(Sign::Positive, p)
        .expect("p >= 1")
        .copysign(src)
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
    fn eq(v: &BigFloat, n: i64) -> bool {
        matches!(v.partial_cmp(&bf(n)).0, Some(Ordering::Equal))
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

    #[test]
    fn csqrt_of_gaussian_integer_is_exact() {
        // csqrt(3 + 4i) = 2 + i, exact (no forced INEXACT).
        let (re, im, s) = csqrt_big(&bf(3), &bf(4), 64, NE);
        assert!(eq(&re, 2) && eq(&im, 1));
        assert!(!s.inexact(), "a Gaussian-integer root is exact");
    }

    #[test]
    fn csqrt_exact_roots_across_all_sign_branches() {
        // Gaussian-integer roots, independently checkable, one per Kahan
        // branch and v-bracket sign case:
        //   (3 − 2i)² =  5 − 12i  → x>0, y<0   (the b<0 division branch)
        //   (2 − 3i)² = −5 − 12i  → x<0, y<0   (the x<0, b<0 copysign branch)
        //   (2 + 3i)² = −5 + 12i  → x<0, y>0   (the x<0, b>0 copysign branch)
        for (zr, zi, wr, wi) in [(5, -12, 3, -2), (-5, -12, 2, -3), (-5, 12, 2, 3)] {
            let (re, im, s) = csqrt_big(&bf(zr), &bf(zi), 64, NE);
            assert!(eq(&re, wr) && eq(&im, wi), "csqrt({zr}+{zi}i) = {re}+{im}i");
            assert!(!s.inexact(), "Gaussian-integer root is exact");
        }
    }

    #[test]
    fn csqrt_squared_round_trips() {
        // csqrt(z)² = z for a non-square z (to rounding): square the root and
        // compare. z = 2 + 3i; w = csqrt(z); w² should be ≈ 2 + 3i.
        let (u, v, _) = csqrt_big(&bf(2), &bf(3), 200, NE);
        // (u + vi)² = (u² − v²) + (2uv)i.
        let re = u.mul_sub_mul(&u, &v, &v, NE).0;
        let two_uv = u.mul(&v, NE).0.scale_by_pow2(1).0;
        // Within a few ulps of (2, 3).
        let near = |x: &BigFloat, n: i64| {
            let d = x.sub(&bf(n), NE).0.abs();
            matches!(
                d.partial_cmp(&bf(1).scale_by_pow2(-180).0).0,
                Some(Ordering::Less)
            )
        };
        assert!(near(&re, 2), "re(w²) ≈ 2, got {re}");
        assert!(near(&two_uv, 3), "im(w²) ≈ 3, got {two_uv}");
    }

    #[test]
    fn csqrt_negative_real_axis_signed_zero_branch() {
        // The branch cut: csqrt(−4 + 0i) = +2i, csqrt(−4 − 0i) = −2i.
        let (re_u, im_u, _) = csqrt_big(&bf(-4), &pz(), 64, NE);
        assert!(re_u.is_zero() && re_u.is_sign_positive());
        assert!(eq(&im_u, 2) && im_u.is_sign_positive(), "upper cut → +2i");
        let (re_l, im_l, _) = csqrt_big(&bf(-4), &nz(), 64, NE);
        assert!(re_l.is_zero() && re_l.is_sign_positive());
        assert!(eq(&im_l, -2) && im_l.is_sign_negative(), "lower cut → −2i");
    }

    #[test]
    fn csqrt_positive_real_axis_signed_zero_imaginary() {
        // csqrt(4 + 0i) = 2 + 0i; csqrt(4 − 0i) = 2 − 0i. The imaginary zero's
        // sign follows y (the defect the explicit copysign fixes; routing it
        // through a division would give +0 in four of five modes).
        for mode in [
            RoundingMode::NearestEven,
            RoundingMode::TowardZero,
            RoundingMode::TowardPositive,
            RoundingMode::TowardNegative,
            RoundingMode::NearestAway,
        ] {
            let (re, im, _) = csqrt_big(&bf(4), &nz(), 64, mode);
            assert!(eq(&re, 2));
            assert!(
                im.is_zero() && im.is_sign_negative(),
                "csqrt(4 − 0i) imaginary part must be −0 under {mode:?}"
            );
        }
    }

    #[test]
    fn csqrt_origin_is_signed_zero() {
        let (re, im, s) = csqrt_big(&pz(), &nz(), 64, NE);
        assert!(re.is_zero() && re.is_sign_positive());
        assert!(im.is_zero() && im.is_sign_negative());
        assert!(!s.inexact());
    }

    #[test]
    fn csqrt_imaginary_axis() {
        // csqrt(2i) = 1 + i (x = +0, y = 2 uses the x ≥ 0 branch).
        let (re, im, _) = csqrt_big(&pz(), &bf(2), 200, NE);
        let near1 = |x: &BigFloat| {
            let d = x.sub(&bf(1), NE).0.abs();
            matches!(
                d.partial_cmp(&bf(1).scale_by_pow2(-180).0).0,
                Some(Ordering::Less)
            )
        };
        assert!(near1(&re) && near1(&im), "csqrt(2i) ≈ 1 + i");
    }

    #[test]
    fn csqrt_infinite_imaginary_dominates() {
        // csqrt(x ± ∞i) = +∞ ± ∞i for every x, including −∞ and NaN.
        for x in [&bf(3), &ninf(), &qnan()] {
            let (re, im, _) = csqrt_big(x, &pinf(), 64, NE);
            assert!(re.is_infinite() && re.is_sign_positive());
            assert!(im.is_infinite() && im.is_sign_positive());
            let (_, im2, _) = csqrt_big(x, &ninf(), 64, NE);
            assert!(im2.is_infinite() && im2.is_sign_negative());
        }
    }

    #[test]
    fn csqrt_positive_infinity_real() {
        // csqrt(+∞ + 3i) = +∞ + 0i; csqrt(+∞ + NaN·i) = +∞ + NaN·i.
        let (re, im, _) = csqrt_big(&pinf(), &bf(3), 64, NE);
        assert!(re.is_infinite() && re.is_sign_positive() && im.is_zero());
        let (re2, im2, _) = csqrt_big(&pinf(), &qnan(), 64, NE);
        assert!(re2.is_infinite() && im2.is_nan());
    }

    #[test]
    fn csqrt_negative_infinity_real() {
        // csqrt(−∞ + 3i) = +0 + ∞i; csqrt(−∞ + NaN·i) = NaN + ∞i [REP sign].
        let (re, im, _) = csqrt_big(&ninf(), &bf(3), 64, NE);
        assert!(re.is_zero() && re.is_sign_positive());
        assert!(im.is_infinite() && im.is_sign_positive());
        let (re2, im2, _) = csqrt_big(&ninf(), &qnan(), 64, NE);
        assert!(re2.is_nan() && im2.is_infinite());
    }

    #[test]
    fn csqrt_nan_without_infinity_is_nan() {
        let (re, im, s) = csqrt_big(&qnan(), &bf(2), 64, NE);
        assert!(re.is_nan() && im.is_nan());
        assert!(!s.invalid(), "a quiet NaN propagates without INVALID");
    }

    #[test]
    fn csqrt_signaling_nan_raises_invalid() {
        let snan = BigFloat::try_new_signaling_nan(Sign::Positive, 64, &[]).unwrap();
        let (_, _, s) = csqrt_big(&snan, &bf(2), 64, NE);
        assert!(s.invalid());
    }

    #[test]
    fn csqrt_directed_modes_bracket_the_enclosure() {
        // csqrt(2 + 3i)'s real part is irrational, computed by the Q16
        // enclosure: directed rounding must give TN ≤ NE ≤ TP with TN < TP.
        let (re_ne, _, _) = csqrt_big(&bf(2), &bf(3), 53, NE);
        let (re_tn, _, _) = csqrt_big(&bf(2), &bf(3), 53, RoundingMode::TowardNegative);
        let (re_tp, _, _) = csqrt_big(&bf(2), &bf(3), 53, RoundingMode::TowardPositive);
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
