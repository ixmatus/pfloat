//! C99/C11 Annex G §G.5.1 infinity and NaN recovery for complex multiply and
//! divide (ADR-0091).
//!
//! Annex G §G.3 defines a complex value with at least one infinite part as a
//! *complex infinity*, even when the other part is NaN. The textbook
//! componentwise formulas do not preserve this: `(1 + 0i)·(∞ + ∞i)` evaluates
//! the cross products `1·∞ − 0·∞` and `1·∞ + 0·∞`, and the `0·∞` terms are
//! NaN, so the naive product collapses to `(NaN, NaN)` where Annex G mandates
//! an infinity. §G.5.1 specifies a recovery: when the naive result is
//! `(NaN, NaN)` but an operand is a complex infinity, *box* the infinite
//! operand (its infinite parts become signed ones, its finite or NaN parts
//! become signed zeros) and recompute against a literal infinity.
//!
//! Division has the dual problem at a zero divisor: `z / 0` for a finite
//! nonzero `z` is a complex infinity, not the componentwise `0/0 = NaN` that
//! C3 shipped. §G.5.1 specifies the directed-infinity recovery for the zero
//! divisor (D1), the infinite dividend over a finite divisor (D2), and the
//! finite dividend over an infinite divisor (D3).
//!
//! These routines are *derived from the Annex G algorithm text*, not copied
//! from any implementation; the helper decomposition and naming are this
//! crate's own. They run as an exact pre-dispatch in `BigFloat`, ahead of the
//! C3 directed-pair Ziv divide and the C3 fused multiply, so those kernels
//! only ever see finite operands (a nonzero finite divisor, for divide). Per
//! §G.5.1p5 and footnote 377 the *values* are normative; the raised flags are
//! not, so this code pins the values and lets the flags fall out best-effort.
//!
//! Overflow recovery (§G.5.1's third multiply branch, where a finite-operand
//! cross product overflows to an infinity that should be recovered) cannot
//! arise here: `BigFloat`'s exponent is an `i64` with saturating arithmetic,
//! so a product of two finite operands never overflows the exponent field
//! before the final round. The branch is named here and intentionally absent;
//! the exemption is the same one `hypot` relies on (ADR-0032).

use pfloat::{BigFloat, RoundingMode, Sign, Status};

/// A complex value is a complex infinity (§G.3) when at least one part is
/// infinite, regardless of the other part.
fn is_complex_inf(re: &BigFloat, im: &BigFloat) -> bool {
    re.is_infinite() || im.is_infinite()
}

/// A complex value is a complex zero when both parts are zero (either sign).
fn is_complex_zero(re: &BigFloat, im: &BigFloat) -> bool {
    re.is_zero() && im.is_zero()
}

/// `+∞` at precision `p`.
fn pos_inf(p: u32) -> BigFloat {
    BigFloat::try_new_infinity(Sign::Positive, p).expect("precision >= 1")
}

/// A quiet NaN at precision `p`.
fn nan(p: u32) -> BigFloat {
    BigFloat::try_new_quiet_nan(Sign::Positive, p, &[]).expect("precision >= 1")
}

/// `(NaN, NaN)` with `INVALID` raised — the unrecovered Annex G result for a
/// genuine `∞/∞`, `0/0`, or NaN-without-recoverable-infinity.
fn nan_pair(p: u32) -> (BigFloat, BigFloat, Status) {
    (nan(p), nan(p), Status::INVALID)
}

/// Round a recovery value (an infinity, NaN, or signed zero) to the output
/// precision `p`. The funnel only adjusts the precision field for these
/// non-normal classes, so it never rounds a real digit.
fn at_p(v: BigFloat, p: u32) -> BigFloat {
    v.round_to_precision(p, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}

/// Box an operand part for the §G.5.1 recovery: an infinite part becomes a
/// signed `1`, a finite or NaN part becomes a signed `0`. The sign is carried
/// from the original part (`copysign`), which is what makes the recovered
/// infinity land on the correct side.
fn box_part(x: &BigFloat, p: u32) -> BigFloat {
    let mag = if x.is_infinite() {
        BigFloat::try_from_i64_exact(1, p).expect("precision >= 1")
    } else {
        BigFloat::try_new_zero(Sign::Positive, p).expect("precision >= 1")
    };
    mag.copysign(x)
}

/// `(re, im) = (ac + bd, bc − ad)` (the division numerators) for already-boxed
/// or finite operands, each a single fused two-product (exact for boxed
/// operands).
fn div_numerators(a: &BigFloat, b: &BigFloat, c: &BigFloat, d: &BigFloat) -> (BigFloat, BigFloat) {
    let m = RoundingMode::NearestEven;
    (
        a.mul_add_mul(c, b, d, m).0, // ac + bd
        b.mul_sub_mul(c, a, d, m).0, // bc − ad
    )
}

/// `Status::INVALID` iff any of the four operand parts is a signaling NaN.
/// IEEE 754 §7.2 raises the invalid-operation exception for a signaling NaN
/// operand of any arithmetic operation, so this floor rides on top of every
/// §G.5.1 recovery row: it is the one flag the recovery must never drop, even
/// where an infinity or a boxed zero overrides the *value* (pf-hdq1). A quiet
/// NaN raises nothing here (it propagates silently).
fn signaling_invalid(a: &BigFloat, b: &BigFloat, c: &BigFloat, d: &BigFloat) -> Status {
    if a.is_signaling_nan() || b.is_signaling_nan() || c.is_signaling_nan() || d.is_signaling_nan()
    {
        Status::INVALID
    } else {
        Status::OK
    }
}

/// Annex G §G.5.1 pre-dispatch for `(a + bi)/(c + di)` at output precision
/// `p`. Returns `Some` for the recovered infinity/NaN rows, or `None` when the
/// operands are finite with a nonzero finite divisor (the case the C3
/// directed-pair Ziv divide handles). The signaling-NaN INVALID floor
/// (`signaling_invalid`) rides on every recovered row.
pub(crate) fn complex_div_special(
    a: &BigFloat,
    b: &BigFloat,
    c: &BigFloat,
    d: &BigFloat,
    p: u32,
) -> Option<(BigFloat, BigFloat, Status)> {
    let inv = signaling_invalid(a, b, c, d);
    complex_div_special_core(a, b, c, d, p).map(|(re, im, s)| (re, im, s | inv))
}

fn complex_div_special_core(
    a: &BigFloat,
    b: &BigFloat,
    c: &BigFloat,
    d: &BigFloat,
    p: u32,
) -> Option<(BigFloat, BigFloat, Status)> {
    let z_inf = is_complex_inf(a, b);
    let w_inf = is_complex_inf(c, d);
    let w_zero = is_complex_zero(c, d);

    // D1: any dividend over a complex-zero divisor. The directed infinity's
    // sign comes from the divisor's real part `c` ONLY (§G.5.1); `d` is never
    // consulted. A zero dividend part yields `∞·0 = NaN` for that part, still
    // a complex infinity by §G.3, and `0/0` reaches here and falls out as
    // `(NaN, NaN)` naturally.
    if w_zero {
        let dir = BigFloat::try_new_infinity(c.sign(), p).expect("precision >= 1");
        let m = RoundingMode::NearestEven;
        let (re, s_re) = dir.mul(a, m);
        let (im, s_im) = dir.mul(b, m);
        return Some((at_p(re, p), at_p(im, p), special_flags(s_re | s_im)));
    }

    // ∞ / ∞ is not recovered (Annex G leaves it `(NaN, NaN)`).
    if z_inf && w_inf {
        return Some(nan_pair(p));
    }

    // D2: complex-infinite dividend over a finite (nonzero) divisor. Box the
    // dividend, then scale the finite numerators by a literal infinity.
    if z_inf && c.is_finite() && d.is_finite() {
        let a = box_part(a, p);
        let b = box_part(b, p);
        let (n_re, n_im) = div_numerators(&a, &b, c, d);
        let (re, s_re) = pos_inf(p).mul(&n_re, RoundingMode::NearestEven);
        let (im, s_im) = pos_inf(p).mul(&n_im, RoundingMode::NearestEven);
        return Some((at_p(re, p), at_p(im, p), special_flags(s_re | s_im)));
    }

    // D3: finite dividend over a complex-infinite divisor. Box the divisor,
    // then scale the finite numerators by a literal zero. Fires whenever the
    // divisor has an infinite part, so an `∞ + NaN·i` divisor is resolved here
    // (its NaN part is boxed to a signed zero).
    if w_inf && a.is_finite() && b.is_finite() {
        let c = box_part(c, p);
        let d = box_part(d, p);
        let (n_re, n_im) = div_numerators(a, b, &c, &d);
        let zero = BigFloat::try_new_zero(Sign::Positive, p).expect("precision >= 1");
        let (re, s_re) = zero.mul(&n_re, RoundingMode::NearestEven);
        let (im, s_im) = zero.mul(&n_im, RoundingMode::NearestEven);
        return Some((at_p(re, p), at_p(im, p), special_flags(s_re | s_im)));
    }

    // A NaN operand with no recoverable infinity or zero divisor: not
    // recovered, `(NaN, NaN)`. The base status is OK: a *quiet* NaN propagates
    // silently (no INVALID). The `signaling_invalid` floor in the wrapper adds
    // INVALID iff the NaN was signaling (pf-hdq1). (Reaching here means at least
    // one part is NaN, since every all-finite-nonzero-divisor case returns
    // `None`.)
    if a.is_nan() || b.is_nan() || c.is_nan() || d.is_nan() {
        return Some((nan(p), nan(p), Status::OK));
    }

    // All finite, nonzero finite divisor: the C3 Ziv divide handles it.
    None
}

/// Annex G §G.5.1 infinity recovery for `(a + bi)·(c + di)`, applied only when
/// the C3 fused product already collapsed to `(NaN, NaN)`. Returns `Some` when
/// an operand is a complex infinity to recover, or `None` to keep the naive
/// `(NaN, NaN)` (a genuine NaN-without-infinity, including `∞·0`).
pub(crate) fn recover_mul(
    a: &BigFloat,
    b: &BigFloat,
    c: &BigFloat,
    d: &BigFloat,
    p: u32,
) -> Option<(BigFloat, BigFloat, Status)> {
    // The signaling-NaN INVALID floor (pf-hdq1): the boxing that recovers a
    // complex infinity turns a signaling-NaN part into a signed zero, silently
    // discarding the INVALID that sNaN operand must raise. Reinstate it on every
    // recovered row. (When no infinity is present `recover_mul` returns `None`
    // and the naive fused product already carries the sNaN INVALID.)
    let inv = signaling_invalid(a, b, c, d);
    recover_mul_core(a, b, c, d, p).map(|(re, im, s)| (re, im, s | inv))
}

fn recover_mul_core(
    a: &BigFloat,
    b: &BigFloat,
    c: &BigFloat,
    d: &BigFloat,
    p: u32,
) -> Option<(BigFloat, BigFloat, Status)> {
    let z_inf = is_complex_inf(a, b);
    let w_inf = is_complex_inf(c, d);
    if !z_inf && !w_inf {
        // No infinity to recover (finite operands cannot overflow the i64
        // exponent before the round, so the M-OVF branch never fires); the
        // naive `(NaN, NaN)` stands.
        return None;
    }
    // A complex infinity times a complex zero is a genuine `∞·0`: not
    // recovered, Annex G leaves it `(NaN, NaN)`.
    if (z_inf && is_complex_zero(c, d)) || (w_inf && is_complex_zero(a, b)) {
        return Some(nan_pair(p));
    }

    // Box the infinite operand(s) and the other operand's NaN parts.
    let mut a = a.clone();
    let mut b = b.clone();
    let mut c = c.clone();
    let mut d = d.clone();
    if z_inf {
        a = box_part(&a, p);
        b = box_part(&b, p);
        c = box_nan_to_zero(&c, p);
        d = box_nan_to_zero(&d, p);
    }
    if w_inf {
        c = box_part(&c, p);
        d = box_part(&d, p);
        a = box_nan_to_zero(&a, p);
        b = box_nan_to_zero(&b, p);
    }

    let m = RoundingMode::NearestEven;
    let n_re = a.mul_sub_mul(&c, &b, &d, m).0; // ac − bd
    let n_im = a.mul_add_mul(&d, &b, &c, m).0; // ad + bc
    let (re, s_re) = pos_inf(p).mul(&n_re, m);
    let (im, s_im) = pos_inf(p).mul(&n_im, m);
    Some((at_p(re, p), at_p(im, p), special_flags(s_re | s_im)))
}

/// A NaN part of the finite operand becomes a signed zero during recovery; a
/// non-NaN part is left as is.
fn box_nan_to_zero(x: &BigFloat, p: u32) -> BigFloat {
    if x.is_nan() {
        BigFloat::try_new_zero(Sign::Positive, p)
            .expect("precision >= 1")
            .copysign(x)
    } else {
        x.clone()
    }
}

/// Keep only the soundness-bearing flags from the recovery arithmetic (the
/// `∞·0`/`0·∞` `INVALID`); the recovery products never round, so there is no
/// `INEXACT` to carry, and §G.5.1 flags are best-effort regardless.
fn special_flags(s: Status) -> Status {
    let mut out = Status::OK;
    if s.invalid() {
        out |= Status::INVALID;
    }
    if s.div_by_zero() {
        out |= Status::DIV_BY_ZERO;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bf(n: i64) -> BigFloat {
        BigFloat::try_from_i64_exact(n, 64).unwrap()
    }
    fn pinf() -> BigFloat {
        BigFloat::try_new_infinity(Sign::Positive, 64).unwrap()
    }
    fn ninf() -> BigFloat {
        BigFloat::try_new_infinity(Sign::Negative, 64).unwrap()
    }
    fn pzero() -> BigFloat {
        BigFloat::try_new_zero(Sign::Positive, 64).unwrap()
    }
    fn is_pos_inf(v: &BigFloat) -> bool {
        v.is_infinite() && v.is_sign_positive()
    }
    fn is_neg_inf(v: &BigFloat) -> bool {
        v.is_infinite() && v.is_sign_negative()
    }

    #[test]
    fn div_finite_over_zero_is_directed_infinity() {
        // D1: (1 + 1i)/(0 + 0i). Direction from c = +0, so both parts are +∞.
        let (re, im, s) = complex_div_special(&bf(1), &bf(1), &pzero(), &pzero(), 64).unwrap();
        assert!(is_pos_inf(&re) && is_pos_inf(&im));
        assert!(!s.invalid(), "1·∞ is a clean infinity, no INVALID");
    }

    #[test]
    fn div_zero_real_part_over_zero_is_nan_component() {
        // D1: (0 + 1i)/(0 + 0i). re = ∞·0 = NaN (still a complex infinity by
        // §G.3 via the imaginary +∞); im = ∞·1 = +∞.
        let (re, im, s) = complex_div_special(&pzero(), &bf(1), &pzero(), &pzero(), 64).unwrap();
        assert!(re.is_nan());
        assert!(is_pos_inf(&im));
        assert!(s.invalid(), "∞·0 raises INVALID");
    }

    #[test]
    fn div_zero_over_zero_is_nan() {
        // 0/0 reaches D1 and falls out as (NaN, NaN).
        let (re, im, _) = complex_div_special(&pzero(), &pzero(), &pzero(), &pzero(), 64).unwrap();
        assert!(re.is_nan() && im.is_nan());
    }

    #[test]
    fn div_direction_from_divisor_real_sign() {
        // D1: divisor (−0 + 0i) flips the directed infinity to −∞ via c = −0.
        let neg_zero = BigFloat::try_new_zero(Sign::Negative, 64).unwrap();
        let (re, _, _) = complex_div_special(&bf(1), &bf(1), &neg_zero, &pzero(), 64).unwrap();
        assert!(is_neg_inf(&re), "direction follows c = −0");
    }

    #[test]
    fn div_infinite_over_finite_is_infinity() {
        // D2: (∞ + 0i)/(2 + 0i) = ∞. Boxed a = 1, b = +0; n_re = 1·2 = 2 > 0,
        // so re = +∞; n_im = 0·2 − 1·0 = 0, im = ∞·0 = NaN.
        let (re, im, _) = complex_div_special(&pinf(), &pzero(), &bf(2), &pzero(), 64).unwrap();
        assert!(is_pos_inf(&re));
        assert!(im.is_nan());
    }

    #[test]
    fn div_finite_over_infinite_is_zero() {
        // D3: (2 + 0i)/(∞ + 0i) = +0. Boxed c = 1, d = +0; numerators finite,
        // scaled by zero, so both parts are signed zeros.
        let (re, im, _) = complex_div_special(&bf(2), &pzero(), &pinf(), &pzero(), 64).unwrap();
        assert!(re.is_zero() && im.is_zero());
    }

    #[test]
    fn div_infinite_over_infinite_is_nan() {
        let (re, im, _) = complex_div_special(&pinf(), &pinf(), &pinf(), &pinf(), 64).unwrap();
        assert!(re.is_nan() && im.is_nan());
    }

    #[test]
    fn div_finite_over_inf_plus_nan_divisor_fires_d3() {
        // D3 fires on an (∞ + NaN·i) divisor: the infinite real part triggers
        // the branch and the NaN imaginary part is boxed to a signed zero, so
        // (2 + 0i)/(∞ + NaN·i) resolves to a signed-zero pair, not NaN.
        let (re, im, _) = complex_div_special(&bf(2), &pzero(), &pinf(), &nan(64), 64).unwrap();
        assert!(re.is_zero() && im.is_zero());
    }

    #[test]
    fn div_all_finite_nonzero_divisor_falls_through() {
        // The case the Ziv divide owns: no special dispatch.
        assert!(complex_div_special(&bf(1), &bf(2), &bf(3), &bf(4), 64).is_none());
    }

    #[test]
    fn div_nan_divisor_without_infinity_is_nan() {
        let q = nan(64);
        let (re, im, _) = complex_div_special(&bf(1), &bf(1), &q, &bf(2), 64).unwrap();
        assert!(re.is_nan() && im.is_nan());
    }

    #[test]
    fn mul_recovers_one_times_complex_infinity() {
        // M1: (1 + 0i)·(∞ + ∞i). Naive cross products give (NaN, NaN); recover
        // boxes (∞, ∞) → (1, 1), keeps (1, 0), and recomputes against ∞:
        // n_re = 1·1 − 0·1 = 1 → +∞, n_im = 1·1 + 0·1 = 1 → +∞.
        let (re, im, _) = recover_mul(&bf(1), &pzero(), &pinf(), &pinf(), 64).unwrap();
        assert!(is_pos_inf(&re) && is_pos_inf(&im));
    }

    #[test]
    fn mul_infinity_times_zero_is_not_recovered() {
        // (∞ + ∞i)·(0 + 0i) is a genuine ∞·0: stays (NaN, NaN).
        let (re, im, s) = recover_mul(&pinf(), &pinf(), &pzero(), &pzero(), 64).unwrap();
        assert!(re.is_nan() && im.is_nan());
        assert!(s.invalid());
    }

    #[test]
    fn mul_finite_operands_are_not_recovered() {
        // No infinity to recover: returns None so the naive (NaN, NaN) stands.
        assert!(recover_mul(&bf(1), &bf(2), &bf(3), &bf(4), 64).is_none());
    }

    #[test]
    fn mul_recovery_sign_follows_the_boxed_parts() {
        // (1 + 0i)·(−∞ + 0i): boxed (−∞, +0) → (−1, +0); n_re = 1·(−1) − 0 =
        // −1 → −∞; n_im = 1·0 + 0·(−1) = 0 → ∞·0 = NaN.
        let (re, im, _) = recover_mul(&bf(1), &pzero(), &ninf(), &pzero(), 64).unwrap();
        assert!(is_neg_inf(&re));
        assert!(im.is_nan());
    }

    #[test]
    fn at_p_normalizes_precision() {
        // A recovery infinity built at precision 64 normalizes to the output
        // precision (the funnel only adjusts the precision field here).
        let v = at_p(pinf(), 128);
        assert!(v.is_infinite() && v.is_sign_positive());
        assert_eq!(v.precision(), 128);
    }
}
