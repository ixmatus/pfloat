//! Componentwise correctly-rounded complex division via a directed-pair
//! enclosure Ziv loop (ADR-0090).
//!
//! `(a + bi)/(c + di) = [(ac + bd) + (bc − ad)i] / (c² + d²)`. Each
//! component is a ratio of fused two-product expressions. Dividing a
//! separately-rounded numerator by a separately-rounded denominator carries
//! two roundings and is not correctly rounded; and forming the exact
//! numerator is infeasible when the two products differ wildly in magnitude
//! (the exact sum can exceed any representable precision). So bracket the
//! true numerator and denominator with their directed
//! (`TowardNegative` / `TowardPositive`) fused-two-product pairs at a
//! working precision above the output precision, form the quotient
//! interval, round both ends to the output precision under the target mode,
//! and accept when they agree — otherwise grow the working precision.
//!
//! The denominator `c² + d²` is non-negative; the quotient-interval sign
//! handling below relies on that.

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Sign, Status};

/// Working-precision guards, following pfloat's Ziv schedule (`ziv.rs`):
/// each iteration uses `working = p + guard`, capped at five iterations
/// (the MPFR measure-zero caveat for the rare unresolved hard-to-round
/// input in a directed mode).
const GUARDS: [u32; 5] = [64, 128, 256, 512, 1024];

/// Componentwise correctly-rounded `(a + bi)/(c + di)` at output precision
/// `p`, returning `(re, im, status)`. The generic complex divide bridges
/// through this with `RealScalar::to_big` / `from_big`, because the working
/// precision exceeds what a `FixedFloat<PREC>` can hold.
pub(crate) fn complex_div_big(
    a: &BigFloat,
    b: &BigFloat,
    c: &BigFloat,
    d: &BigFloat,
    p: u32,
    mode: RoundingMode,
) -> (BigFloat, BigFloat, Status) {
    // Annex G §G.5.1 infinity/NaN recovery runs first: a zero or infinite
    // divisor, an infinite dividend, or a NaN operand is resolved to its
    // mandated value here (ADR-0091), so the directed-pair Ziv loop below only
    // ever sees finite operands with a nonzero finite divisor.
    if let Some(special) = crate::specials::complex_div_special(a, b, c, d, p) {
        return special;
    }

    let mut last: Option<(Resolved, Resolved)> = None;
    for &guard in &GUARDS {
        let w = p.saturating_add(guard);
        // Denominator c² + d² (≥ 0), bracketed; shared by re and im.
        let den_lo = mul_add_mul(c, c, d, d, w, RoundingMode::TowardNegative);
        let den_hi = mul_add_mul(c, c, d, d, w, RoundingMode::TowardPositive);
        // re numerator ac + bd; im numerator bc − ad.
        let re = resolve(
            &mul_add_mul(a, c, b, d, w, RoundingMode::TowardNegative),
            &mul_add_mul(a, c, b, d, w, RoundingMode::TowardPositive),
            &den_lo,
            &den_hi,
            p,
            mode,
        );
        let im = resolve(
            &mul_sub_mul(b, c, a, d, w, RoundingMode::TowardNegative),
            &mul_sub_mul(b, c, a, d, w, RoundingMode::TowardPositive),
            &den_lo,
            &den_hi,
            p,
            mode,
        );
        if re.converged && im.converged {
            return (re.value, im.value, re.status | im.status);
        }
        last = Some((re, im));
    }
    // Cap exhausted: best-effort with the highest-precision candidate.
    let (re, im) = last.expect("GUARDS is non-empty");
    (re.value, im.value, re.status | im.status)
}

#[inline]
fn mul_add_mul(
    a: &BigFloat,
    b: &BigFloat,
    c: &BigFloat,
    d: &BigFloat,
    w: u32,
    mode: RoundingMode,
) -> BigFloat {
    a.mul_add_mul_round(b, c, d, w, mode).expect("w >= 1").0
}

#[inline]
fn mul_sub_mul(
    a: &BigFloat,
    b: &BigFloat,
    c: &BigFloat,
    d: &BigFloat,
    w: u32,
    mode: RoundingMode,
) -> BigFloat {
    a.mul_sub_mul_round(b, c, d, w, mode).expect("w >= 1").0
}

struct Resolved {
    value: BigFloat,
    converged: bool,
    status: Status,
}

/// One component: bracket the quotient `N / D` (with `D = c² + d² ≥ 0`),
/// round both ends to `p` under `mode`, and decide convergence and status.
fn resolve(
    n_lo: &BigFloat,
    n_hi: &BigFloat,
    d_lo: &BigFloat,
    d_hi: &BigFloat,
    p: u32,
    mode: RoundingMode,
) -> Resolved {
    // Exact-zero numerator over a positive denominator (e.g. the imaginary
    // part of z/z, where bc − ad cancels to exactly 0): the directed
    // numerator pair brackets ±0, and the sign-aware Ziv test never agrees on
    // the zero's sign, so handle it directly. The quotient is an exact signed
    // zero — positive except in TowardNegative, matching the IEEE sign of a
    // cancelling difference divided by a positive value. (A zero numerator
    // over a zero denominator is 0/0 = NaN, left to the division below.)
    if n_lo.is_zero() && n_hi.is_zero() && !d_hi.is_zero() {
        let sign = if matches!(mode, RoundingMode::TowardNegative) {
            Sign::Negative
        } else {
            Sign::Positive
        };
        return Resolved {
            value: BigFloat::try_new_zero(sign, p).expect("p >= 1"),
            converged: true,
            status: Status::OK,
        };
    }

    // Sign-aware quotient enclosure (D ≥ 0). The most-negative quotient
    // divides the lower numerator by the smaller denominator when the
    // numerator is negative and by the larger when it is non-negative; the
    // most-positive is dual. The naive `[n_lo/d_hi, n_hi/d_lo]` holds only
    // for a non-negative numerator, but the imaginary numerator `bc − ad`
    // is routinely negative.
    let w = d_lo.precision().max(n_lo.precision());
    let den_for_lo = if n_lo.is_sign_negative() { d_lo } else { d_hi };
    let den_for_hi = if n_hi.is_sign_negative() { d_hi } else { d_lo };
    let (q_lo, s_lo) = n_lo
        .div_round(den_for_lo, w, RoundingMode::TowardNegative)
        .expect("w >= 1");
    let (q_hi, s_hi) = n_hi
        .div_round(den_for_hi, w, RoundingMode::TowardPositive)
        .expect("w >= 1");
    let div_flags = s_lo | s_hi;

    let (lo_r, _) = q_lo.round_to_precision(p, mode).expect("p >= 1");
    let (hi_r, _) = q_hi.round_to_precision(p, mode).expect("p >= 1");

    // A NaN component (NaN input, or a 0/0 denominator) yields NaN;
    // converged iff both ends are NaN.
    if lo_r.is_nan() || hi_r.is_nan() {
        return Resolved {
            converged: lo_r.is_nan() && hi_r.is_nan(),
            status: special_flags(div_flags),
            value: lo_r,
        };
    }

    // Convergence needs equal value AND equal sign (to separate ±0, which
    // IEEE comparison treats as equal).
    let value_equal = matches!(lo_r.partial_cmp(&hi_r).0, Some(Ordering::Equal));
    let sign_equal = lo_r.is_sign_negative() == hi_r.is_sign_negative();
    if value_equal && sign_equal {
        let exact = bracket_is_exact(&q_lo, &q_hi, &lo_r);
        let mut status = if exact { Status::OK } else { Status::INEXACT };
        status |= special_flags(div_flags);
        Resolved {
            value: lo_r,
            converged: true,
            status,
        }
    } else {
        Resolved {
            value: lo_r,
            converged: false,
            status: Status::INEXACT | special_flags(div_flags),
        }
    }
}

/// The quotient is exact iff the directed bracket collapsed to one value and
/// that value is representable at the output precision (no rounding).
fn bracket_is_exact(q_lo: &BigFloat, q_hi: &BigFloat, lo_r: &BigFloat) -> bool {
    let collapsed = matches!(q_lo.partial_cmp(q_hi).0, Some(Ordering::Equal))
        && q_lo.is_sign_negative() == q_hi.is_sign_negative();
    let representable = matches!(lo_r.partial_cmp(q_lo).0, Some(Ordering::Equal))
        && lo_r.is_sign_negative() == q_lo.is_sign_negative();
    collapsed && representable
}

/// Keep only the DIV_BY_ZERO / INVALID flags from the directed divisions;
/// their INEXACT is about the bracket rounding, not the final result.
fn special_flags(div_flags: Status) -> Status {
    let mut s = Status::OK;
    if div_flags.div_by_zero() {
        s |= Status::DIV_BY_ZERO;
    }
    if div_flags.invalid() {
        s |= Status::INVALID;
    }
    s
}
