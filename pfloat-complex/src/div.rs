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

    // The exact operands the exactness residual recomputes from, per component
    // (re numerator `ac + bd`, im numerator `bc − ad`, shared denom `c² + d²`).
    let re_inputs = ExactInputs {
        a,
        b,
        c,
        d,
        is_real: true,
    };
    let im_inputs = ExactInputs {
        a,
        b,
        c,
        d,
        is_real: false,
    };

    let mut last: Option<(Resolved, Resolved, Status)> = None;
    for &guard in &GUARDS {
        let w = p.saturating_add(guard);
        use RoundingMode::{TowardNegative as TN, TowardPositive as TP};
        // Denominator c² + d² (≥ 0), bracketed; shared by re and im.
        let (den_lo, sd_lo) = mul_add_mul(c, c, d, d, w, TN);
        let (den_hi, sd_hi) = mul_add_mul(c, c, d, d, w, TP);
        // re numerator ac + bd; im numerator bc − ad.
        let (rn_lo, srn_lo) = mul_add_mul(a, c, b, d, w, TN);
        let (rn_hi, srn_hi) = mul_add_mul(a, c, b, d, w, TP);
        let (in_lo, sin_lo) = mul_sub_mul(b, c, a, d, w, TN);
        let (in_hi, sin_hi) = mul_sub_mul(b, c, a, d, w, TP);

        // A finite product of finite operands can still saturate the i64
        // exponent when a squared tiny divisor underflows below `i64::MIN`
        // (`c² + d²`, pf-bv2i): the true value is then unrepresentable, so the
        // directed pair no longer brackets it and the enclosure is unsound.
        // Carry the OVERFLOW/UNDERFLOW so the result is honestly not-OK rather
        // than a wrong value with a clean flag. Normal-magnitude operands never
        // saturate, so this is `OK` on the ordinary path (no spurious flag).
        let den_sat = saturation_flags(sd_lo | sd_hi);
        let re_sat = saturation_flags(srn_lo | srn_hi) | den_sat;
        let im_sat = saturation_flags(sin_lo | sin_hi) | den_sat;

        let re = resolve(&rn_lo, &rn_hi, &den_lo, &den_hi, p, mode, &re_inputs);
        let im = resolve(&in_lo, &in_hi, &den_lo, &den_hi, p, mode, &im_inputs);
        if re.converged && im.converged {
            return (re.value, im.value, re.status | im.status | re_sat | im_sat);
        }
        last = Some((re, im, re_sat | im_sat));
    }
    // Cap exhausted: best-effort with the highest-precision candidate.
    let (re, im, sat) = last.expect("GUARDS is non-empty");
    (re.value, im.value, re.status | im.status | sat)
}

/// The exact operands defining one quotient component, for the exactness
/// residual `N − r·D` (`N = ac+bd` when `is_real`, else `bc−ad`; `D = c²+d²`).
struct ExactInputs<'a> {
    a: &'a BigFloat,
    b: &'a BigFloat,
    c: &'a BigFloat,
    d: &'a BigFloat,
    is_real: bool,
}

/// Keep only the exponent-saturation flags (`OVERFLOW` / `UNDERFLOW`) from a
/// bracket-formation status; the INEXACT of a directed bracket rounding is
/// expected and must not ride into the result.
fn saturation_flags(s: Status) -> Status {
    let mut out = Status::OK;
    if s.overflow() {
        out |= Status::OVERFLOW;
    }
    if s.underflow() {
        out |= Status::UNDERFLOW;
    }
    out
}

#[inline]
fn mul_add_mul(
    a: &BigFloat,
    b: &BigFloat,
    c: &BigFloat,
    d: &BigFloat,
    w: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    a.mul_add_mul_round(b, c, d, w, mode).expect("w >= 1")
}

#[inline]
fn mul_sub_mul(
    a: &BigFloat,
    b: &BigFloat,
    c: &BigFloat,
    d: &BigFloat,
    w: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    a.mul_sub_mul_round(b, c, d, w, mode).expect("w >= 1")
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
    exact: &ExactInputs,
) -> Resolved {
    // Zero numerator over a positive denominator: the quotient is an exact
    // signed zero, but the zero's sign has two sources that must be told apart
    // (pf-pz9r).
    //
    // - The directed pair AGREES in sign (`[−0, −0]` or `[+0, +0]`): the
    //   numerator is an exactly-representable signed zero fixed by the inputs
    //   (e.g. a zero dividend part, `ac + bd` with `a = b = 0`). Since
    //   `D = c² + d² ≥ 0`, IEEE `(±0)/(+D) = ±0` takes the numerator's sign,
    //   NOT the rounding mode.
    // - The pair STRADDLES (`[−0, +0]`): the numerator is a genuine cancelling
    //   difference (the imaginary part of z/z, `bc − ad = 0`), whose zero sign
    //   is mode-determined — `−0` under TowardNegative, `+0` otherwise.
    //
    // (A zero numerator over a zero denominator is 0/0 = NaN, left to the
    // division below.)
    if n_lo.is_zero() && n_hi.is_zero() && !d_hi.is_zero() {
        let sign = if n_lo.is_sign_negative() == n_hi.is_sign_negative() {
            if n_lo.is_sign_negative() {
                Sign::Negative
            } else {
                Sign::Positive
            }
        } else if matches!(mode, RoundingMode::TowardNegative) {
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
        // The quotient is exact when the directed quotient bracket collapsed to
        // one representable value; but that under-reports exactness, because an
        // exact quotient can have a non-collapsed bracket when N and D are
        // separately inexact at the working precision yet their ratio is exact
        // (z/z, where N = D; pf-bv2i part a). Fall back to the sound residual
        // certificate `N − r·D = 0` in that case.
        let is_exact = bracket_is_exact(&q_lo, &q_hi, &lo_r)
            || quotient_is_exact(exact, &lo_r, n_lo, n_hi, d_lo, d_hi);
        let mut status = if is_exact {
            Status::OK
        } else {
            Status::INEXACT
        };
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

/// Sound exactness certificate for a converged component whose directed
/// quotient bracket did not collapse. The rounded value `r` is the exact
/// quotient iff the residual `N − r·D = 0`, with `N` the component numerator
/// (`ac+bd` or `bc−ad`) and `D = c² + d²`.
///
/// A cheap screen first brackets the residual at the working precision from the
/// numerator/denominator pairs already in hand and bails unless it can contain
/// zero (the common inexact case, where `|N − r·D| ≈ ulp(r)·D` is far from
/// zero). Only then is the exact residual formed at a generous bounded
/// precision; if any step cannot be represented exactly (INEXACT or exponent
/// saturation), the check conservatively returns `false` — a spurious INEXACT,
/// never a false OK. That conservatism is the same measure-zero caveat
/// ADR-0090 already accepts for the value on a wild exponent gap.
fn quotient_is_exact(
    exact: &ExactInputs,
    r: &BigFloat,
    n_lo: &BigFloat,
    n_hi: &BigFloat,
    d_lo: &BigFloat,
    d_hi: &BigFloat,
) -> bool {
    // Only a finite nonzero `r` can be an exact quotient here (the zero case is
    // the short-circuit above; NaN/Inf are handled by the caller's dispatch).
    if !r.is_finite() || r.is_zero() {
        return false;
    }
    let w = n_lo.precision().max(d_lo.precision());
    if residual_bracket_excludes_zero(r, n_lo, n_hi, d_lo, d_hi, w) {
        return false;
    }
    exact_residual_is_zero(exact, r)
}

/// The residual `R = N − r·D` bracketed at precision `w` from the numerator
/// pair `[n_lo, n_hi] ∋ N` and the denominator pair `[d_lo, d_hi] ∋ D` (with
/// `D ≥ 0`). Returns `true` when the bracket wholly excludes zero, which proves
/// `R ≠ 0` (the quotient is not exact) without the expensive exact residual.
fn residual_bracket_excludes_zero(
    r: &BigFloat,
    n_lo: &BigFloat,
    n_hi: &BigFloat,
    d_lo: &BigFloat,
    d_hi: &BigFloat,
    w: u32,
) -> bool {
    use RoundingMode::{TowardNegative as TN, TowardPositive as TP};
    // r·D bounds. D ∈ [d_lo, d_hi] with D ≥ 0, so multiplying by a negative r
    // flips the endpoints.
    let (rd_lo, rd_hi) = if r.is_sign_negative() {
        (
            r.mul_round(d_hi, w, TN).expect("w >= 1").0,
            r.mul_round(d_lo, w, TP).expect("w >= 1").0,
        )
    } else {
        (
            r.mul_round(d_lo, w, TN).expect("w >= 1").0,
            r.mul_round(d_hi, w, TP).expect("w >= 1").0,
        )
    };
    // R = N − r·D ∈ [n_lo − rd_hi, n_hi − rd_lo].
    let r_lo = n_lo.sub_round(&rd_hi, w, TN).expect("w >= 1").0;
    let r_hi = n_hi.sub_round(&rd_lo, w, TP).expect("w >= 1").0;
    // Wholly positive (lower bound > 0) or wholly negative (upper bound < 0).
    (r_lo.is_sign_positive() && !r_lo.is_zero()) || (r_hi.is_sign_negative() && !r_hi.is_zero())
}

/// Form the exact residual `N − r·D` and test it against zero. Every step runs
/// at a generous bounded precision and is checked for exactness; any INEXACT or
/// exponent saturation makes the certificate bail (`false`). A `true` result
/// means `N − r·D = 0` exactly, hence `N/D = r` exactly.
fn exact_residual_is_zero(exact: &ExactInputs, r: &BigFloat) -> bool {
    use RoundingMode::NearestEven as NE;
    let base = exact
        .a
        .precision()
        .max(exact.b.precision())
        .max(exact.c.precision())
        .max(exact.d.precision())
        .max(r.precision());
    let w = base.saturating_mul(2).saturating_add(4096);

    let (n, s_n) = if exact.is_real {
        exact.a.mul_add_mul_round(exact.c, exact.b, exact.d, w, NE)
    } else {
        exact.b.mul_sub_mul_round(exact.c, exact.a, exact.d, w, NE)
    }
    .expect("w >= 1");
    let (den, s_d) = exact
        .c
        .mul_add_mul_round(exact.c, exact.d, exact.d, w, NE)
        .expect("w >= 1");
    if !is_exact_status(s_n) || !is_exact_status(s_d) {
        return false;
    }
    let (rd, s_rd) = r.mul_round(&den, w, NE).expect("w >= 1");
    if !is_exact_status(s_rd) {
        return false;
    }
    let (residual, s_res) = n.sub_round(&rd, w, NE).expect("w >= 1");
    if !is_exact_status(s_res) {
        return false;
    }
    residual.is_zero()
}

/// A fused/round step was exact when it neither rounded nor saturated the
/// exponent; only then does its value certify the exact residual.
fn is_exact_status(s: Status) -> bool {
    !s.inexact() && !s.overflow() && !s.underflow() && !s.invalid()
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
