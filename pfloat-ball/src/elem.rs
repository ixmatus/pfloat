//! Ball elementary functions: a thin sound enclosure over each pfloat
//! scalar kernel, gated per pfloat family.
//!
//! Three enclosure shapes cover the surface:
//!
//! - **Monotonic** (`exp`, `ln`, `sinh`, `atan`, `cbrt`, …): evaluate the
//!   correctly-rounded kernel at the outward-rounded input-interval
//!   endpoints and build the ball from that output interval. Tight and
//!   sound for any monotonic function.
//! - **Lipschitz-1** (`sin`, `cos`): the midpoint is the kernel result
//!   and the radius is `rad_mid + ra`, because `|f'| ≤ 1` makes the
//!   propagated error at most the input radius.
//! - **Composed** (`tan = sin / cos`, `cosh` via its even symmetry):
//!   built from the above plus ball arithmetic.
//!
//! Domain edges are handled soundly: where the function is unbounded at a
//! boundary the ball reaches (`ln` at `0`, `atanh` at `±1`) the result is
//! the entire real line; where it is bounded but the input strays out of
//! domain (`asin`/`acos` past `±1`, `acosh` below `1`) the result is
//! sound over the in-domain part and raises [`Status::INVALID`]. An
//! overflowing endpoint collapses to entire rather than panicking.

use pfloat::{RoundingMode, Status};

use crate::arith::half_spread;
use crate::ball::Ball;
use crate::mag::Mag;
use crate::scalar::RealScalar;

const NE: RoundingMode = RoundingMode::NearestEven;
const TP: RoundingMode = RoundingMode::TowardPositive;
const TN: RoundingMode = RoundingMode::TowardNegative;

type Kernel<T> = fn(&T, RoundingMode) -> (T, Status);

/// `(alo, ahi)`: the input ball's interval endpoints, rounded outward.
/// Precondition: the ball is not entire (so both are finite).
fn input_endpoints<T: RealScalar>(b: &Ball<T>) -> (T, T) {
    let a = b.midpoint();
    let ra = T::radius_to_scalar(b.radius());
    (a.sub(&ra, TN).0, a.add(&ra, TP).0)
}

/// Assemble the output ball from a directed endpoint pair `[lo, hi]`,
/// threading the kernel `status` (the OR of the two directed evaluations)
/// and flagging a *degenerate* enclosure with [`Status::OVERFLOW`].
///
/// A directed endpoint kernel raises `OVERFLOW` only at pfloat's `i64`
/// exponent rim (the no-emax contract, ADR-0099): the returned value is a
/// *finite* saturation whose true image is larger, so a bounded ball
/// built from it would UNDER-cover — a Law-1 hole of the same family
/// ADR-0099 closed for `mul`/`div`. The only sound representable response
/// is to widen the enclosure to the entire line. A non-finite (`±∞`)
/// endpoint reaches the same conclusion through [`Ball::from_interval`]'s
/// rejection. Either way the enclosure degenerates to unbounded and
/// `OVERFLOW` is surfaced; the `INEXACT`/`UNDERFLOW` bits pass through the
/// OR-monoid untouched (Law 5). Domain-driven `entire` results are
/// flagged `INVALID` by their callers before reaching here and are never
/// reclassified.
fn finish_enclosure<T: RealScalar>(prec: u32, lo: &T, hi: &T, status: Status) -> (Ball<T>, Status) {
    if status.overflow() {
        return (Ball::entire(prec), status | Status::OVERFLOW);
    }
    match Ball::from_interval(lo, hi) {
        Ok(ball) => (ball, status),
        Err(_) => (Ball::entire(prec), status | Status::OVERFLOW),
    }
}

/// The ball enclosing an increasing function's image over `[alo, ahi]`:
/// `[f(alo)↓, f(ahi)↑]`, with the kernel status threaded and a degenerate
/// (unbounded) enclosure flagged (see [`finish_enclosure`]).
fn enclose_increasing<T: RealScalar>(
    prec: u32,
    alo: &T,
    ahi: &T,
    k: Kernel<T>,
) -> (Ball<T>, Status) {
    let (lo, st_lo) = k(alo, TN);
    let (hi, st_hi) = k(ahi, TP);
    finish_enclosure(prec, &lo, &hi, st_lo | st_hi)
}

/// The ball enclosing a decreasing function's image over `[alo, ahi]`:
/// `[f(ahi)↓, f(alo)↑]`, with the kernel status threaded.
#[cfg(feature = "trig")]
fn enclose_decreasing<T: RealScalar>(
    prec: u32,
    alo: &T,
    ahi: &T,
    k: Kernel<T>,
) -> (Ball<T>, Status) {
    let (lo, st_lo) = k(ahi, TN);
    let (hi, st_hi) = k(alo, TP);
    finish_enclosure(prec, &lo, &hi, st_lo | st_hi)
}

/// A total monotonic-increasing function, enclosed by its endpoints.
fn increasing<T: RealScalar>(b: &Ball<T>, k: Kernel<T>) -> (Ball<T>, Status) {
    let prec = b.precision();
    if b.is_entire() {
        return (Ball::entire(prec), Status::OK);
    }
    let (alo, ahi) = input_endpoints(b);
    enclose_increasing(prec, &alo, &ahi, k)
}

/// A 1-Lipschitz function (`|f'| ≤ 1`): midpoint plus `rad_mid + ra`.
///
/// The output is bounded (`sin`/`cos` land in `[−1, 1]`), so the kernel
/// never overflows: the midpoint's round-to-nearest status is the whole
/// story and is threaded straight through (Law 5).
#[cfg(feature = "trig")]
fn lipschitz1<T: RealScalar>(b: &Ball<T>, k: Kernel<T>) -> (Ball<T>, Status) {
    let prec = b.precision();
    if b.is_entire() {
        return (Ball::entire(prec), Status::OK);
    }
    let a = b.midpoint();
    let (mid, status) = k(a, NE);
    let rad_mid = half_spread(&k(a, TN).0, &k(a, TP).0);
    let ra = T::radius_to_scalar(b.radius());
    let rad = rad_mid.add(&ra, TP).0.magnitude_to_mag();
    (Ball::from_parts(mid, rad), status)
}

/// `+1` and `-1` at the ball's precision, for domain clamping.
fn one_and_neg_one<T: RealScalar>() -> (T, T) {
    let one = T::radius_to_scalar(Mag::from_pow2(0));
    let neg_one = one.negated();
    (one, neg_one)
}

fn lt<T: RealScalar>(x: &T, y: &T) -> bool {
    x.compare(y).0 == Some(core::cmp::Ordering::Less)
}

// ---------- exp-log family ----------

impl<T: RealScalar> Ball<T> {
    /// `e^self`. Increasing and total.
    #[must_use]
    pub fn exp(&self) -> (Self, Status) {
        increasing(self, T::exp)
    }

    /// `e^self − 1`. Increasing and total.
    #[must_use]
    pub fn expm1(&self) -> (Self, Status) {
        increasing(self, T::expm1)
    }

    /// `2^self`. Increasing and total.
    #[must_use]
    pub fn exp2(&self) -> (Self, Status) {
        increasing(self, T::exp2)
    }

    /// `10^self`. Increasing and total.
    #[must_use]
    pub fn exp10(&self) -> (Self, Status) {
        increasing(self, T::exp10)
    }

    /// `ln(self)`. Increasing on `(0, ∞)`; a ball reaching `0` is
    /// unbounded below, so the result is entire with [`Status::INVALID`].
    #[must_use]
    pub fn ln(&self) -> (Self, Status) {
        self.log_family(T::ln)
    }

    /// `log2(self)`. See [`ln`](Self::ln).
    #[must_use]
    pub fn log2(&self) -> (Self, Status) {
        self.log_family(T::log2)
    }

    /// `log10(self)`. See [`ln`](Self::ln).
    #[must_use]
    pub fn log10(&self) -> (Self, Status) {
        self.log_family(T::log10)
    }

    fn log_family(&self, k: Kernel<T>) -> (Self, Status) {
        let prec = self.precision();
        if self.is_entire() {
            return (Ball::entire(prec), Status::INVALID);
        }
        let (alo, ahi) = input_endpoints(self);
        // Domain (0, ∞): alo ≤ 0 means the ball touches or crosses 0,
        // where ln → −∞, so the enclosure is unbounded below.
        if alo.is_zero() || alo.is_sign_negative() {
            return (Ball::entire(prec), Status::INVALID);
        }
        enclose_increasing(prec, &alo, &ahi, k)
    }

    /// `ln(1 + self)`. Increasing on `(−1, ∞)`; a ball reaching `−1` is
    /// unbounded below → entire + [`Status::INVALID`].
    #[must_use]
    pub fn log1p(&self) -> (Self, Status) {
        let prec = self.precision();
        if self.is_entire() {
            return (Ball::entire(prec), Status::INVALID);
        }
        let (alo, ahi) = input_endpoints(self);
        // Domain x > −1: alo + 1 ≤ 0 ⇔ alo ≤ −1.
        let neg_one = one_and_neg_one::<T>().1;
        if !lt(&neg_one, &alo) {
            return (Ball::entire(prec), Status::INVALID);
        }
        enclose_increasing(prec, &alo, &ahi, T::log1p)
    }

    /// `sinh(self)`. Increasing and total.
    #[must_use]
    pub fn sinh(&self) -> (Self, Status) {
        increasing(self, T::sinh)
    }

    /// `tanh(self)`. Increasing and total (output in `(−1, 1)`).
    #[must_use]
    pub fn tanh(&self) -> (Self, Status) {
        increasing(self, T::tanh)
    }

    /// `asinh(self)`. Increasing and total.
    #[must_use]
    pub fn asinh(&self) -> (Self, Status) {
        increasing(self, T::asinh)
    }

    /// `cosh(self)`. Even with its minimum `1` at `0`: enclosed over the
    /// magnitude interval `[min|x|, max|x|]`.
    #[must_use]
    pub fn cosh(&self) -> (Self, Status) {
        let prec = self.precision();
        if self.is_entire() {
            return (Ball::entire(prec), Status::OK);
        }
        let (alo, ahi) = input_endpoints(self);
        let neg = |v: &T| v.is_sign_negative() && !v.is_zero();
        let straddles = neg(&alo) && !neg(&ahi); // alo < 0 ≤ ahi
        let abs_lo = alo.abs();
        let abs_hi = ahi.abs();
        // max magnitude is the larger of the two absolute endpoints.
        let (max_mag, other) = if lt(&abs_lo, &abs_hi) {
            (abs_hi, abs_lo)
        } else {
            (abs_lo, abs_hi)
        };
        // min magnitude is 0 when the interval straddles 0, else the
        // smaller absolute endpoint.
        let min_mag = if straddles { T::zero(prec) } else { other };
        // cosh is increasing in |x|.
        enclose_increasing(prec, &min_mag, &max_mag, T::cosh)
    }

    /// `acosh(self)`. Increasing on `[1, ∞)`; the in-domain part is
    /// enclosed and a dip below `1` raises [`Status::INVALID`].
    #[must_use]
    pub fn acosh(&self) -> (Self, Status) {
        let prec = self.precision();
        if self.is_entire() {
            return (Ball::entire(prec), Status::INVALID);
        }
        let (alo, ahi) = input_endpoints(self);
        let one = one_and_neg_one::<T>().0;
        // Fully below the domain: ahi < 1.
        if lt(&ahi, &one) {
            return (Ball::entire(prec), Status::INVALID);
        }
        // Clamp the lower end up to 1; flag if it was below.
        let (clamped_lo, status) = if lt(&alo, &one) {
            (one, Status::INVALID)
        } else {
            (alo, Status::OK)
        };
        let (ball, kst) = enclose_increasing(prec, &clamped_lo, &ahi, T::acosh);
        (ball, status | kst)
    }

    /// `atanh(self)`. Increasing on `(−1, 1)`; a ball reaching `±1` is
    /// unbounded → entire + [`Status::INVALID`].
    #[must_use]
    pub fn atanh(&self) -> (Self, Status) {
        let prec = self.precision();
        if self.is_entire() {
            return (Ball::entire(prec), Status::INVALID);
        }
        let (alo, ahi) = input_endpoints(self);
        let (one, neg_one) = one_and_neg_one::<T>();
        // Domain (−1, 1): touching either bound makes the result unbounded.
        if !lt(&neg_one, &alo) || !lt(&ahi, &one) {
            return (Ball::entire(prec), Status::INVALID);
        }
        enclose_increasing(prec, &alo, &ahi, T::atanh)
    }

    /// `hypot(self, other) = √(self² + other²)`. 1-Lipschitz in each
    /// argument, so the radius is `rad_mid + ra + rb`.
    #[must_use]
    pub fn hypot(&self, other: &Self) -> (Self, Status) {
        let prec = self.precision().max(other.precision());
        if self.is_entire() || other.is_entire() {
            return (Ball::entire(prec), Status::OK);
        }
        let (a, b) = (self.midpoint(), other.midpoint());
        let (mid, status) = a.hypot(b, NE);
        let rad_mid = half_spread(&a.hypot(b, TN).0, &a.hypot(b, TP).0);
        let ra = T::radius_to_scalar(self.radius());
        let rb = T::radius_to_scalar(other.radius());
        let acc = rad_mid.add(&ra, TP).0.add(&rb, TP).0;
        (Ball::from_parts(mid, acc.magnitude_to_mag()), status)
    }
}

// ---------- trig family ----------

#[cfg(feature = "trig")]
impl<T: RealScalar> Ball<T> {
    /// `sin(self)`. 1-Lipschitz (the enclosure may exceed `[−1, 1]` but
    /// stays sound).
    #[must_use]
    pub fn sin(&self) -> (Self, Status) {
        lipschitz1(self, T::sin)
    }

    /// `cos(self)`. 1-Lipschitz.
    #[must_use]
    pub fn cos(&self) -> (Self, Status) {
        lipschitz1(self, T::cos)
    }

    /// `tan(self) = sin(self) / cos(self)`. Sound (the decorrelated
    /// quotient over-covers); a ball whose cosine straddles zero (a pole
    /// in range) yields the entire line with [`Status::DIV_BY_ZERO`].
    #[must_use]
    pub fn tan(&self) -> (Self, Status) {
        let (s, ss) = self.sin();
        let (c, sc) = self.cos();
        let (q, sq) = s.div(&c);
        (q, ss | sc | sq)
    }

    /// `asin(self)`. Increasing on `[−1, 1]`; out-of-domain input is
    /// clamped and raises [`Status::INVALID`].
    #[must_use]
    pub fn asin(&self) -> (Self, Status) {
        let prec = self.precision();
        if self.is_entire() {
            return (Ball::entire(prec), Status::INVALID);
        }
        let (alo, ahi) = input_endpoints(self);
        let (lo, hi, status) = match clamp_unit(&alo, &ahi) {
            Some(v) => v,
            None => return (Ball::entire(prec), Status::INVALID),
        };
        let (ball, kst) = enclose_increasing(prec, &lo, &hi, T::asin);
        (ball, status | kst)
    }

    /// `acos(self)`. Decreasing on `[−1, 1]`; out-of-domain input is
    /// clamped and raises [`Status::INVALID`].
    #[must_use]
    pub fn acos(&self) -> (Self, Status) {
        let prec = self.precision();
        if self.is_entire() {
            return (Ball::entire(prec), Status::INVALID);
        }
        let (alo, ahi) = input_endpoints(self);
        let (lo, hi, status) = match clamp_unit(&alo, &ahi) {
            Some(v) => v,
            None => return (Ball::entire(prec), Status::INVALID),
        };
        let (ball, kst) = enclose_decreasing(prec, &lo, &hi, T::acos);
        (ball, status | kst)
    }

    /// `atan(self)`. Increasing and total (output in `(−π/2, π/2)`).
    #[must_use]
    pub fn atan(&self) -> (Self, Status) {
        increasing(self, T::atan)
    }

    /// `atan2(self, x)`: the angle of the point `(x, self)`, output in
    /// `(−π, π]`. The radius uses the gradient magnitude `1/r` with
    /// `r = √(x² + self²)`: `rad = rad_mid + (rx + ry)/r_lo`, where `r_lo`
    /// is a lower bound on `r` over the input box. The box reaching the
    /// origin (`r_lo ≤ 0`) or straddling the branch cut (the negative
    /// `x`-axis) yields the entire line, the sound enclosure of the
    /// discontinuity.
    #[must_use]
    pub fn atan2(&self, x: &Self) -> (Self, Status) {
        let prec = self.precision().max(x.precision());
        if self.is_entire() || x.is_entire() {
            return (Ball::entire(prec), Status::INVALID);
        }
        let zero = T::zero(prec);

        // r_lo: a lower bound on √(x² + y²) over the box (the hypot ball's
        // lower endpoint). r_lo ≤ 0 ⇒ the box can reach the origin.
        let (r_ball, _) = self.hypot(x);
        let r_lo = r_ball.lower();
        if r_lo.is_zero() || r_lo.is_sign_negative() {
            return (Ball::entire(prec), Status::OK);
        }

        // Branch cut: x wholly negative AND y straddling 0 makes the angle
        // wrap ±π. The connected sound enclosure is the entire line.
        let (y_lo, y_hi) = (self.lower(), self.upper());
        let y_straddles_zero = !lt(&zero, &y_lo) && !lt(&y_hi, &zero);
        let x_wholly_negative = lt(&x.upper(), &zero);
        if x_wholly_negative && y_straddles_zero {
            return (Ball::entire(prec), Status::OK);
        }

        let (y0, x0) = (self.midpoint(), x.midpoint());
        let (mid, status) = y0.atan2(x0, NE);
        let rad_mid = half_spread(&y0.atan2(x0, TN).0, &y0.atan2(x0, TP).0);
        // (rx + ry)/r_lo: a point d from the centre changes the angle by
        // ≤ d/r_lo, and d ≤ √(rx²+ry²) ≤ rx + ry over the box.
        let ry = T::radius_to_scalar(self.radius());
        let rx = T::radius_to_scalar(x.radius());
        let num = ry.add(&rx, TP).0;
        let prop = num.div(&r_lo, TP).0;
        let rad = rad_mid.add(&prop, TP).0.magnitude_to_mag();
        (Ball::from_parts(mid, rad), status)
    }
}

/// Clamp `[alo, ahi]` to the unit domain `[−1, 1]`. Returns
/// `(clamped_lo, clamped_hi, status)`, or `None` when the interval lies
/// wholly outside `[−1, 1]`.
#[cfg(feature = "trig")]
fn clamp_unit<T: RealScalar>(alo: &T, ahi: &T) -> Option<(T, T, Status)> {
    let (one, neg_one) = one_and_neg_one::<T>();
    // Wholly outside: ahi < −1 or alo > 1.
    if lt(ahi, &neg_one) || lt(&one, alo) {
        return None;
    }
    let mut status = Status::OK;
    let lo = if lt(alo, &neg_one) {
        status = Status::INVALID;
        neg_one
    } else {
        alo.clone()
    };
    let hi = if lt(&one, ahi) {
        status = Status::INVALID;
        one
    } else {
        ahi.clone()
    };
    Some((lo, hi, status))
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;
    use pfloat::BigFloat;

    type B = Ball<BigFloat>;

    fn bf(n: i64, p: u32) -> BigFloat {
        BigFloat::try_from_i64_exact(n, p).unwrap()
    }
    fn ptb(n: i64, p: u32) -> B {
        Ball::point(bf(n, p)).unwrap()
    }
    fn ball(mid: i64, radius: Mag, p: u32) -> B {
        Ball::new(bf(mid, p), radius).unwrap()
    }
    fn contains(b: &B, x: &BigFloat) -> bool {
        b.lower().partial_cmp(x).0 != Some(Ordering::Greater)
            && b.upper().partial_cmp(x).0 != Some(Ordering::Less)
    }
    // High-precision scalar reference.
    fn refv(f: impl Fn(&BigFloat, RoundingMode) -> (BigFloat, Status), x: i64) -> BigFloat {
        f(&bf(x, 400), NE).0
    }

    #[test]
    fn exp_ln_round_trip_and_enclose() {
        let (e, _) = ptb(1, 64).exp();
        assert!(contains(&e, &refv(BigFloat::exp, 1))); // e ≈ 2.718
                                                        // exp is increasing: [0±ε] → around 1.
        let (e2, _) = ball(0, Mag::from_pow2(-30), 64).exp();
        assert!(contains(&e2, &bf(1, 64)));
        // ln of a positive ball encloses the truth. ln(10) is irrational,
        // so the kernel rounds: the (now-threaded, ADR-0116) status is
        // INEXACT — the normal correct outcome for an inexact ball op
        // (Law 5) — not OK and not a domain error.
        let (l, st) = ptb(10, 64).ln();
        assert!(st.inexact() && !st.invalid() && contains(&l, &refv(BigFloat::ln, 10)));
    }

    #[cfg(feature = "exp-log")]
    #[test]
    fn ln_of_zero_crossing_is_entire_invalid() {
        let (l, st) = ball(0, Mag::from_pow2(0), 64).ln(); // [-1, 1] crosses 0
        assert!(st.invalid() && l.is_entire());
        let (l2, st2) = ptb(0, 64).ln(); // exactly 0
        assert!(st2.invalid() && l2.is_entire());
    }

    #[cfg(feature = "exp-log")]
    #[test]
    fn cosh_straddling_zero_has_min_one() {
        // [-2, 2]: cosh min is cosh(0)=1, max cosh(2)≈3.76.
        let (c, _) = ball(0, Mag::from_pow2(1), 64).cosh();
        // lower must be ≤ 1 (contains the minimum), upper ≥ cosh(2).
        assert!(c.lower().partial_cmp(&bf(1, 64)).0 != Some(Ordering::Greater));
        assert!(contains(&c, &refv(BigFloat::cosh, 2)));
        assert!(contains(&c, &bf(1, 64))); // cosh(0)=1 is in range
    }

    #[cfg(feature = "exp-log")]
    #[test]
    fn acosh_domain_clamp() {
        // [0, 2] dips below 1: INVALID but sound over [1, 2].
        let (r, st) = ball(1, Mag::from_pow2(0), 64).acosh(); // [0, 2]
        assert!(st.invalid());
        assert!(contains(&r, &refv(BigFloat::acosh, 2)));
        // Wholly below 1: entire + INVALID.
        let (r2, st2) = ptb(0, 64).acosh();
        assert!(st2.invalid() && r2.is_entire());
    }

    #[cfg(feature = "exp-log")]
    #[test]
    fn hypot_encloses() {
        // hypot(3,4)=5 exactly.
        let (h, _) = ptb(3, 64).hypot(&ptb(4, 64));
        assert!(h.midpoint().partial_cmp(&bf(5, 64)).0 == Some(Ordering::Equal));
        // with radii, contains 5.
        let (h2, _) = ball(3, Mag::from_pow2(-2), 64).hypot(&ball(4, Mag::from_pow2(-2), 64));
        assert!(contains(&h2, &bf(5, 64)));
    }

    #[cfg(feature = "trig")]
    #[test]
    fn sin_cos_lipschitz_enclose() {
        // sin(0)=0; [0 ± 0.5] → sin in [-0.48, 0.48]-ish, contains 0.
        let (s, _) = ball(0, Mag::from_pow2(-1), 64).sin();
        assert!(contains(&s, &bf(0, 64)));
        // cos(0)=1; small ball around 0 contains 1.
        let (c, _) = ball(0, Mag::from_pow2(-10), 64).cos();
        assert!(contains(&c, &bf(1, 64)));
        // sin is 1-Lipschitz: the radius is ≥ the input radius.
        let big_r = ball(0, Mag::from_pow2(2), 64).sin(); // radius 4
        assert!(!big_r.0.is_exact());
    }

    #[cfg(feature = "trig")]
    #[test]
    fn asin_domain_clamp_and_acos() {
        // asin over [-2, 2] clamps to [-1,1] with INVALID.
        let (r, st) = ball(0, Mag::from_pow2(1), 64).asin(); // [-2, 2]
        assert!(st.invalid());
        assert!(contains(&r, &bf(0, 64))); // asin(0)=0
                                           // acos(1)=0, acos(-1)=π. acos over [-1,1] (radius 1 at mid 0).
        let (ac, _) = ball(0, Mag::from_pow2(0), 64).acos();
        assert!(contains(&ac, &refv(BigFloat::acos, 0))); // acos(0)=π/2
    }

    #[cfg(feature = "trig")]
    #[test]
    fn atan_encloses_and_total() {
        let (r, st) = ptb(1, 64).atan();
        // atan(1) = π/4 is irrational; the threaded kernel status
        // (ADR-0116) is INEXACT, not OK.
        assert!(st.inexact() && !st.invalid());
        assert!(contains(&r, &refv(BigFloat::atan, 1))); // atan(1)=π/4
    }
}
