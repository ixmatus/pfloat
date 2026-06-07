//! Ball arithmetic: `add`, `sub`, `mul`, `div`, `sqrt` on the
//! directed-pair radius route. ADR-0077.
//!
//! Each binary kernel computes the midpoint with the correctly-rounded
//! pfloat scalar kernel (round-to-nearest) and bounds the radius as
//!
//! ```text
//!     rad = rad_mid + propagated_input_error
//! ```
//!
//! where `rad_mid = round_up((hi − lo)/2)` is the midpoint's own
//! rounding error read off the directed pair `lo = a op↓ b`,
//! `hi = a op↑ b`, and `propagated_input_error` bounds how far the true
//! result `f(x, y)` (for `x` in `[a]`, `y` in `[b]`) can stray from
//! `f(a.mid, b.mid)`. Soundness:
//!
//! ```text
//!     |f(x,y) − mid| ≤ |f(x,y) − f(a.mid,b.mid)| + |f(a.mid,b.mid) − mid|
//!                    ≤ propagated_input_error    + rad_mid    = rad
//! ```
//!
//! Every radius scalar operation rounds *outward* (numerators and the
//! whole radius toward `+∞`; the division denominator's `|y|` lower bound
//! toward `−∞`), and the final radius narrows up to [`Mag`]. Exactness in
//! produces exactness out: an exact result with exact inputs leaves
//! `rad = 0`.
//!
//! `sqrt` is unary and monotonic, so it is enclosed by evaluating the
//! correctly-rounded kernel at the input-interval endpoints (outward)
//! and building the ball from those.

use pfloat::{RoundingMode, Status};

use crate::ball::Ball;
use crate::scalar::RealScalar;

const NE: RoundingMode = RoundingMode::NearestEven;
const TP: RoundingMode = RoundingMode::TowardPositive;
const TN: RoundingMode = RoundingMode::TowardNegative;

/// `(hi − lo)/2` rounded up: an upper bound on the nearest-rounding error
/// `|mid − exact|`, where `lo`/`hi` are the directed pair bracketing the
/// exact result. `hi − lo` is a small exact multiple of the ulp, so the
/// halving is exact and the only over-estimate is at most one ulp.
fn half_spread<T: RealScalar>(lo: &T, hi: &T) -> T {
    let spread = hi.sub(lo, TP).0; // ≥ hi − lo ≥ 0
    spread.scale_by_pow2(-1).0 // exact ÷2
}

/// `a +↑ b` etc.: directed scalar helpers that keep the radius
/// accumulation an upper bound.
fn up<T: RealScalar>(x: &T, y: &T, op: fn(&T, &T, RoundingMode) -> (T, Status)) -> T {
    op(x, y, TP).0
}

impl<T: RealScalar> Ball<T> {
    /// `self + other`. The result encloses `{x + y : x ∈ self, y ∈ other}`.
    #[must_use]
    pub fn add(&self, other: &Self) -> (Self, Status) {
        let (a, b) = (self.midpoint(), other.midpoint());
        let (mid, status) = a.add(b, NE);
        let lo = a.add(b, TN).0;
        let hi = a.add(b, TP).0;
        let rad_mid = half_spread(&lo, &hi);
        // Propagated error for add: ra + rb.
        let ra = T::radius_to_scalar(self.radius());
        let rb = T::radius_to_scalar(other.radius());
        let acc = up(&rad_mid, &ra, T::add);
        let acc = up(&acc, &rb, T::add);
        (Self::from_parts(mid, acc.magnitude_to_mag()), status)
    }

    /// `self − other`. Encloses `{x − y : x ∈ self, y ∈ other}`.
    #[must_use]
    pub fn sub(&self, other: &Self) -> (Self, Status) {
        let (a, b) = (self.midpoint(), other.midpoint());
        let (mid, status) = a.sub(b, NE);
        let lo = a.sub(b, TN).0;
        let hi = a.sub(b, TP).0;
        let rad_mid = half_spread(&lo, &hi);
        // Propagated error for sub: ra + rb (|∂/∂x| = |∂/∂y| = 1).
        let ra = T::radius_to_scalar(self.radius());
        let rb = T::radius_to_scalar(other.radius());
        let acc = up(&rad_mid, &ra, T::add);
        let acc = up(&acc, &rb, T::add);
        (Self::from_parts(mid, acc.magnitude_to_mag()), status)
    }

    /// `self · other`. Encloses `{x · y : x ∈ self, y ∈ other}`.
    #[must_use]
    pub fn mul(&self, other: &Self) -> (Self, Status) {
        let (a, b) = (self.midpoint(), other.midpoint());
        let (mid, status) = a.mul(b, NE);
        let lo = a.mul(b, TN).0;
        let hi = a.mul(b, TP).0;
        let rad_mid = half_spread(&lo, &hi);
        // Propagated error for mul: |a|·rb + |b|·ra + ra·rb.
        let ra = T::radius_to_scalar(self.radius());
        let rb = T::radius_to_scalar(other.radius());
        let abs_a = a.abs();
        let abs_b = b.abs();
        let t1 = up(&abs_a, &rb, T::mul);
        let t2 = up(&abs_b, &ra, T::mul);
        let t3 = up(&ra, &rb, T::mul);
        let prop = up(&up(&t1, &t2, T::add), &t3, T::add);
        let rad = up(&rad_mid, &prop, T::add).magnitude_to_mag();
        (Self::from_parts(mid, rad), status)
    }

    /// `self / other`. Encloses `{x / y : x ∈ self, y ∈ other}`.
    ///
    /// If the divisor ball contains zero, the quotient is unbounded:
    /// returns the entire real line with [`Status::DIV_BY_ZERO`].
    #[must_use]
    pub fn div(&self, other: &Self) -> (Self, Status) {
        let (a, b) = (self.midpoint(), other.midpoint());
        let prec = a.precision().max(b.precision());
        let ra = T::radius_to_scalar(self.radius());
        let rb = T::radius_to_scalar(other.radius());
        let abs_b = b.abs();

        // blo = |b| − rb, a lower bound on |y| for y ∈ other (round down).
        // blo ≤ 0 means the divisor ball straddles zero.
        let blo = abs_b.sub(&rb, TN).0;
        if blo.is_zero() || blo.is_sign_negative() {
            return (Self::entire(prec), Status::DIV_BY_ZERO);
        }

        let (mid, status) = a.div(b, NE);
        let lo = a.div(b, TN).0;
        let hi = a.div(b, TP).0;
        let rad_mid = half_spread(&lo, &hi);

        // Propagated error: (|b|·ra + |a|·rb) / (blo·|b|), numerator up,
        // denominator down — both push the fraction to an upper bound.
        let abs_a = a.abs();
        let num = up(&up(&abs_b, &ra, T::mul), &up(&abs_a, &rb, T::mul), T::add);
        let den = blo.mul(&abs_b, TN).0; // lower bound on the denominator
        let prop = num.div(&den, TP).0; // upper bound on the quotient
        let rad = up(&rad_mid, &prop, T::add).magnitude_to_mag();
        (Self::from_parts(mid, rad), status)
    }

    /// `√self`. Encloses `{√x : x ∈ self, x ≥ 0}`.
    ///
    /// If the ball dips below zero it is sound over the in-domain part
    /// `[max(0, lo), hi]` and raises [`Status::INVALID`]; a wholly
    /// negative ball returns the entire real line with `INVALID`.
    #[must_use]
    pub fn sqrt(&self) -> (Self, Status) {
        let a = self.midpoint();
        let prec = a.precision();

        // An entire ball (rad = +∞, e.g. from dividing by a
        // zero-containing divisor) spans all reals: it includes negatives
        // (INVALID) and is unbounded above, so √ of it is again entire.
        // Guarding here also keeps the endpoint kernels below finite (the
        // radius scalar would otherwise be +∞ and make `ahi`/`alo`
        // non-finite, which `from_interval` rejects).
        if self.is_entire() {
            return (Self::entire(prec), Status::INVALID);
        }

        let ra = T::radius_to_scalar(self.radius());

        // Input interval endpoints, outward.
        let ahi = a.add(&ra, TP).0;
        let alo = a.sub(&ra, TN).0;

        let neg = |v: &T| v.is_sign_negative() && !v.is_zero();
        if neg(&ahi) {
            // Whole ball is negative: √ is defined nowhere on it.
            return (Self::entire(prec), Status::INVALID);
        }
        let mut status = Status::OK;
        // Clamp the lower end to the √ domain [0, ∞); a dip below zero is
        // INVALID but the enclosure of the in-domain part stays sound.
        let domain_lo = if neg(&alo) {
            status = Status::INVALID;
            T::zero(prec)
        } else {
            alo
        };

        let lo_out = domain_lo.sqrt(TN).0; // √lo rounded down
        let hi_out = ahi.sqrt(TP).0; // √hi rounded up
        let ball =
            Self::from_interval(&lo_out, &hi_out).expect("√ endpoints are finite and ordered");
        (ball, status)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mag::Mag;
    use core::cmp::Ordering;
    use pfloat::BigFloat;

    type B = Ball<BigFloat>;

    fn bf(n: i64, p: u32) -> BigFloat {
        BigFloat::try_from_i64_exact(n, p).unwrap()
    }
    fn pt(n: i64, p: u32) -> B {
        Ball::point(bf(n, p)).unwrap()
    }
    fn ball(mid: i64, radius: Mag, p: u32) -> B {
        Ball::new(bf(mid, p), radius).unwrap()
    }
    fn contains(b: &B, x: &BigFloat) -> bool {
        // lower ≤ x ≤ upper
        b.lower().partial_cmp(x).0 != Some(Ordering::Greater)
            && b.upper().partial_cmp(x).0 != Some(Ordering::Less)
    }

    // ---------- exact-in-exact-out ----------

    #[test]
    fn exact_add_sub_mul_are_exact() {
        let (s, st) = pt(2, 53).add(&pt(3, 53));
        assert!(st.is_ok() && s.is_exact());
        assert!(s.midpoint().partial_cmp(&bf(5, 53)).0 == Some(Ordering::Equal));

        let (d, _) = pt(10, 53).sub(&pt(3, 53));
        assert!(d.is_exact() && d.midpoint().partial_cmp(&bf(7, 53)).0 == Some(Ordering::Equal));

        let (m, _) = pt(6, 53).mul(&pt(7, 53));
        assert!(m.is_exact() && m.midpoint().partial_cmp(&bf(42, 53)).0 == Some(Ordering::Equal));
    }

    #[test]
    fn exact_div_when_divisible() {
        let (q, st) = pt(20, 53).div(&pt(4, 53));
        assert!(st.is_ok());
        assert!(q.is_exact());
        assert!(q.midpoint().partial_cmp(&bf(5, 53)).0 == Some(Ordering::Equal));
    }

    #[test]
    fn exact_sqrt_perfect_square() {
        let (r, st) = pt(9, 53).sqrt();
        assert!(st.is_ok() && r.is_exact());
        assert!(r.midpoint().partial_cmp(&bf(3, 53)).0 == Some(Ordering::Equal));
    }

    // ---------- inexact midpoint gets a positive radius ----------

    #[test]
    fn inexact_div_has_positive_radius_and_contains_truth() {
        // 1/3 at p=8 is not representable: the ball must be inexact and
        // contain the true 1/3 (checked via a 200-bit reference).
        let (q, _) = pt(1, 8).div(&pt(3, 8));
        assert!(!q.is_exact());
        let third = {
            let (t, _) = bf(1, 200).div(&bf(3, 200), NE);
            t
        };
        assert!(contains(&q, &third), "ball must enclose the true 1/3");
    }

    // ---------- propagated input radius (FTIA spot checks) ----------

    #[test]
    fn mul_propagates_input_radius() {
        // [10 ± 1] · [10 ± 1] must contain 9·9=81 .. 11·11=121, and 10·10.
        let a = ball(10, Mag::from_pow2(0), 53);
        let (p, _) = a.mul(&a);
        for &v in &[81i64, 100, 121, 99, 110] {
            assert!(contains(&p, &bf(v, 53)), "product ball must contain {v}");
        }
    }

    #[test]
    fn add_propagates_input_radius() {
        let a = ball(5, Mag::from_pow2(-1), 53); // [5 ± 0.5]
        let b = ball(7, Mag::from_pow2(-1), 53); // [7 ± 0.5]
        let (s, _) = a.add(&b); // [12 ± 1]
        for &v in &[11i64, 12, 13] {
            assert!(contains(&s, &bf(v, 53)));
        }
    }

    // ---------- division by a zero-containing ball ----------

    #[test]
    fn div_by_zero_containing_ball_is_entire() {
        let num = pt(1, 53);
        let denom = ball(0, Mag::from_pow2(0), 53); // [0 ± 1] contains 0
        let (q, st) = num.div(&denom);
        assert!(st.div_by_zero());
        assert!(q.is_entire());

        // Exact zero divisor too.
        let (q2, st2) = pt(1, 53).div(&pt(0, 53));
        assert!(st2.div_by_zero() && q2.is_entire());
    }

    // ---------- sqrt domain edges ----------

    #[test]
    fn sqrt_straddling_zero_is_invalid_but_sound() {
        // [1 ± 4] = [-3, 5]: dips below 0. INVALID, but encloses √[0,5].
        let b = ball(1, Mag::from_pow2(2), 53); // radius 4
        let (r, st) = b.sqrt();
        assert!(st.invalid());
        // √4 = 2 ∈ [0,5] must be enclosed; lower must be ≤ 0.
        assert!(contains(&r, &bf(2, 53)));
        assert!(r.lower().partial_cmp(&bf(0, 53)).0 != Some(Ordering::Greater));
    }

    #[test]
    fn sqrt_wholly_negative_is_entire_invalid() {
        let b = ball(-5, Mag::from_pow2(0), 53); // [-6, -4]
        let (r, st) = b.sqrt();
        assert!(st.invalid() && r.is_entire());
    }

    #[test]
    fn sqrt_of_entire_ball_does_not_panic() {
        // Regression (FTIA review): sqrt of an entire ball (rad = +∞) used
        // to panic in from_interval on the infinite endpoint.
        let entire = Ball::new(bf(0, 53), Mag::INFINITY).unwrap();
        let (r, st) = entire.sqrt();
        assert!(r.is_entire() && st.invalid());

        // Reachable path: divide by a zero-containing ball, then sqrt.
        let (q, _) = pt(1, 53).div(&ball(0, Mag::from_pow2(0), 53));
        assert!(q.is_entire());
        let (r2, _) = q.sqrt();
        assert!(r2.is_entire());
    }

    // ---------- radius soundness: enclose the true result for many points ----------

    #[test]
    fn ftia_spot_check_mul_over_grid() {
        // For balls a, b, every product x·y (x ∈ a, y ∈ b sampled at the
        // endpoints and midpoint) must lie in a.mul(b).
        let a = ball(3, Mag::from_pow2(-2), 53); // [3 ± 0.25]
        let b = ball(-4, Mag::from_pow2(-1), 53); // [-4 ± 0.5]
        let (p, _) = a.mul(&b);
        let xs = [a.lower(), a.midpoint().clone(), a.upper()];
        let ys = [b.lower(), b.midpoint().clone(), b.upper()];
        for x in &xs {
            for y in &ys {
                let (prod, _) = x.mul(y, NE);
                assert!(contains(&p, &prod), "mul ball must contain {prod:?}");
            }
        }
    }
}
