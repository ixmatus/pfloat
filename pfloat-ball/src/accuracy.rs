//! Accuracy introspection: reading the certified accuracy off a ball's
//! radius, and a precision-increasing retry that separates *effort* (the
//! working precision) from *certified accuracy* (the radius).
//!
//! On a ball the radius is the primary accuracy channel — a small
//! positive radius is the normal correct outcome of an inexact
//! operation, not a failure — so these are the methods a caller reaches
//! for instead of the IEEE status flags.

use pfloat::{BigFloat, RoundingMode, Status};

use crate::ball::Ball;
use crate::mag::Mag;
use crate::scalar::RealScalar;

const TP: RoundingMode = RoundingMode::TowardPositive;

impl<T: RealScalar> Ball<T> {
    /// The certified relative accuracy in bits, `≈ log2(|mid| / rad)`
    /// (the exponent of the midpoint minus the exponent of the radius,
    /// matching Arb's `arb_rel_accuracy_bits`).
    ///
    /// `i64::MAX` for an exact ball (`rad = 0`); `i64::MIN` for an entire
    /// ball or one whose midpoint is zero with a positive radius (no
    /// relative accuracy). The value is an estimate good to about one
    /// bit, not a strict bound.
    #[must_use]
    pub fn rel_accuracy_bits(&self) -> i64 {
        let rad_exp = match self.radius() {
            Mag::Zero => return i64::MAX,     // exact: infinite accuracy
            Mag::Infinity => return i64::MIN, // entire: none
            Mag::Finite { exponent, .. } => exponent,
        };
        match self.midpoint().exponent() {
            // Zero midpoint with a positive radius: no relative accuracy.
            None => i64::MIN,
            Some(mid_exp) => (i128::from(mid_exp) - i128::from(rad_exp))
                .clamp(i128::from(i64::MIN), i128::from(i64::MAX))
                as i64,
        }
    }

    /// The absolute error bound: the radius. Every point the ball denotes
    /// is within this magnitude of the midpoint.
    #[inline]
    #[must_use]
    pub fn abs_error(&self) -> Mag {
        self.radius()
    }

    /// An upper bound on the relative error `rad / |mid|`, as a [`Mag`].
    ///
    /// `0` for an exact ball; `+∞` for an entire ball or a zero-midpoint
    /// ball with positive radius (the relative error is unbounded).
    #[must_use]
    pub fn rel_error_bound(&self) -> Mag {
        match self.radius() {
            Mag::Zero => Mag::ZERO,
            Mag::Infinity => Mag::INFINITY,
            rad => {
                let abs_mid = self.midpoint().abs();
                if abs_mid.is_zero() {
                    return Mag::INFINITY;
                }
                // rad / |mid| rounded up: radius scalar up, divide toward +∞.
                let rad_s = T::radius_to_scalar(rad);
                rad_s.div(&abs_mid, TP).0.magnitude_to_mag()
            }
        }
    }
}

/// Re-evaluates a precision-parameterized ball computation at growing
/// working precision until it certifies at least `target_bits` of
/// relative accuracy, or `max_precision` is reached.
///
/// This is the effort-vs-accuracy separation: the caller writes the
/// computation once as a function of the working precision, and the
/// driver supplies whatever precision the certified radius demands.
/// `compute(p)` should build its inputs at precision `p` and return the
/// resulting ball and status. The returned ball is the first to reach
/// the target, or the `max_precision` result if the target is
/// unreachable (e.g. a genuinely entire result). An exact ball
/// short-circuits immediately.
///
/// Precision grows geometrically (×1.5, at least +32 bits) so the loop
/// terminates in `O(log(max_precision))` evaluations.
pub fn refine_to_accuracy<F>(
    target_bits: i64,
    start_precision: u32,
    max_precision: u32,
    mut compute: F,
) -> (Ball<BigFloat>, Status)
where
    F: FnMut(u32) -> (Ball<BigFloat>, Status),
{
    let max_precision = max_precision.max(1);
    let mut prec = start_precision.clamp(1, max_precision);
    loop {
        let (ball, status) = compute(prec);
        if ball.rel_accuracy_bits() >= target_bits || prec >= max_precision {
            return (ball, status);
        }
        let grow = prec.saturating_add(prec / 2).max(prec.saturating_add(32));
        prec = grow.min(max_precision);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    fn bf(n: i64, p: u32) -> BigFloat {
        BigFloat::try_from_i64_exact(n, p).unwrap()
    }

    #[test]
    fn exact_ball_has_max_accuracy() {
        let b = Ball::point(bf(5, 53)).unwrap();
        assert_eq!(b.rel_accuracy_bits(), i64::MAX);
        assert_eq!(b.abs_error(), Mag::ZERO);
        assert_eq!(b.rel_error_bound(), Mag::ZERO);
    }

    #[test]
    fn entire_ball_has_min_accuracy() {
        let b = Ball::new(bf(1, 53), Mag::INFINITY).unwrap();
        assert_eq!(b.rel_accuracy_bits(), i64::MIN);
        assert_eq!(b.rel_error_bound(), Mag::INFINITY);
    }

    #[test]
    fn rel_accuracy_tracks_radius_position() {
        // mid = 1 (exponent 0), rad = 2^-50: ~50 accurate bits.
        let b = Ball::new(bf(1, 53), Mag::from_pow2(-50)).unwrap();
        assert_eq!(b.rel_accuracy_bits(), 50);
        // mid = 1024 (exponent 10), rad = 2^-50: ~60 accurate bits.
        let b2 = Ball::new(bf(1024, 53), Mag::from_pow2(-50)).unwrap();
        assert_eq!(b2.rel_accuracy_bits(), 60);
    }

    #[test]
    fn zero_midpoint_positive_radius_has_no_relative_accuracy() {
        let b = Ball::new(bf(0, 53), Mag::from_pow2(-10)).unwrap();
        assert_eq!(b.rel_accuracy_bits(), i64::MIN);
        assert_eq!(b.rel_error_bound(), Mag::INFINITY);
    }

    #[test]
    fn rel_error_bound_is_sound_upper_bound() {
        // [4 ± 1]: relative error 1/4 = 0.25. The Mag bound must be ≥ 0.25.
        let b = Ball::new(bf(4, 53), Mag::from_pow2(0)).unwrap();
        let rel = b.rel_error_bound();
        let quarter = Mag::from_pow2(-2); // 0.25
        assert!(rel >= quarter, "rel_error_bound must be ≥ true 0.25");
    }

    #[test]
    fn refine_drives_precision_until_accurate() {
        // Compute 1/3 at precision p: the radius shrinks ~ 2^-p, so the
        // accuracy grows with p. Ask for 100 bits.
        let (ball, _) = refine_to_accuracy(100, 16, 4096, |p| {
            let one = Ball::point(bf(1, p)).unwrap();
            let three = Ball::point(bf(3, p)).unwrap();
            one.div(&three)
        });
        assert!(ball.rel_accuracy_bits() >= 100);
        // And it actually encloses the true 1/3.
        let third = bf(1, 400).div(&bf(3, 400), RoundingMode::NearestEven).0;
        assert!(ball.lower().partial_cmp(&third).0 != Some(Ordering::Greater));
        assert!(ball.upper().partial_cmp(&third).0 != Some(Ordering::Less));
    }

    #[test]
    fn refine_caps_at_max_precision_for_unreachable_target() {
        // Entire result can never reach the target; the driver returns at
        // max_precision without looping forever.
        let (ball, _) = refine_to_accuracy(1000, 16, 256, |p| {
            let one = Ball::point(bf(1, p)).unwrap();
            let zero = Ball::point(bf(0, p)).unwrap();
            one.div(&zero) // entire + DIV_BY_ZERO
        });
        assert!(ball.is_entire());
    }

    #[test]
    fn refine_short_circuits_on_exact() {
        // An exact computation returns at the start precision.
        let mut calls = 0;
        let (ball, _) = refine_to_accuracy(1000, 16, 4096, |p| {
            calls += 1;
            Ball::point(bf(5, p))
                .unwrap()
                .add(&Ball::point(bf(0, p)).unwrap())
        });
        assert!(ball.is_exact());
        assert_eq!(calls, 1, "exact ball should not trigger a retry");
    }
}
