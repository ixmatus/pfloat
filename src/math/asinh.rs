//! `asinh(x) = ln(x + sqrt(x² + 1))`: inverse hyperbolic sine.
//!
//! Identity used: `asinh(x) = sign(x) · log1p(|x| + |x|² /
//! (sqrt(|x|² + 1) + 1))`. The formula has no cancellation for
//! `|x| ≥ 0`: every term in the `log1p` argument is non-negative,
//! and `sqrt(|x|² + 1) + 1 ≥ 2 > 0` for any `x`. Naively evaluating
//! `ln(x + sqrt(x² + 1))` for large negative `x` cancels leading
//! bits because `x + sqrt(x² + 1) → 0+`; computing on `|x|` and
//! applying the sign avoids the issue.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `asinh(±0) = ±0`.
//! - `asinh(±∞) = ±∞`.
//! - `asinh(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::ziv::ziv_round_with_depth;
use super::ziv_calibration::ASINH_ERROR_GUARD;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `asinh(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn asinh(&self, mode: RoundingMode) -> (Self, Status) {
        self.asinh_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `asinh(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.27, ADR-0038).
    pub fn asinh_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(asinh_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `asinh(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::asinh`].
    #[must_use]
    pub fn asinh(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().asinh(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn asinh_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    match &x.class {
        Class::Nan {
            quiet,
            sign,
            payload,
        } => {
            if !*quiet {
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            let nan = BigFloat::try_new_quiet_nan(*sign, target_precision, payload)
                .expect("precision >= 1");
            return (nan, Status::OK);
        }
        Class::Zero { sign } => {
            let z = BigFloat::try_new_zero(*sign, target_precision).expect("precision >= 1");
            return (z, Status::OK);
        }
        Class::Infinity { sign } => {
            return (
                BigFloat::try_new_infinity(*sign, target_precision).expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Normal { .. } => {}
    }

    // Tiny x: asinh(x) = x − x³/6 + … shrinks toward x in magnitude
    // (the −x³/6 correction opposes x's sign), so round x with that
    // opposite-sign infinitesimal directly, bypassing the log1p
    // composition that otherwise grinds at moderate working precision
    // for moderately tiny x (ADR-0059).
    let e = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => unreachable!(),
    };
    // The depth must also clear the INPUT's grid (pf-fbjn, ADR-0104):
    // a high-precision x parked next to a rounding-change point puts
    // the cubic series correction (position 3e) across a boundary the
    // residue (position e − p − 2) never reaches. Arm-failing inputs
    // go to the driver, whose ADR-0103 deep rung takes the input at
    // full precision and certifies the true boundary side.
    if e <= -(i64::from(target_precision) + 2)
        && e.saturating_mul(-2) >= i64::from(x.precision).saturating_add(6)
    {
        return crate::rounding::round_with_infinitesimal(
            x,
            x.sign(),
            true, // magnitude shrinks: the −x³/6 correction opposes x's sign
            target_precision,
            mode,
        );
    }

    // Ziv-driven correct rounding under every IEEE mode. The eval
    // closure runs the cancellation-resistant identity
    // `asinh(|x|) = log1p(|x| + |x|²/(sqrt(|x|²+1)+1))` on |x| at
    // working precision `w`, with sign reapplied. The composition
    // has no cancellation regime (every term in the log1p argument
    // is non-negative, divisor ≥ 2). log1p is Ziv-driven (slice
    // p1.24).
    let sign = x.sign();
    let abs_x = x.abs();
    let (result, status) = ziv_round_with_depth(
        |w| {
            let x_w = abs_x
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let one = BigFloat::try_from_i64_exact(1, w).expect("precision >= 1");
            let (x_sq, _) = x_w.mul(&x_w, RoundingMode::NearestEven);
            let (x_sq_plus_one, _) = x_sq.add(&one, RoundingMode::NearestEven);
            let (s, _) = x_sq_plus_one.sqrt(RoundingMode::NearestEven);
            let (s_plus_one, _) = s.add(&one, RoundingMode::NearestEven);
            let (correction, _) = x_sq.div(&s_plus_one, RoundingMode::NearestEven);
            let (arg, _) = x_w.add(&correction, RoundingMode::NearestEven);
            let (lp, _) = arg.log1p(RoundingMode::NearestEven);
            if matches!(sign, Sign::Negative) {
                lp.negated()
            } else {
                lp
            }
        },
        target_precision,
        mode,
        ASINH_ERROR_GUARD,
        // Parked-input certification depth (pf-fbjn, ADR-0104):
        // arm-rejected tiny inputs resolve at the deep rung, which
        // must reach both the input's precision and the series
        // correction's depth. Lazy: free unless the schedule exhausts.
        || {
            if e < 0 {
                u32::try_from(e.saturating_mul(-5))
                    .unwrap_or(u32::MAX)
                    .max(x.precision)
                    .saturating_add(64)
            } else {
                0
            }
        },
    );
    // asinh(x) for finite normal x ≠ 0 is transcendental (Lindemann–
    // Weierstrass), hence irrational, hence INEXACT even where it rounds
    // onto a grid value (pf-uqd1, ADR-0063). asinh(±0) = ±0 is dispatched
    // above; the tiny-x fast path sets INEXACT via round_with_infinitesimal.
    let status = super::force_transcendental_inexact(&result, status);
    auto_raise(status);
    (result, status)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    fn close_at(v: &BigFloat, expected: &BigFloat, bits: u32) -> bool {
        let (diff, _) = v.sub(expected, RoundingMode::NearestEven);
        let abs_diff = diff.abs();
        if abs_diff.is_zero() {
            return true;
        }
        let p = v.precision().max(expected.precision());
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let abs_b = expected.abs();
        let mut bound = if abs_b.is_zero() { one } else { abs_b };
        for _ in 0..bits {
            bound = bound.div(&two, RoundingMode::NearestEven).0;
        }
        matches!(
            abs_diff.partial_cmp(&bound).0,
            Some(Ordering::Less | Ordering::Equal)
        )
    }

    #[test]
    fn asinh_zero() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.asinh(RoundingMode::NearestEven);
            assert!(r.is_zero());
            assert_eq!(r.is_sign_negative(), matches!(s, Sign::Negative));
        }
    }

    #[test]
    fn asinh_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.asinh(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn asinh_neg_inf() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, _) = ni.asinh(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative());
    }

    #[test]
    fn asinh_sinh_round_trip() {
        // asinh(sinh(x)) = x for moderate x.
        let p = 113u32;
        for n in &[-3i64, -1, 1, 3, 5] {
            let x = BigFloat::try_from_i64_exact(*n, p).unwrap();
            let (s, _) = x.sinh(RoundingMode::NearestEven);
            let (back, _) = s.asinh(RoundingMode::NearestEven);
            assert!(
                close_at(&back, &x, p - 16),
                "asinh(sinh({n})) = {back}, expected {x}"
            );
        }
    }

    #[test]
    fn asinh_one() {
        // asinh(1) = ln(1 + sqrt(2)) ≈ 0.8813735870195429
        let p = 113u32;
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (r, _) = one.asinh(RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let (sqrt2, _) = two.sqrt(RoundingMode::NearestEven);
        let (arg, _) = one.add(&sqrt2, RoundingMode::NearestEven);
        let (expected, _) = arg.ln(RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, p - 12));
    }

    #[test]
    fn asinh_negation() {
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let neg_two = BigFloat::try_from_i64_exact(-2, 113).unwrap();
        let (a, _) = two.asinh(RoundingMode::NearestEven);
        let (b, _) = neg_two.asinh(RoundingMode::NearestEven);
        let neg_a = a.negated();
        assert!(close_at(&b, &neg_a, 113 - 8));
    }

    #[test]
    fn asinh_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.asinh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn asinh_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.asinh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn asinh_tiny_input_directed_modes() {
        // Tiny-x short-circuit (ADR-0059). asinh(x) = x − x³/6 + …
        // shrinks toward zero, so for tiny x the result is x under the
        // three modes that round toward or onto x and x's toward-zero
        // neighbour under the two that round inward. Pins the shrink
        // direction (the signature of subtracts_magnitude = true): a
        // flipped flag or a bare "return x" would round the inward modes
        // to x.
        use core::cmp::Ordering;
        use RoundingMode::{NearestAway, NearestEven, TowardNegative, TowardPositive, TowardZero};
        let two = BigFloat::try_from_i64_exact(2, 200).unwrap();
        let mut x = BigFloat::try_from_i64_exact(1, 200).unwrap();
        for _ in 0..100 {
            x = x.div(&two, RoundingMode::NearestEven).0;
        }
        let x53 = x
            .round_to_precision(53, RoundingMode::NearestEven)
            .unwrap()
            .0;

        // Positive tiny x: TowardZero and TowardNegative round below x.
        for &m in &[NearestEven, NearestAway, TowardPositive] {
            let (r, _) = x53.asinh(m);
            assert_eq!(
                r.partial_cmp(&x53).0,
                Some(Ordering::Equal),
                "asinh +tiny {m:?}"
            );
        }
        for &m in &[TowardZero, TowardNegative] {
            let (r, _) = x53.asinh(m);
            assert_eq!(
                r.partial_cmp(&x53).0,
                Some(Ordering::Less),
                "asinh +tiny {m:?} shrinks"
            );
        }

        // Negative tiny x: TowardZero and TowardPositive round above x.
        let neg = x53.negated();
        for &m in &[NearestEven, NearestAway, TowardNegative] {
            let (r, _) = neg.asinh(m);
            assert_eq!(
                r.partial_cmp(&neg).0,
                Some(Ordering::Equal),
                "asinh -tiny {m:?}"
            );
        }
        for &m in &[TowardZero, TowardPositive] {
            let (r, _) = neg.asinh(m);
            assert_eq!(
                r.partial_cmp(&neg).0,
                Some(Ordering::Greater),
                "asinh -tiny {m:?} shrinks"
            );
        }
    }
}
