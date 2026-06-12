//! `atanh(x) = (1/2) · ln((1 + x)/(1 − x))`: inverse hyperbolic
//! tangent, defined for `|x| < 1`.
//!
//! Identity used: `atanh(x) = (log1p(x) − log1p(−x)) / 2`. Both
//! `log1p` calls handle their small-argument cancellation regimes
//! internally; for `x → 0` the subtraction `log1p(x) − log1p(−x)`
//! reduces to `2x` in real math without leading-bit cancellation
//! in the `BigFloat` representation.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `atanh(±0) = ±0`.
//! - `atanh(±1) = ±∞ + DIV_BY_ZERO`.
//! - `atanh(x) = qNaN + INVALID` for `|x| > 1`, including `±∞`.
//! - `atanh(NaN) = NaN`; `sNaN` raises `INVALID`.

use core::cmp::Ordering;

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::ziv::ziv_round_with_depth;
use super::ziv_calibration::ATANH_ERROR_GUARD;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `atanh(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn atanh(&self, mode: RoundingMode) -> (Self, Status) {
        self.atanh_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `atanh(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.27, ADR-0038).
    pub fn atanh_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(atanh_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `atanh(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::atanh`].
    #[must_use]
    pub fn atanh(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().atanh(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn atanh_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
        Class::Infinity { .. } => {
            // |∞| > 1: domain error.
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        Class::Normal { .. } => {}
    }

    // Domain dispatch on |x| vs 1.
    let sign = x.sign();
    let abs_x = x.abs();
    let one_at_input = BigFloat::try_from_i64_exact(1, x.precision).expect("precision >= 1");
    match abs_x.partial_cmp(&one_at_input).0 {
        Some(Ordering::Equal) => {
            // atanh(±1) = ±∞ + DIV_BY_ZERO.
            let inf = BigFloat::try_new_infinity(sign, target_precision).expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            return (inf, Status::DIV_BY_ZERO);
        }
        Some(Ordering::Greater) => {
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        _ => {}
    }

    // Tiny x: atanh(x) = x + x³/3 + … grows away from x in magnitude
    // (every term of the odd series shares x's sign), so round x with
    // that same-sign infinitesimal directly, bypassing the Ziv loop that
    // otherwise drives the full log1p(x) − log1p(−x) identity at high
    // working precision for a value that is x to within rounding. The
    // infinitesimal carries the directed-mode information a bare
    // "x rounded under mode" return would drop — the same trap log1p
    // documents at its own tiny-x short-circuit (ADR-0059).
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
            sign,
            false, // magnitude grows: the +x³/3 correction shares x's sign
            target_precision,
            mode,
        );
    }

    // Ziv-driven correct rounding under every IEEE mode. The eval
    // closure runs the cancellation-resistant identity
    // `atanh(x) = (log1p(x) − log1p(-x)) / 2` at working precision
    // `w`. log1p is Ziv-driven (slice p1.24) and each log1p
    // handles its own small-argument cancellation; the outer
    // subtraction adds the leading terms (no cancellation).
    let (result, status) = ziv_round_with_depth(
        |w| {
            let x_w = x
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let neg_x = x_w.negated();
            let (lp_pos, _) = x_w.log1p(RoundingMode::NearestEven);
            let (lp_neg, _) = neg_x.log1p(RoundingMode::NearestEven);
            let (diff, _) = lp_pos.sub(&lp_neg, RoundingMode::NearestEven);
            let two = BigFloat::try_from_i64_exact(2, w).expect("precision >= 1");
            diff.div(&two, RoundingMode::NearestEven).0
        },
        target_precision,
        mode,
        ATANH_ERROR_GUARD,
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
    // atanh(x) for finite normal x with 0 < |x| < 1 is transcendental
    // (Lindemann–Weierstrass), hence irrational, hence INEXACT even where
    // it rounds onto a grid value (pf-uqd1, ADR-0063). atanh(±0) = ±0 is
    // dispatched above, the poles at ±1 too; the tiny-x fast path sets
    // INEXACT via round_with_infinitesimal.
    let status = if matches!(result.class, Class::Normal { .. }) {
        status | Status::INEXACT
    } else {
        status
    };
    auto_raise(status);
    (result, status)
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn atanh_zero() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.atanh(RoundingMode::NearestEven);
            assert!(r.is_zero());
            assert_eq!(r.is_sign_negative(), matches!(s, Sign::Negative));
        }
    }

    #[test]
    fn atanh_one_is_pos_inf_div() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, status) = one.atanh(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
        assert!(status.div_by_zero());
    }

    #[test]
    fn atanh_neg_one_is_neg_inf_div() {
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        let (r, status) = neg_one.atanh(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative());
        assert!(status.div_by_zero());
    }

    #[test]
    fn atanh_above_one_is_invalid() {
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, status) = two.atanh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn atanh_pos_inf_is_invalid() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, status) = pi.atanh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn atanh_half() {
        // atanh(0.5) = (1/2)·ln(3) ≈ 0.5493061443340549
        let p = 113u32;
        let half = BigFloat::parse_str("0.5", p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = half.atanh(RoundingMode::NearestEven);
        let three = BigFloat::try_from_i64_exact(3, p).unwrap();
        let (ln3, _) = three.ln(RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let (expected, _) = ln3.div(&two, RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, p - 12));
    }

    #[test]
    fn atanh_tanh_round_trip() {
        // atanh(tanh(x)) = x for moderate x.
        let p = 113u32;
        for n in &[-3i64, -1, 1, 2, 4] {
            let x = BigFloat::try_from_i64_exact(*n, p).unwrap();
            let (t, _) = x.tanh(RoundingMode::NearestEven);
            let (back, _) = t.atanh(RoundingMode::NearestEven);
            assert!(
                close_at(&back, &x, p - 16),
                "atanh(tanh({n})) = {back}, expected {x}"
            );
        }
    }

    #[test]
    fn atanh_negation() {
        let p = 113u32;
        let half = BigFloat::parse_str("0.5", p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let neg_half = half.negated();
        let (a, _) = half.atanh(RoundingMode::NearestEven);
        let (b, _) = neg_half.atanh(RoundingMode::NearestEven);
        let neg_a = a.negated();
        assert!(close_at(&b, &neg_a, p - 12));
    }

    #[test]
    fn atanh_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.atanh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn atanh_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.atanh(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn atanh_tiny_input_directed_modes() {
        // Tiny-x short-circuit (ADR-0059). atanh(x) = x + x³/3 + …
        // grows away from x: for x = 2⁻¹⁰⁰ at p=53 the cubic correction
        // ≈ 2⁻³⁰⁰ sits far below the half-ULP ≈ 2⁻¹⁵³, so the result is
        // x under the four modes that round toward x and x's
        // away-from-zero neighbour under the single mode that rounds
        // away. This is the visible signature of `subtracts_magnitude =
        // false`; a regression to a bare "return x" or a flipped flag
        // would round the away mode to x, which this asserts against.
        // Unlike tanh (which shrinks and defers directed modes to the
        // sweep), the grow direction is known, so all five modes are
        // pinned here.
        let two = BigFloat::try_from_i64_exact(2, 200).unwrap();
        let mut x = BigFloat::try_from_i64_exact(1, 200).unwrap();
        for _ in 0..100 {
            x = x.div(&two, RoundingMode::NearestEven).0;
        }
        let x53 = x
            .round_to_precision(53, RoundingMode::NearestEven)
            .unwrap()
            .0;

        // x's away-from-zero neighbour at p=53: x + 2⁻¹⁵². The two set
        // bits span 52 places, so the add is exact at p=53.
        let mut ulp = BigFloat::try_from_i64_exact(1, 53).unwrap();
        for _ in 0..152 {
            ulp = ulp.div(&two, RoundingMode::NearestEven).0;
        }
        let (x_up, _) = x53.add(&ulp, RoundingMode::NearestEven);

        use RoundingMode::{NearestAway, NearestEven, TowardNegative, TowardPositive, TowardZero};

        // Positive tiny x: only TowardPositive rounds away (up).
        for &mode in &[NearestEven, NearestAway, TowardZero, TowardNegative] {
            let (r, _) = x53.atanh(mode);
            assert_eq!(
                r.partial_cmp(&x53).0,
                Some(Ordering::Equal),
                "atanh(2^-100) under {mode:?} = {r}, expected {x53}"
            );
        }
        let (r_tp, _) = x53.atanh(TowardPositive);
        assert_eq!(
            r_tp.partial_cmp(&x_up).0,
            Some(Ordering::Equal),
            "atanh(2^-100) under TowardPositive = {r_tp}, expected {x_up}"
        );

        // Negative tiny x: mirror, only TowardNegative rounds away (down).
        let neg = x53.negated();
        let neg_down = x_up.negated();
        for &mode in &[NearestEven, NearestAway, TowardZero, TowardPositive] {
            let (r, _) = neg.atanh(mode);
            assert_eq!(
                r.partial_cmp(&neg).0,
                Some(Ordering::Equal),
                "atanh(-2^-100) under {mode:?} = {r}, expected {neg}"
            );
        }
        let (r_tn, _) = neg.atanh(TowardNegative);
        assert_eq!(
            r_tn.partial_cmp(&neg_down).0,
            Some(Ordering::Equal),
            "atanh(-2^-100) under TowardNegative = {r_tn}, expected {neg_down}"
        );
    }
}
