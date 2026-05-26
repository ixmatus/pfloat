//! `log1p(x) = ln(1 + x)`: natural logarithm of one plus the
//! argument.
//!
//! Naïvely `ln(1 + x)` loses precision near zero: for `x ≈ 2^−n`,
//! the addition `1 + x` cancels the leading bits of `x`. The kernel
//! boosts working precision by `−exponent(x)` bits, computes
//! `1 + x` and `ln` at that precision, then rounds back.
//!
//! For `|x|` so small that `log1p(x)` rounds to `x` at target
//! precision (`x.exponent ≤ −target − 8`), the kernel short-circuits
//! and returns `x` directly.
//!
//! Correctly rounded under every IEEE 754-2019 rounding mode via
//! the shared [`crate::math::ziv::ziv_round`] driver (slice p1.24,
//! ADR-0038). The cancellation-boost composition runs inside the
//! eval closure at each Ziv working precision; the outer envelope
//! certifies the rounding-mode interval test on the final round.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `log1p(±0) = ±0`.
//! - `log1p(−1) = −∞ + DIV_BY_ZERO` (§7.3).
//! - `log1p(x) = qNaN + INVALID` for `x < −1`.
//! - `log1p(+∞) = +∞`, `log1p(−∞) = qNaN + INVALID`.
//! - `log1p(NaN) = NaN`; `sNaN` raises `INVALID`.

use core::cmp::Ordering;

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::ziv::ziv_round;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `ln(1 + self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn log1p(&self, mode: RoundingMode) -> (Self, Status) {
        self.log1p_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `log1p(self)` with explicit result precision.
    pub fn log1p_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(log1p_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `log1p(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::log1p`].
    #[must_use]
    pub fn log1p(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().log1p(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn log1p_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
        Class::Infinity {
            sign: Sign::Positive,
        } => {
            return (
                BigFloat::try_new_infinity(Sign::Positive, target_precision)
                    .expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Infinity {
            sign: Sign::Negative,
        } => {
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        Class::Normal { .. } => {}
    }

    // x is Normal. Check x ≤ −1.
    let neg_one = BigFloat::try_from_i64_exact(-1, x.precision).expect("precision >= 1");
    match x.partial_cmp(&neg_one).0 {
        Some(Ordering::Equal) => {
            let ninf = BigFloat::try_new_infinity(Sign::Negative, target_precision)
                .expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            return (ninf, Status::DIV_BY_ZERO);
        }
        Some(Ordering::Less) => {
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        _ => {}
    }

    let e = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => unreachable!(),
    };

    // Cancellation boost: `1 + x` cancels ~|exponent(x)| leading
    // bits when x is small. The boost moves INSIDE the Ziv eval
    // closure so each working-precision retry inherits it.
    //
    // The pre-Phase-1f kernel had a short-circuit at
    // `e ≤ -target - 8` that returned `x rounded under mode`.
    // That shortcut was NE-correct but produced wrong
    // directed-mode results (same shape as expm1 — for positive
    // small x under TP, true log1p(x) < x but rounds up under TP
    // to x's neighbour ABOVE; for negative small x under TZ, true
    // log1p(x) > x and rounds toward zero away from x). The Ziv
    // driver converges in 1-2 iterations on tiny-x inputs and
    // certifies the correct rounding-mode behaviour.
    let cancellation: u32 = if e < 0 {
        u32::try_from(-e).unwrap_or(u32::MAX)
    } else {
        0
    };

    // Ziv-driven correct rounding under every IEEE mode. The eval
    // closure runs the existing composition `ln(1 + x_w)` at a
    // working precision boosted by `cancellation` above the Ziv
    // driver's requested working precision `w`, then rounds the
    // composition's result to `w` under NE so the Ziv interval
    // test sees a w-precision value with the cancellation absorbed.
    let (result, status) = ziv_round(
        |w| {
            let inner_w = w.saturating_add(cancellation).min(w.saturating_add(1024));
            let x_w = x
                .round_to_precision(inner_w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let one = BigFloat::try_from_i64_exact(1, inner_w).expect("precision >= 1");
            let (one_plus_x, _) = one.add(&x_w, RoundingMode::NearestEven);
            let (ln_val, _) = one_plus_x.ln(RoundingMode::NearestEven);
            ln_val
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0
        },
        target_precision,
        mode,
    );
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
    fn log1p_zero_is_zero() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.log1p(RoundingMode::NearestEven);
            assert!(r.is_zero());
            assert_eq!(r.is_sign_negative(), matches!(s, Sign::Negative));
        }
    }

    #[test]
    fn log1p_neg_one_is_neg_inf() {
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        let (r, status) = neg_one.log1p(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative());
        assert!(status.div_by_zero());
    }

    #[test]
    fn log1p_below_neg_one_is_invalid() {
        let neg_two = BigFloat::try_from_i64_exact(-2, 53).unwrap();
        let (r, status) = neg_two.log1p(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn log1p_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.log1p(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn log1p_neg_inf_is_invalid() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, status) = ni.log1p(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn log1p_at_one_is_ln_two() {
        // log1p(1) = ln(2)
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.log1p(RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (expected, _) = two.ln(RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, 113 - 12));
    }

    #[test]
    fn log1p_small_x_round_trip() {
        // log1p(small) ≈ small. Round-trip via expm1: expm1(log1p(x)) ≈ x.
        let p = 113u32;
        let x = BigFloat::parse_str("1e-30", p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (lp, _) = x.log1p(RoundingMode::NearestEven);
        let (back, _) = lp.expm1(RoundingMode::NearestEven);
        assert!(close_at(&back, &x, p - 12));
    }

    #[test]
    fn log1p_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.log1p(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn log1p_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.log1p(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
