//! `exp10(x) = 10^x`: decimal exponential.
//!
//! Composition: `10^x = exp(x · ln(10))`. The kernel computes the
//! product at working precision against `ln(10)` (from
//! [`super::ln_10_at`], which evaluates `ln(10)` lazily per call),
//! then dispatches to `exp`. Special cases compose through the same
//! pattern as `exp2`.
//!
//! Correctly rounded under every IEEE 754-2019 rounding mode via
//! the shared [`crate::math::ziv::ziv_round`] driver (slice p1.24,
//! ADR-0038). The `exp · ln(10)` composition has no cancellation
//! regime; the Ziv envelope's working-precision growth certifies
//! the rounding-mode interval test on the final round.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `exp10(±0) = 1`.
//! - `exp10(+∞) = +∞`, `exp10(−∞) = +0`.
//! - `exp10(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::ln_10_at;
use super::ziv::ziv_round;
use super::ziv_calibration::EXP10_ERROR_GUARD;

impl BigFloat {
    /// `10^self` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn exp10(&self, mode: RoundingMode) -> (Self, Status) {
        self.exp10_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `exp10(self)` with explicit result precision.
    pub fn exp10_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(exp10_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `exp10(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::exp10`].
    #[must_use]
    pub fn exp10(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().exp10(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn exp10_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
        Class::Zero { .. } => {
            return (
                BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1"),
                Status::OK,
            );
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
            return (
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Normal { .. } => {}
    }

    // Exact-input dispatch (pf-njs5, ADR-0060). 10^x is exactly
    // representable iff x is a non-negative integer and 10^x fits the
    // target precision (10^x = 5^x·2^x; the odd 5^x factor sets the
    // significant-bit count, the 2^x factor is a free exponent
    // shift). For non-integer x the value is irrational (Lindemann–
    // Weierstrass via 10^x = exp(x·ln 10)); for negative integer x it
    // is 1/10^|x|, not dyadic. Both fall through and force INEXACT.
    if let Some(v) = exp10_exact_if_small_nonneg_int(x, target_precision) {
        return (v, Status::OK);
    }

    // Ziv-driven correct rounding under every IEEE mode. The
    // composition `exp(x · ln(10))` has no cancellation regime; the
    // Ziv driver's working-precision growth handles the rounding-
    // mode interval test at the final round to target.
    let (result, status) = ziv_round(
        |w| {
            let x_w = x
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let ln_10 = ln_10_at(w);
            let (product, _) = x_w.mul(&ln_10, RoundingMode::NearestEven);
            let (e_val, _) = product.exp(RoundingMode::NearestEven);
            e_val
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0
        },
        target_precision,
        mode,
        EXP10_ERROR_GUARD,
    );
    // x is not in the exact-input set ⟹ 10^x is irrational or its
    // exact value does not fit ⟹ INEXACT, even where the working-
    // precision evaluation rounds onto a grid value.
    let status = status | Status::INEXACT;
    auto_raise(status);
    (result, status)
}

/// `Some(10^x)` at `target_precision` if `x` is a non-negative
/// integer whose exact `10^x` fits the target precision, else `None`.
/// The dispatch returns this value with `Status::OK`; soundness (no
/// wrongly cleared flag) rests on [`super::pow::ten_pow_if_fits`]
/// returning only exact, representable powers of ten.
fn exp10_exact_if_small_nonneg_int(x: &BigFloat, target_precision: u32) -> Option<BigFloat> {
    let k = super::pow::integer_exponent(x)?;
    if k < 0 {
        return None;
    }
    super::pow::ten_pow_if_fits(k as u64, target_precision)
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
    fn exp10_zero_is_one() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, _) = z.exp10(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    #[test]
    fn exp10_one_is_ten() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.exp10(RoundingMode::NearestEven);
        let ten = BigFloat::try_from_i64_exact(10, 113).unwrap();
        assert!(close_at(&r, &ten, 113 - 12));
    }

    #[test]
    fn exp10_three_is_thousand() {
        let three = BigFloat::try_from_i64_exact(3, 113).unwrap();
        let (r, _) = three.exp10(RoundingMode::NearestEven);
        let k = BigFloat::try_from_i64_exact(1000, 113).unwrap();
        assert!(close_at(&r, &k, 113 - 12));
    }

    #[test]
    fn exp10_negative() {
        // 10^-2 = 0.01
        let neg_two = BigFloat::try_from_i64_exact(-2, 113).unwrap();
        let (r, _) = neg_two.exp10(RoundingMode::NearestEven);
        let hundredth = BigFloat::parse_str("0.01", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        assert!(close_at(&r, &hundredth, 113 - 12));
    }

    #[test]
    fn exp10_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.exp10(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn exp10_neg_inf_is_zero() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, _) = ni.exp10(RoundingMode::NearestEven);
        assert!(r.is_zero() && r.is_sign_positive());
    }

    #[test]
    fn exp10_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.exp10(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn exp10_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.exp10(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
