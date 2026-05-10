//! `log2(x)`: base-2 logarithm.
//!
//! Composition: `log2(x) = ln(x) / ln(2)` at working precision
//! `target + 64`. Special cases flow through composition: `ln(x)`
//! handles the entire IEEE 754-2019 §9.2 dispatch (NaN, `±0` →
//! `−∞ + DIV_BY_ZERO`, negative finite → `qNaN + INVALID`, `±∞`),
//! and the subsequent division by the positive `ln(2)` preserves
//! sign, infinity, and NaN.
//!
//! pfloat stores normal `BigFloat` as `m × 2^e` with the top bit of
//! the mantissa set. An exact-power-of-two fast path would return
//! `e` directly without rounding; slice 3d skips that optimization
//! and lets the composition produce the same value via `ln`. The
//! result on a power-of-two input rounds to the exact integer at
//! the 64-bit guard precision.

use crate::big::{BigFloat, BuildError};
use crate::rounding::RoundingMode;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::ln_2_at;

impl BigFloat {
    /// `log2(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn log2(&self, mode: RoundingMode) -> (Self, Status) {
        self.log2_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `log2(self)` with explicit result precision.
    pub fn log2_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(log2_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `log2(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::log2`].
    #[must_use]
    pub fn log2(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().log2(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn log2_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    let working_prec = target_precision.saturating_add(64).min(1024);
    let (ln_x, ln_status) = x
        .ln_round(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1");
    let ln_2 = ln_2_at(working_prec);
    let (ratio, div_status) = ln_x.div(&ln_2, RoundingMode::NearestEven);
    let (rounded, round_status) = ratio
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    auto_raise(round_status);
    (rounded, ln_status | div_status | round_status)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sign::Sign;
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
    fn log2_one_is_zero() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.log2(RoundingMode::NearestEven);
        assert!(r.is_zero());
    }

    #[test]
    fn log2_two_is_one() {
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (r, _) = two.log2(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        assert!(close_at(&r, &one, 113 - 12));
    }

    #[test]
    fn log2_eight_is_three() {
        let eight = BigFloat::try_from_i64_exact(8, 113).unwrap();
        let (r, _) = eight.log2(RoundingMode::NearestEven);
        let three = BigFloat::try_from_i64_exact(3, 113).unwrap();
        assert!(close_at(&r, &three, 113 - 12));
    }

    #[test]
    fn log2_zero_is_neg_inf() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, status) = z.log2(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative());
        assert!(status.div_by_zero());
    }

    #[test]
    fn log2_negative_is_invalid() {
        let neg = BigFloat::try_from_i64_exact(-5, 53).unwrap();
        let (r, status) = neg.log2(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn log2_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.log2(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn log2_neg_inf_is_invalid() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, status) = ni.log2(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn log2_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.log2(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn log2_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.log2(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn log2_round_trip_with_exp2() {
        for n in &[2i64, 3, 5, 7, 10, 100] {
            let x = BigFloat::try_from_i64_exact(*n, 113).unwrap();
            let (l, _) = x.log2(RoundingMode::NearestEven);
            let (back, _) = l.exp2(RoundingMode::NearestEven);
            assert!(close_at(&back, &x, 113 - 12), "log2/exp2 round-trip at {n}");
        }
    }
}
