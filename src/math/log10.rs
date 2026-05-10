//! `log10(x)`: base-10 logarithm.
//!
//! Composition: `log10(x) = ln(x) / ln(10)` at working precision
//! `target + 64`. Special cases compose through `ln` and the
//! division (see [`super::log2`] for the full reasoning).
//!
//! Powers of 10 are not exact in binary floating-point, so unlike
//! `log2` there is no fast path for `log10(10^k) = k`. The
//! composition produces a target-precision integer up to one ULP
//! away from `k`, which rounds to the exact integer at the 64-bit
//! guard precision for moderate `k`.

use crate::big::{BigFloat, BuildError};
use crate::rounding::RoundingMode;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::ln_10_at;

impl BigFloat {
    /// `log10(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn log10(&self, mode: RoundingMode) -> (Self, Status) {
        self.log10_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `log10(self)` with explicit result precision.
    pub fn log10_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(log10_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `log10(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::log10`].
    #[must_use]
    pub fn log10(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().log10(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn log10_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    let working_prec = target_precision.saturating_add(64).min(1024);
    let (ln_x, ln_status) = x
        .ln_round(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1");
    let ln_10 = ln_10_at(working_prec);
    let (ratio, div_status) = ln_x.div(&ln_10, RoundingMode::NearestEven);
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
    fn log10_one_is_zero() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.log10(RoundingMode::NearestEven);
        assert!(r.is_zero());
    }

    #[test]
    fn log10_ten_is_one() {
        let ten = BigFloat::try_from_i64_exact(10, 113).unwrap();
        let (r, _) = ten.log10(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        assert!(close_at(&r, &one, 113 - 12));
    }

    #[test]
    fn log10_thousand_is_three() {
        let k = BigFloat::try_from_i64_exact(1000, 113).unwrap();
        let (r, _) = k.log10(RoundingMode::NearestEven);
        let three = BigFloat::try_from_i64_exact(3, 113).unwrap();
        assert!(close_at(&r, &three, 113 - 12));
    }

    #[test]
    fn log10_zero_is_neg_inf() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, status) = z.log10(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative());
        assert!(status.div_by_zero());
    }

    #[test]
    fn log10_negative_is_invalid() {
        let neg = BigFloat::try_from_i64_exact(-5, 53).unwrap();
        let (r, status) = neg.log10(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn log10_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.log10(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn log10_neg_inf_is_invalid() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, status) = ni.log10(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn log10_round_trip_with_exp10() {
        for n in &[2i64, 5, 10, 50, 100] {
            let x = BigFloat::try_from_i64_exact(*n, 113).unwrap();
            let (l, _) = x.log10(RoundingMode::NearestEven);
            let (back, _) = l.exp10(RoundingMode::NearestEven);
            assert!(
                close_at(&back, &x, 113 - 12),
                "log10/exp10 round-trip at {n}"
            );
        }
    }

    #[test]
    fn log10_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.log10(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn log10_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.log10(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
