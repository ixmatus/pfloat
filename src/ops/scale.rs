//! Exact power-of-two scaling for [`BigFloat`] (IEEE 754-2019 §5.3.3
//! `scaleB`, restricted to an exact `2^k` scale).
//!
//! Multiplying a binary float by `2^k` is exact: it shifts the
//! unbiased exponent by `k` and leaves the mantissa, sign, and
//! precision untouched. pfloat already mutates the `Normal` exponent
//! this way inside the arithmetic kernels; this module surfaces the
//! operation through the validated public API without exposing a
//! raw-parts constructor. ADR-0072.
//!
//! Special cases per IEEE 754-2019 §6.2 and §6.3:
//!
//! - NaN: a quiet NaN propagates (sign and payload preserved); a
//!   signaling NaN raises [`Status::INVALID`] and is quieted, matching
//!   every other general-computational kernel (`scaleB` is
//!   general-computational per §5.3.3).
//! - `±0` and `±∞`: returned unchanged (`±0 × 2^k = ±0`,
//!   `±∞ × 2^k = ±∞`); the sign is preserved and no flag is raised.
//! - Finite non-zero: the exponent shifts by `k`. pfloat has no `emax`
//!   or `emin`, so the exponent is an `i64` that *saturates* rather
//!   than producing `±∞` or `±0`. A shift that would carry the
//!   exponent past `i64::MAX` clamps to `i64::MAX` and raises
//!   [`Status::OVERFLOW`]; a shift past `i64::MIN` clamps to
//!   `i64::MIN` and raises [`Status::UNDERFLOW`]. This is the same
//!   saturating contract [`crate::ops::mul`] applies (see its
//!   `mul_extreme_exponent_saturates_not_panics` regression): a
//!   saturated exponent is a finite clamped value, never `±∞`. Every
//!   non-saturating shift is exact (`Status::OK`).

use crate::big::BigFloat;
use crate::class::Class;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

impl BigFloat {
    /// Returns `self × 2^k`, exact whenever the shifted exponent stays
    /// within the `i64` range.
    ///
    /// The mantissa, sign, and precision are unchanged; only the
    /// unbiased binary exponent shifts by `k`. The returned
    /// [`Status`] is [`Status::OK`] for every exact shift,
    /// [`Status::OVERFLOW`] when the exponent saturates `i64::MAX`,
    /// and [`Status::UNDERFLOW`] when it saturates `i64::MIN`. A
    /// signaling-NaN operand additionally raises
    /// [`Status::INVALID`]; `±0`, `±∞`, and quiet NaN pass through
    /// with [`Status::OK`].
    ///
    /// Unlike the rounding kernels this never allocates a new mantissa
    /// or changes precision, so it is `O(1)` in the limb count
    /// (modulo the unavoidable clone of `self`).
    ///
    /// # Examples
    ///
    /// ```
    /// # #[cfg(feature = "big")] {
    /// use pfloat::BigFloat;
    /// let three = BigFloat::try_from_i64_exact(3, 53).unwrap();
    /// let (twelve, status) = three.scale_by_pow2(2); // 3 × 2^2 = 12
    /// assert!(status.is_ok());
    /// let expected = BigFloat::try_from_i64_exact(12, 53).unwrap();
    /// assert_eq!(twelve.partial_cmp(&expected).0, Some(core::cmp::Ordering::Equal));
    /// # }
    /// ```
    #[must_use]
    pub fn scale_by_pow2(&self, k: i64) -> (Self, Status) {
        scale_kernel(self, k)
    }

    /// [`scale_by_pow2`](Self::scale_by_pow2) accumulating into a
    /// caller-supplied flag bag (the `no_std`-friendly variant).
    #[must_use]
    pub fn scale_by_pow2_with_flags(&self, k: i64, flags: &mut Status) -> Self {
        let (value, status) = self.scale_by_pow2(k);
        *flags |= status;
        value
    }
}

fn scale_kernel(a: &BigFloat, k: i64) -> (BigFloat, Status) {
    match &a.class {
        Class::Nan {
            quiet,
            sign,
            payload,
        } => {
            if *quiet {
                let nan = BigFloat::try_new_quiet_nan(*sign, a.precision, payload)
                    .expect("BigFloat invariant: precision >= 1");
                (nan, Status::OK)
            } else {
                // sNaN raises INVALID and propagates a quiet NaN.
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, a.precision, &[])
                    .expect("BigFloat invariant: precision >= 1");
                auto_raise(Status::INVALID);
                (nan, Status::INVALID)
            }
        }
        // ±0 × 2^k = ±0; ±∞ × 2^k = ±∞. Sign preserved, exact.
        Class::Zero { .. } | Class::Infinity { .. } => (a.clone(), Status::OK),
        Class::Normal {
            sign,
            exponent,
            mantissa,
        } => {
            // pfloat has no emax/emin: the exponent saturates inside
            // the i64 range instead of producing ±∞/±0. Compute the
            // shifted exponent in i128 (the sum of two i64s cannot
            // overflow i128) and clamp, flagging OVERFLOW/UNDERFLOW on
            // saturation. This mirrors `ops::mul`'s no-emax contract.
            let wide = i128::from(*exponent) + i128::from(k);
            let (new_exponent, status) = if wide > i128::from(i64::MAX) {
                (i64::MAX, Status::OVERFLOW)
            } else if wide < i128::from(i64::MIN) {
                (i64::MIN, Status::UNDERFLOW)
            } else {
                (wide as i64, Status::OK)
            };
            let scaled = BigFloat {
                class: Class::Normal {
                    sign: *sign,
                    exponent: new_exponent,
                    mantissa: mantissa.clone(),
                },
                precision: a.precision,
            };
            auto_raise(status);
            (scaled, status)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rounding::RoundingMode;
    use core::cmp::Ordering;

    fn from_i64(n: i64, p: u32) -> BigFloat {
        BigFloat::try_from_i64_exact(n, p).unwrap()
    }

    fn eq(a: &BigFloat, b: &BigFloat) -> bool {
        a.partial_cmp(b).0 == Some(Ordering::Equal) && a.precision() == b.precision()
    }

    #[test]
    fn scale_up_is_exact() {
        let three = from_i64(3, 53);
        let (twelve, s) = three.scale_by_pow2(2);
        assert!(s.is_ok());
        assert!(eq(&twelve, &from_i64(12, 53)));
    }

    #[test]
    fn scale_down_is_exact() {
        let twelve = from_i64(12, 53);
        let (three, s) = twelve.scale_by_pow2(-2);
        assert!(s.is_ok());
        assert!(eq(&three, &from_i64(3, 53)));
    }

    #[test]
    fn scale_by_zero_is_identity() {
        let x = from_i64(-7, 113);
        let (y, s) = x.scale_by_pow2(0);
        assert!(s.is_ok());
        assert_eq!(x, y); // structurally identical, not merely equal
    }

    #[test]
    fn scale_preserves_precision_and_mantissa() {
        // 5 at precision 53: mantissa 5<<61, exponent 2. After ×2^100
        // the mantissa is byte-for-byte identical, only the exponent
        // moves.
        let five = from_i64(5, 53);
        let (scaled, s) = five.scale_by_pow2(100);
        assert!(s.is_ok());
        match (&five.class, &scaled.class) {
            (
                Class::Normal {
                    mantissa: m0,
                    exponent: e0,
                    ..
                },
                Class::Normal {
                    mantissa: m1,
                    exponent: e1,
                    ..
                },
            ) => {
                assert_eq!(m0, m1, "mantissa must be untouched");
                assert_eq!(*e1, *e0 + 100);
            }
            _ => panic!("expected Normal"),
        }
        assert_eq!(scaled.precision(), 53);
    }

    #[test]
    fn round_trip_up_then_down() {
        let x = from_i64(42, 200);
        let (up, _) = x.scale_by_pow2(1_000_000);
        let (back, _) = up.scale_by_pow2(-1_000_000);
        assert_eq!(x, back);
    }

    #[test]
    fn overflow_saturates_to_i64_max_finite() {
        let x = from_i64(3, 53);
        let (y, s) = x.scale_by_pow2(i64::MAX);
        assert!(s.overflow());
        assert!(!s.underflow());
        // Saturated finite, never ±∞ (pfloat has no emax).
        assert!(!y.is_infinite());
        assert!(!y.is_nan());
        match &y.class {
            Class::Normal { exponent, .. } => assert_eq!(*exponent, i64::MAX),
            _ => panic!("expected Normal"),
        }
    }

    #[test]
    fn underflow_saturates_to_i64_min_finite() {
        // 0.5 has exponent -1, so shifting by i64::MIN carries the
        // exponent below i64::MIN and must saturate (a value with a
        // non-negative exponent would land at i64::MIN exactly with no
        // saturation).
        let (x, _) = from_i64(1, 53).scale_by_pow2(-1);
        let (y, s) = x.scale_by_pow2(i64::MIN);
        assert!(s.underflow());
        assert!(!s.overflow());
        assert!(!y.is_zero());
        assert!(!y.is_nan());
        match &y.class {
            Class::Normal { exponent, .. } => assert_eq!(*exponent, i64::MIN),
            _ => panic!("expected Normal"),
        }
    }

    #[test]
    fn no_saturation_just_below_the_edge() {
        // exponent of `1` is 0; shifting by i64::MAX lands exactly on
        // i64::MAX with no saturation.
        let one = from_i64(1, 53);
        let (y, s) = one.scale_by_pow2(i64::MAX);
        assert!(s.is_ok());
        match &y.class {
            Class::Normal { exponent, .. } => assert_eq!(*exponent, i64::MAX),
            _ => panic!("expected Normal"),
        }
    }

    #[test]
    fn zero_passes_through_with_sign() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (a, sa) = pz.scale_by_pow2(1000);
        let (b, sb) = nz.scale_by_pow2(-1000);
        assert!(sa.is_ok() && sb.is_ok());
        assert!(a.is_zero() && a.is_sign_positive());
        assert!(b.is_zero() && b.is_sign_negative());
    }

    #[test]
    fn infinity_passes_through_with_sign() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (y, s) = ni.scale_by_pow2(i64::MAX);
        assert!(s.is_ok());
        assert!(y.is_infinite() && y.is_sign_negative());
    }

    #[test]
    fn quiet_nan_propagates_without_flag() {
        let q = BigFloat::try_new_quiet_nan(Sign::Negative, 53, &[7]).unwrap();
        let (y, s) = q.scale_by_pow2(3);
        assert!(s.is_ok());
        assert!(y.is_quiet_nan());
        assert!(y.is_sign_negative());
    }

    #[test]
    fn signaling_nan_raises_invalid_and_quiets() {
        let s_nan = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (y, s) = s_nan.scale_by_pow2(3);
        assert!(s.invalid());
        assert!(y.is_quiet_nan());
    }

    #[test]
    fn with_flags_accumulates() {
        let x = from_i64(3, 53);
        let mut flags = Status::OK;
        let _ = x.scale_by_pow2_with_flags(i64::MAX, &mut flags);
        assert!(flags.overflow());
        // A subsequent exact scale does not clear the prior flag.
        let _ = x.scale_by_pow2_with_flags(1, &mut flags);
        assert!(flags.overflow());
    }

    #[test]
    fn scaling_agrees_with_multiplication() {
        // scale_by_pow2(k) must equal multiplication by the BigFloat
        // 2^k for any in-range k.
        for &k in &[1i64, 5, 13, 64, -1, -5, -64] {
            let x = from_i64(7, 64);
            let (scaled, _) = x.scale_by_pow2(k);
            // 2^k as a BigFloat: 1.0 with its exponent shifted by k
            // (exact for any in-range k, including negative).
            let (factor, _) = from_i64(1, 64).scale_by_pow2(k);
            let (product, _) = x.mul(&factor, RoundingMode::NearestEven);
            assert!(eq(&scaled, &product), "k={k}");
        }
    }
}
