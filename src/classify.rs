//! IEEE 754-2019 classification predicates and the [`IeeeClass`]
//! ten-variant enum (§5.7.2 / §6.2 / §6.3).
//!
//! pfloat's `IeeeClass` mirrors ferrodec's at
//! `ferrodec/src/classify.rs::IeeeClass` so the differential lane in
//! Phase 5 has a 1:1 enum translation when comparing pfloat values
//! against MPFR (via `gmp-mpfr-sys`) or against ferrodec's
//! `Decimal128` for cross-format sanity.
//!
//! pfloat does **not** have subnormals at arbitrary precision: there
//! is no implicit minimum exponent, so every finite non-zero value
//! is "normal" in IEEE terms. The
//! [`NegativeSubnormal`](IeeeClass::NegativeSubnormal) and
//! [`PositiveSubnormal`](IeeeClass::PositiveSubnormal) variants
//! exist for ABI parity with the ferrodec / MPFR / `f64::classify`
//! conventions but are never produced by pfloat values. Predicate
//! [`is_subnormal`](crate::big::BigFloat::is_subnormal) returns
//! `false` always.

use core::num::FpCategory;

#[cfg(feature = "big")]
use crate::big::BigFloat;
#[cfg(feature = "big")]
use crate::class::Class;
#[cfg(feature = "big")]
use crate::sign::Sign;

/// Ten-variant classification of an IEEE 754-2019 value.
///
/// The variant order matches ferrodec's `IeeeClass` for differential
/// translation. The two subnormal variants are unused by pfloat
/// values (pfloat has no subnormals); they are kept for ABI parity.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum IeeeClass {
    SignalingNaN,
    QuietNaN,
    NegativeInfinity,
    NegativeNormal,
    /// Unused for pfloat: pfloat has no subnormals at arbitrary
    /// precision. Retained for ABI parity with ferrodec's enum.
    NegativeSubnormal,
    NegativeZero,
    PositiveZero,
    /// Unused for pfloat. See [`NegativeSubnormal`](Self::NegativeSubnormal).
    PositiveSubnormal,
    PositiveNormal,
    PositiveInfinity,
}

#[cfg(feature = "big")]
impl BigFloat {
    /// `true` iff the value is a NaN (quiet or signaling).
    #[inline]
    #[must_use]
    pub fn is_nan(&self) -> bool {
        matches!(self.class, Class::Nan { .. })
    }

    /// `true` iff the value is a signaling NaN.
    #[inline]
    #[must_use]
    pub fn is_signaling_nan(&self) -> bool {
        matches!(self.class, Class::Nan { quiet: false, .. })
    }

    /// `true` iff the value is a quiet NaN.
    #[inline]
    #[must_use]
    pub fn is_quiet_nan(&self) -> bool {
        matches!(self.class, Class::Nan { quiet: true, .. })
    }

    /// `true` iff the value is `+∞` or `-∞`.
    #[inline]
    #[must_use]
    pub fn is_infinite(&self) -> bool {
        matches!(self.class, Class::Infinity { .. })
    }

    /// `true` iff the value is finite (zero or normal). False for
    /// infinities and NaNs.
    #[inline]
    #[must_use]
    pub fn is_finite(&self) -> bool {
        matches!(self.class, Class::Zero { .. } | Class::Normal { .. })
    }

    /// `true` iff the value is `+0` or `-0`.
    #[inline]
    #[must_use]
    pub fn is_zero(&self) -> bool {
        matches!(self.class, Class::Zero { .. })
    }

    /// `true` iff the value is a finite non-zero number.
    ///
    /// Every finite non-zero pfloat value is "normal" in IEEE terms.
    /// pfloat has no subnormals at arbitrary precision.
    #[inline]
    #[must_use]
    pub fn is_normal(&self) -> bool {
        matches!(self.class, Class::Normal { .. })
    }

    /// Always `false` for pfloat values.
    ///
    /// pfloat has no implicit minimum exponent, so there is no
    /// "subnormal" range. Retained for API parity with `f64`'s
    /// [`is_subnormal`](f64::is_subnormal); calling it makes sense
    /// when generic code wants the `f64`-shape interface.
    #[inline]
    #[must_use]
    pub fn is_subnormal(&self) -> bool {
        false
    }

    /// `true` iff the value's sign is [`Sign::Negative`].
    ///
    /// Includes `-0`, `-∞`, and negative-signed NaNs.
    #[inline]
    #[must_use]
    pub fn is_sign_negative(&self) -> bool {
        matches!(self.sign(), Sign::Negative)
    }

    /// `true` iff the value's sign is [`Sign::Positive`].
    ///
    /// Includes `+0`, `+∞`, and positive-signed NaNs.
    #[inline]
    #[must_use]
    pub fn is_sign_positive(&self) -> bool {
        matches!(self.sign(), Sign::Positive)
    }

    /// Returns the value's sign attribute.
    #[inline]
    #[must_use]
    pub fn sign(&self) -> Sign {
        match &self.class {
            Class::Zero { sign }
            | Class::Infinity { sign }
            | Class::Nan { sign, .. }
            | Class::Normal { sign, .. } => *sign,
        }
    }

    /// Coarse five-category classification matching
    /// [`f64::classify`](f64::classify).
    ///
    /// Quiet and signaling NaNs both map to [`FpCategory::Nan`].
    /// pfloat has no subnormals; this method never returns
    /// [`FpCategory::Subnormal`].
    #[inline]
    #[must_use]
    pub fn classify(&self) -> FpCategory {
        match &self.class {
            Class::Nan { .. } => FpCategory::Nan,
            Class::Infinity { .. } => FpCategory::Infinite,
            Class::Zero { .. } => FpCategory::Zero,
            Class::Normal { .. } => FpCategory::Normal,
        }
    }

    /// Returns the absolute value (sign forced to
    /// [`Sign::Positive`]).
    ///
    /// Preserves NaN payloads and quiet/signaling state per IEEE
    /// 754-2019 §6.3 ("`abs` is a non-arithmetic sign-bit
    /// operation"). Does not raise [`Status::INVALID`] on signaling
    /// NaN, unlike most arithmetic.
    ///
    /// [`Status::INVALID`]: crate::status::Status::INVALID
    #[inline]
    #[must_use]
    pub fn abs(&self) -> Self {
        Self {
            class: with_sign(&self.class, Sign::Positive),
            precision: self.precision,
        }
    }

    /// Returns the negation (sign flipped).
    ///
    /// Preserves NaN payloads and quiet/signaling state per IEEE
    /// 754-2019 §6.3.
    #[inline]
    #[must_use]
    pub fn negated(&self) -> Self {
        let new_sign = self.sign().flip();
        Self {
            class: with_sign(&self.class, new_sign),
            precision: self.precision,
        }
    }

    /// Returns a copy of `self` with the sign of `src`, per IEEE
    /// 754-2019 §6.3.
    ///
    /// Preserves NaN payloads and quiet/signaling state. Does not
    /// raise [`Status::INVALID`] on signaling NaN inputs.
    ///
    /// [`Status::INVALID`]: crate::status::Status::INVALID
    #[inline]
    #[must_use]
    pub fn copysign(&self, src: &Self) -> Self {
        Self {
            class: with_sign(&self.class, src.sign()),
            precision: self.precision,
        }
    }

    /// Returns `±1` matching `self`'s sign for every non-NaN value
    /// (including `±0` and `±∞`). For NaN, returns the input
    /// unchanged (preserving the payload, sign, and quiet/signaling
    /// state).
    ///
    /// Matches [`f64::signum`](f64::signum)'s shape: `±0.signum()`
    /// is `±1.0`, `±∞.signum()` is `±1.0`. Differs from MPFR's
    /// `mpfr_sgn`, which returns `0` for `±0`. Use
    /// [`is_zero`](Self::is_zero) when the zero-vs-nonzero
    /// distinction matters.
    #[inline]
    #[must_use]
    pub fn signum(&self) -> Self {
        match &self.class {
            Class::Nan { .. } => self.clone(),
            Class::Zero { sign } | Class::Infinity { sign } | Class::Normal { sign, .. } => {
                // `try_from_i64_exact(1, precision)` always
                // succeeds when precision >= 1 (the BigFloat
                // invariant); precision must be at least 1 because
                // every constructor validates it.
                let one = BigFloat::try_from_i64_exact(1, self.precision)
                    .expect("BigFloat invariant: precision >= 1");
                if matches!(sign, Sign::Negative) {
                    one.negated()
                } else {
                    one
                }
            }
        }
    }

    /// Ten-variant IEEE 754-2019 classification.
    #[inline]
    #[must_use]
    pub fn ieee_class(&self) -> IeeeClass {
        match &self.class {
            Class::Nan { quiet: false, .. } => IeeeClass::SignalingNaN,
            Class::Nan { quiet: true, .. } => IeeeClass::QuietNaN,
            Class::Infinity {
                sign: Sign::Negative,
            } => IeeeClass::NegativeInfinity,
            Class::Infinity {
                sign: Sign::Positive,
            } => IeeeClass::PositiveInfinity,
            Class::Zero {
                sign: Sign::Negative,
            } => IeeeClass::NegativeZero,
            Class::Zero {
                sign: Sign::Positive,
            } => IeeeClass::PositiveZero,
            Class::Normal {
                sign: Sign::Negative,
                ..
            } => IeeeClass::NegativeNormal,
            Class::Normal {
                sign: Sign::Positive,
                ..
            } => IeeeClass::PositiveNormal,
        }
    }
}

/// Returns a [`Class`] with the same kind and fields as the input
/// but with the sign replaced.
///
/// Used by [`BigFloat::abs`], [`BigFloat::negated`], and
/// [`BigFloat::copysign`].
#[cfg(feature = "big")]
fn with_sign(class: &Class, new_sign: Sign) -> Class {
    match class {
        Class::Zero { .. } => Class::Zero { sign: new_sign },
        Class::Infinity { .. } => Class::Infinity { sign: new_sign },
        Class::Nan { quiet, payload, .. } => Class::Nan {
            quiet: *quiet,
            sign: new_sign,
            payload: payload.clone(),
        },
        Class::Normal {
            exponent, mantissa, ..
        } => Class::Normal {
            sign: new_sign,
            exponent: *exponent,
            mantissa: mantissa.clone(),
        },
    }
}

#[cfg(test)]
#[cfg(feature = "big")]
mod tests {
    use super::*;

    fn at_each_precision<F: Fn(u32)>(f: F) {
        for &p in &[1, 53, 113, 256] {
            f(p);
        }
    }

    #[test]
    fn nan_classification() {
        at_each_precision(|p| {
            let q = BigFloat::try_new_quiet_nan(Sign::Positive, p, &[]).unwrap();
            let s = BigFloat::try_new_signaling_nan(Sign::Positive, p, &[]).unwrap();

            assert!(q.is_nan());
            assert!(q.is_quiet_nan());
            assert!(!q.is_signaling_nan());
            assert!(!q.is_finite());
            assert!(!q.is_infinite());
            assert!(!q.is_zero());
            assert!(!q.is_normal());

            assert!(s.is_nan());
            assert!(s.is_signaling_nan());
            assert!(!s.is_quiet_nan());
            assert!(!s.is_finite());
        });
    }

    #[test]
    fn infinity_classification() {
        at_each_precision(|p| {
            let pos = BigFloat::try_new_infinity(Sign::Positive, p).unwrap();
            let neg = BigFloat::try_new_infinity(Sign::Negative, p).unwrap();

            assert!(pos.is_infinite());
            assert!(!pos.is_finite());
            assert!(!pos.is_nan());
            assert!(pos.is_sign_positive());

            assert!(neg.is_infinite());
            assert!(neg.is_sign_negative());
        });
    }

    #[test]
    fn zero_classification() {
        at_each_precision(|p| {
            let pos = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
            let neg = BigFloat::try_new_zero(Sign::Negative, p).unwrap();

            assert!(pos.is_zero());
            assert!(pos.is_finite());
            assert!(!pos.is_normal());
            assert!(pos.is_sign_positive());

            assert!(neg.is_zero());
            assert!(neg.is_sign_negative());
        });
    }

    #[test]
    fn one_is_normal() {
        for &p in &[1, 53, 113, 256] {
            let one = BigFloat::try_from_i64_exact(1, p).unwrap();
            assert!(one.is_normal());
            assert!(one.is_finite());
            assert!(!one.is_zero());
            assert!(!one.is_infinite());
            assert!(!one.is_nan());
            assert!(!one.is_subnormal());
            assert!(one.is_sign_positive());
        }
    }

    #[test]
    fn ieee_class_covers_all_ten_variants_pfloat_can_produce() {
        // pfloat does not produce subnormals; the other 8 variants
        // are reachable.
        let qnan = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let snan = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let pinf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let ninf = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let pzero = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let nzero = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let pnormal = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let nnormal = BigFloat::try_from_i64_exact(-1, 53).unwrap();

        assert_eq!(qnan.ieee_class(), IeeeClass::QuietNaN);
        assert_eq!(snan.ieee_class(), IeeeClass::SignalingNaN);
        assert_eq!(pinf.ieee_class(), IeeeClass::PositiveInfinity);
        assert_eq!(ninf.ieee_class(), IeeeClass::NegativeInfinity);
        assert_eq!(pzero.ieee_class(), IeeeClass::PositiveZero);
        assert_eq!(nzero.ieee_class(), IeeeClass::NegativeZero);
        assert_eq!(pnormal.ieee_class(), IeeeClass::PositiveNormal);
        assert_eq!(nnormal.ieee_class(), IeeeClass::NegativeNormal);
    }

    #[test]
    fn classify_matches_fpcategory_coarsening() {
        let qnan = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        assert_eq!(qnan.classify(), FpCategory::Nan);
        let snan = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        assert_eq!(snan.classify(), FpCategory::Nan);
        let pinf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        assert_eq!(pinf.classify(), FpCategory::Infinite);
        let pzero = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        assert_eq!(pzero.classify(), FpCategory::Zero);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(one.classify(), FpCategory::Normal);
    }

    #[test]
    fn pfloat_never_returns_subnormal() {
        // Every shape of value pfloat can produce returns
        // is_subnormal == false.
        for v in [
            BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap(),
            BigFloat::try_new_infinity(Sign::Positive, 53).unwrap(),
            BigFloat::try_new_zero(Sign::Positive, 53).unwrap(),
            BigFloat::try_from_i64_exact(1, 53).unwrap(),
        ] {
            assert!(!v.is_subnormal());
        }
    }

    #[test]
    fn sign_extraction_works_across_kinds() {
        assert_eq!(
            BigFloat::try_new_zero(Sign::Negative, 53).unwrap().sign(),
            Sign::Negative
        );
        assert_eq!(
            BigFloat::try_new_infinity(Sign::Negative, 53)
                .unwrap()
                .sign(),
            Sign::Negative
        );
        assert_eq!(
            BigFloat::try_new_quiet_nan(Sign::Negative, 53, &[])
                .unwrap()
                .sign(),
            Sign::Negative
        );
        assert_eq!(
            BigFloat::try_from_i64_exact(-1, 53).unwrap().sign(),
            Sign::Negative
        );
    }

    #[test]
    fn abs_clears_sign() {
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        assert!(neg_one.abs().is_sign_positive());
        let neg_zero = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        assert!(neg_zero.abs().is_sign_positive());
        let neg_inf = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        assert!(neg_inf.abs().is_sign_positive());
        let neg_nan = BigFloat::try_new_quiet_nan(Sign::Negative, 53, &[7]).unwrap();
        let abs_nan = neg_nan.abs();
        assert!(abs_nan.is_sign_positive());
        assert!(abs_nan.is_quiet_nan());
    }

    #[test]
    fn negated_flips_sign() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert!(one.negated().is_sign_negative());
        let pos_zero = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        assert!(pos_zero.negated().is_sign_negative());
        let qnan = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        assert!(qnan.negated().is_sign_negative());
        // Involution: negated.negated() == original.
        assert_eq!(one.negated().negated(), one);
    }

    #[test]
    fn negated_preserves_signaling_ness() {
        let snan = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        assert!(snan.negated().is_signaling_nan());
    }

    #[test]
    fn copysign_takes_sign_from_arg() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let neg_two = BigFloat::try_from_i64_exact(-2, 53).unwrap();
        let combined = one.copysign(&neg_two);
        assert!(combined.is_sign_negative());
        assert!(combined.is_normal());
        // Magnitude is preserved (still 1, just with negative sign).
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        assert_eq!(combined, neg_one);
    }

    #[test]
    fn copysign_preserves_nan_payload() {
        let n = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[42]).unwrap();
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        let result = n.copysign(&neg_one);
        assert!(result.is_quiet_nan());
        assert!(result.is_sign_negative());
        match result.class {
            Class::Nan { payload, .. } => assert_eq!(payload[0], 42),
            _ => panic!("expected NaN"),
        }
    }

    #[test]
    fn signum_basics() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(one.signum(), one);
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert_eq!(two.signum(), one);
        let neg_three = BigFloat::try_from_i64_exact(-3, 53).unwrap();
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        assert_eq!(neg_three.signum(), neg_one);
        // Infinity → ±1.
        let pinf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        assert_eq!(pinf.signum(), one);
        let ninf = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        assert_eq!(ninf.signum(), neg_one);
    }

    #[test]
    fn signum_zero_returns_signed_one() {
        // f64-shape semantics: signum(±0) == ±1.
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        assert_eq!(pz.signum(), one);
        assert_eq!(nz.signum(), neg_one);
    }

    #[test]
    fn signum_nan_preserves_payload() {
        let n = BigFloat::try_new_quiet_nan(Sign::Negative, 53, &[99]).unwrap();
        let s = n.signum();
        assert!(s.is_quiet_nan());
        assert!(s.is_sign_negative());
        match s.class {
            Class::Nan { payload, .. } => assert_eq!(payload[0], 99),
            _ => panic!("expected NaN"),
        }
    }
}
