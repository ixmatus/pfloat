//! [`RealScalar`]: the sealed trait of real scalars a [`Complex`] component
//! can be. ADR-0089.
//!
//! Sealed — implemented only for [`pfloat::BigFloat`] and
//! [`pfloat::FixedFloat<PREC>`], so a third party cannot instantiate a
//! `Complex` over an unverified scalar through *this* crate's surface. The
//! seal is scoped to pfloat-complex, not universal (see the crate docs).
//!
//! This cut presents the subset the componentwise additive arithmetic
//! needs (`add`, `sub`, `neg`, and the predicates the tests use). The trait
//! is sealed, so later slices extend it (the `mul_add_mul` / `mul_sub_mul`
//! forms for complex multiply and divide, `hypot` and `atan2` for magnitude
//! and phase, the elementary kernels) without a breaking change.
//!
//! [`Complex`]: crate::Complex

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Status};

mod sealed {
    /// Private supertrait that seals [`RealScalar`](super::RealScalar):
    /// an external crate cannot name it, so it cannot add an impl.
    pub trait Sealed {}
}

/// The real scalar types a [`Complex`](crate::Complex) component can be.
///
/// Sealed. Each method delegates to the inherent pfloat kernel of the same
/// name, so a `Complex` operation is correctly rounded componentwise by
/// construction. Every fallible-on-precision pfloat method is presented in
/// its always-valid form (a `RealScalar` value always has precision `>= 1`,
/// so the `BuildError` path is unreachable).
pub trait RealScalar: Clone + core::fmt::Debug + sealed::Sealed {
    /// Precision in bits.
    fn precision(&self) -> u32;
    /// `true` for NaN.
    fn is_nan(&self) -> bool;
    /// `true` for `±0`.
    fn is_zero(&self) -> bool;
    /// Negation (exact, sign-bit flip).
    fn negated(&self) -> Self;
    /// IEEE 754-2019 §5.11 partial comparison: `None` for a NaN operand;
    /// `Status::INVALID` on a signaling-NaN comparand.
    fn compare(&self, other: &Self) -> (Option<Ordering>, Status);
    /// `self + other`, correctly rounded under `mode` to the result
    /// precision (the larger operand precision; for `FixedFloat<PREC>`,
    /// `PREC`).
    fn add(&self, other: &Self, mode: RoundingMode) -> (Self, Status);
    /// `self − other`, correctly rounded under `mode`.
    fn sub(&self, other: &Self, mode: RoundingMode) -> (Self, Status);
}

impl sealed::Sealed for BigFloat {}

impl RealScalar for BigFloat {
    #[inline]
    fn precision(&self) -> u32 {
        BigFloat::precision(self)
    }
    #[inline]
    fn is_nan(&self) -> bool {
        BigFloat::is_nan(self)
    }
    #[inline]
    fn is_zero(&self) -> bool {
        BigFloat::is_zero(self)
    }
    #[inline]
    fn negated(&self) -> Self {
        BigFloat::negated(self)
    }
    #[inline]
    fn compare(&self, other: &Self) -> (Option<Ordering>, Status) {
        BigFloat::partial_cmp(self, other)
    }
    #[inline]
    fn add(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        BigFloat::add(self, other, mode)
    }
    #[inline]
    fn sub(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        BigFloat::sub(self, other, mode)
    }
}

#[cfg(feature = "fixed")]
mod fixed_impl {
    use super::{sealed, Ordering, RealScalar, RoundingMode, Status};
    use pfloat::{limbs_for, FixedFloat};

    impl<const PREC: u32> sealed::Sealed for FixedFloat<PREC> where [(); limbs_for(PREC)]: {}

    impl<const PREC: u32> RealScalar for FixedFloat<PREC>
    where
        [(); limbs_for(PREC)]:,
    {
        #[inline]
        fn precision(&self) -> u32 {
            PREC
        }
        #[inline]
        fn is_nan(&self) -> bool {
            FixedFloat::is_nan(self)
        }
        #[inline]
        fn is_zero(&self) -> bool {
            FixedFloat::is_zero(self)
        }
        #[inline]
        fn negated(&self) -> Self {
            FixedFloat::negated(*self)
        }
        #[inline]
        fn compare(&self, other: &Self) -> (Option<Ordering>, Status) {
            FixedFloat::partial_cmp(self, other)
        }
        #[inline]
        fn add(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
            FixedFloat::add(self, other, mode)
        }
        #[inline]
        fn sub(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
            FixedFloat::sub(self, other, mode)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Exercises the trait through `T: RealScalar`, the way the complex
    // kernels reach the scalar engine. If any delegation recursed (the
    // inherent-vs-trait name clash), this would overflow the stack.
    fn add_sub_round_trips<T: RealScalar>(a: &T, b: &T) -> (Option<Ordering>, Status) {
        let (sum, _) = a.add(b, RoundingMode::NearestEven);
        let (back, _) = sum.sub(b, RoundingMode::NearestEven);
        back.compare(a)
    }

    #[test]
    fn bigfloat_delegates_without_recursion() {
        let a = BigFloat::try_from_i64_exact(7, 64).unwrap();
        let b = BigFloat::try_from_i64_exact(3, 64).unwrap();
        let (ord, st) = add_sub_round_trips(&a, &b);
        assert_eq!(ord, Some(Ordering::Equal));
        assert!(st.is_ok());
        assert_eq!(a.precision(), 64);
        assert!(!a.is_nan() && !a.is_zero());
        assert_eq!(a.negated().compare(&a).0, Some(Ordering::Less));
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixedfloat_delegates_without_recursion() {
        use pfloat::FixedFloat;
        let a = FixedFloat::<64>::try_from_i64_exact(7).unwrap();
        let b = FixedFloat::<64>::try_from_i64_exact(3).unwrap();
        let (ord, st) = add_sub_round_trips(&a, &b);
        assert_eq!(ord, Some(Ordering::Equal));
        assert!(st.is_ok());
        assert_eq!(a.precision(), 64);
    }
}
