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
    /// `self·b + c·d`, the fused two-product sum, correctly rounded under
    /// `mode` with a single rounding (ADR-0088). The complex multiply uses
    /// this for the `a·d + b·c` component.
    fn mul_add_mul(&self, b: &Self, c: &Self, d: &Self, mode: RoundingMode) -> (Self, Status);
    /// `self·b − c·d`, the fused two-product difference, correctly rounded
    /// under `mode` with a single rounding. The complex multiply uses this
    /// for the `a·c − b·d` component.
    fn mul_sub_mul(&self, b: &Self, c: &Self, d: &Self, mode: RoundingMode) -> (Self, Status);
    /// `hypot(self, other) = sqrt(self² + other²)`, correctly rounded under
    /// `mode` to `max(self, other)` precision (IEEE 754-2019 §9.2.1). The
    /// complex magnitude `abs` delegates to this; the scalar kernel's Ziv
    /// driver makes it correctly rounded (not the lossy `sqrt(self·self +
    /// other·other)`) and carries the infinity-dominates-NaN special cases.
    #[cfg(feature = "exp-log")]
    fn hypot(&self, other: &Self, mode: RoundingMode) -> (Self, Status);
    /// `atan2(self, x)`: the polar angle of `(x, self)` in `(−π, π]`,
    /// correctly rounded under `mode` to `max(self, x)` precision. The
    /// complex phase `arg` and the imaginary part of complex `log` delegate
    /// to this; it already carries the C99 Annex G signed-zero branch
    /// convention (the IEEE 754-2019 §9.2.1 dispatch, e.g. `atan2(+0, −0) =
    /// +π`, `atan2(−0, −0) = −π`).
    #[cfg(feature = "trig")]
    fn atan2(&self, x: &Self, mode: RoundingMode) -> (Self, Status);
    /// Convert to a [`pfloat::BigFloat`], exactly. The complex divide runs
    /// its directed-pair enclosure Ziv loop in `BigFloat` at a working
    /// precision *above* the output precision (which `FixedFloat<PREC>`
    /// cannot hold), so it bridges through this conversion.
    fn to_big(&self) -> pfloat::BigFloat;
    /// Convert a [`pfloat::BigFloat`] at this type's output precision back to
    /// this scalar type (for `FixedFloat<PREC>`, exactly `PREC` bits).
    fn from_big(b: &pfloat::BigFloat) -> Self;
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
    #[inline]
    fn mul_add_mul(&self, b: &Self, c: &Self, d: &Self, mode: RoundingMode) -> (Self, Status) {
        BigFloat::mul_add_mul(self, b, c, d, mode)
    }
    #[inline]
    fn mul_sub_mul(&self, b: &Self, c: &Self, d: &Self, mode: RoundingMode) -> (Self, Status) {
        BigFloat::mul_sub_mul(self, b, c, d, mode)
    }
    #[cfg(feature = "exp-log")]
    #[inline]
    fn hypot(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        // `BigFloat::hypot` rounds to `self.precision`; round to the larger
        // operand precision instead, matching `add`/`atan2` so the complex
        // magnitude carries the full working precision of its components.
        let p = BigFloat::precision(self).max(BigFloat::precision(other));
        BigFloat::hypot_round(self, other, p, mode).expect("precision >= 1 by RealScalar invariant")
    }
    #[cfg(feature = "trig")]
    #[inline]
    fn atan2(&self, x: &Self, mode: RoundingMode) -> (Self, Status) {
        BigFloat::atan2(self, x, mode)
    }
    #[inline]
    fn to_big(&self) -> BigFloat {
        self.clone()
    }
    #[inline]
    fn from_big(b: &BigFloat) -> Self {
        b.clone()
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
        #[inline]
        fn mul_add_mul(&self, b: &Self, c: &Self, d: &Self, mode: RoundingMode) -> (Self, Status) {
            FixedFloat::mul_add_mul(self, b, c, d, mode)
        }
        #[inline]
        fn mul_sub_mul(&self, b: &Self, c: &Self, d: &Self, mode: RoundingMode) -> (Self, Status) {
            FixedFloat::mul_sub_mul(self, b, c, d, mode)
        }
        #[cfg(feature = "exp-log")]
        #[inline]
        fn hypot(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
            FixedFloat::hypot(self, other, mode)
        }
        #[cfg(feature = "trig")]
        #[inline]
        fn atan2(&self, x: &Self, mode: RoundingMode) -> (Self, Status) {
            FixedFloat::atan2(self, x, mode)
        }
        #[inline]
        fn to_big(&self) -> pfloat::BigFloat {
            FixedFloat::to_big(self)
        }
        #[inline]
        fn from_big(b: &pfloat::BigFloat) -> Self {
            FixedFloat::try_from_big_exact(b.clone()).expect("BigFloat is at PREC")
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
