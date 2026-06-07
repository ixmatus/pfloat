//! [`RealScalar`]: the sealed trait of midpoint scalars a [`Ball`] can
//! wrap. ADR-0075.
//!
//! `Ball<T>` is generic over `T: RealScalar`, and `RealScalar` is
//! *sealed*: only [`pfloat::BigFloat`] and [`pfloat::FixedFloat<PREC>`]
//! implement it, and a third party cannot add an impl. That makes
//! "the midpoint is a correctly-rounded pfloat scalar" a fact the ball
//! crate's own surface cannot be made to break — a `Ball` can never be
//! instantiated over an unverified or wrongly-rounded scalar type.
//!
//! The seal is scoped, not universal: because Phase 3 already shipped
//! `num_traits::Num` for `FixedFloat<PREC>` (pfloat ADR-0070), a third
//! party can still build, say, a `num_complex::Complex<FixedFloat<P>>`
//! outside this crate. `RealScalar` closes *pfloat-ball's* inhabitant
//! set, not the universe of generic numeric code.
//!
//! [`Ball`]: crate::Ball

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Status};

use crate::mag::Mag;

mod sealed {
    /// Private supertrait that seals [`RealScalar`](super::RealScalar):
    /// an external crate cannot name it, so it cannot add an impl.
    pub trait Sealed {}
}

/// The scalar types a [`Ball`](crate::Ball) can use as its midpoint.
///
/// Sealed — implemented only for [`pfloat::BigFloat`] and
/// [`pfloat::FixedFloat<PREC>`]. The methods are the surface a ball
/// kernel needs: the directed arithmetic that brackets a result, the
/// adjacency / scaling primitives that build and convert radii, and the
/// two bridges to [`Mag`]. Every fallible-on-precision pfloat method is
/// presented in its always-valid form here (a `RealScalar` value always
/// has precision `≥ 1`, so the `BuildError` path is unreachable).
pub trait RealScalar: Clone + core::fmt::Debug + sealed::Sealed {
    /// Precision in bits.
    fn precision(&self) -> u32;

    /// `+0` at the given precision (for `FixedFloat<PREC>` the argument
    /// is ignored and `PREC` is used). The finite placeholder midpoint of
    /// an entire-result ball.
    fn zero(precision: u32) -> Self;

    /// `true` for a finite value (zero or normal).
    fn is_finite(&self) -> bool;
    /// `true` for NaN.
    fn is_nan(&self) -> bool;
    /// `true` for `±0`.
    fn is_zero(&self) -> bool;
    /// `true` for `±∞`.
    fn is_infinite(&self) -> bool;
    /// `true` for a negative sign bit.
    fn is_sign_negative(&self) -> bool;

    /// Negation (exact, sign-bit flip).
    fn negated(&self) -> Self;
    /// Absolute value (exact, sign forced positive).
    fn abs(&self) -> Self;

    /// IEEE 754-2019 §5.11 partial comparison: `None` for a NaN operand;
    /// `Status::INVALID` on a signaling-NaN comparand.
    fn compare(&self, other: &Self) -> (Option<Ordering>, Status);

    /// `self + other`, correctly rounded under `mode` to the result
    /// precision (the larger operand precision; for `FixedFloat<PREC>`,
    /// `PREC`).
    fn add(&self, other: &Self, mode: RoundingMode) -> (Self, Status);
    /// `self − other`, correctly rounded under `mode`.
    fn sub(&self, other: &Self, mode: RoundingMode) -> (Self, Status);
    /// `self · other`, correctly rounded under `mode`.
    fn mul(&self, other: &Self, mode: RoundingMode) -> (Self, Status);
    /// `self / other`, correctly rounded under `mode`.
    fn div(&self, other: &Self, mode: RoundingMode) -> (Self, Status);
    /// `√self`, correctly rounded under `mode`.
    fn sqrt(&self, mode: RoundingMode) -> (Self, Status);

    /// `self · 2^k`, exact (saturating at the `i64` exponent range).
    fn scale_by_pow2(&self, k: i64) -> (Self, Status);
    /// The least representable value greater than `self`.
    fn next_up(&self) -> (Self, Status);
    /// The greatest representable value less than `self`.
    fn next_down(&self) -> (Self, Status);
    /// The unit in the last place at `self`.
    fn ulp(&self) -> (Self, Status);

    /// `|self|` narrowed *up* to a single-limb [`Mag`] (sound upper
    /// bound).
    fn magnitude_to_mag(&self) -> Mag;

    /// A radius [`Mag`] widened to a scalar of this type, rounded *up*
    /// so the result is `≥` the radius (sound for outward endpoints).
    /// Exact when the type's precision is `≥ 64`.
    fn radius_to_scalar(radius: Mag) -> Self;
}

impl sealed::Sealed for BigFloat {}

impl RealScalar for BigFloat {
    #[inline]
    fn precision(&self) -> u32 {
        BigFloat::precision(self)
    }
    #[inline]
    fn zero(precision: u32) -> Self {
        BigFloat::try_new_zero(pfloat::Sign::Positive, precision).expect("precision ≥ 1")
    }
    #[inline]
    fn is_finite(&self) -> bool {
        BigFloat::is_finite(self)
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
    fn is_infinite(&self) -> bool {
        BigFloat::is_infinite(self)
    }
    #[inline]
    fn is_sign_negative(&self) -> bool {
        BigFloat::is_sign_negative(self)
    }
    #[inline]
    fn negated(&self) -> Self {
        BigFloat::negated(self)
    }
    #[inline]
    fn abs(&self) -> Self {
        BigFloat::abs(self)
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
    fn mul(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        BigFloat::mul(self, other, mode)
    }
    #[inline]
    fn div(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        BigFloat::div(self, other, mode)
    }
    #[inline]
    fn sqrt(&self, mode: RoundingMode) -> (Self, Status) {
        BigFloat::sqrt(self, mode)
    }
    #[inline]
    fn scale_by_pow2(&self, k: i64) -> (Self, Status) {
        BigFloat::scale_by_pow2(self, k)
    }
    #[inline]
    fn next_up(&self) -> (Self, Status) {
        BigFloat::next_up(self)
    }
    #[inline]
    fn next_down(&self) -> (Self, Status) {
        BigFloat::next_down(self)
    }
    #[inline]
    fn ulp(&self) -> (Self, Status) {
        BigFloat::ulp(self)
    }
    #[inline]
    fn magnitude_to_mag(&self) -> Mag {
        Mag::from_bigfloat_ceil(self)
    }
    #[inline]
    fn radius_to_scalar(radius: Mag) -> Self {
        // to_bigfloat is exact at precision 64; the radius is positive so
        // no rounding is needed (the "up" contract is vacuously met).
        radius.to_bigfloat()
    }
}

#[cfg(feature = "fixed")]
mod fixed_impl {
    use super::{sealed, Mag, Ordering, RealScalar, RoundingMode, Status};
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
        fn zero(_precision: u32) -> Self {
            FixedFloat::zero()
        }
        #[inline]
        fn is_finite(&self) -> bool {
            FixedFloat::is_finite(self)
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
        fn is_infinite(&self) -> bool {
            FixedFloat::is_infinite(self)
        }
        #[inline]
        fn is_sign_negative(&self) -> bool {
            FixedFloat::is_sign_negative(self)
        }
        #[inline]
        fn negated(&self) -> Self {
            FixedFloat::negated(*self)
        }
        #[inline]
        fn abs(&self) -> Self {
            FixedFloat::abs(*self)
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
        fn mul(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
            FixedFloat::mul(self, other, mode)
        }
        #[inline]
        fn div(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
            FixedFloat::div(self, other, mode)
        }
        #[inline]
        fn sqrt(&self, mode: RoundingMode) -> (Self, Status) {
            FixedFloat::sqrt(self, mode)
        }
        #[inline]
        fn scale_by_pow2(&self, k: i64) -> (Self, Status) {
            FixedFloat::scale_by_pow2(self, k)
        }
        #[inline]
        fn next_up(&self) -> (Self, Status) {
            FixedFloat::next_up(self)
        }
        #[inline]
        fn next_down(&self) -> (Self, Status) {
            FixedFloat::next_down(self)
        }
        #[inline]
        fn ulp(&self) -> (Self, Status) {
            FixedFloat::ulp(self)
        }
        #[inline]
        fn magnitude_to_mag(&self) -> Mag {
            Mag::from_bigfloat_ceil(&self.to_big())
        }
        #[inline]
        fn radius_to_scalar(radius: Mag) -> Self {
            // Round the 64-bit radius up to PREC (TowardPositive on a
            // positive value), so the scalar radius is ≥ the Mag — sound
            // for outward endpoints. Exact when PREC ≥ 64.
            FixedFloat::try_from_big_round(&radius.to_bigfloat(), RoundingMode::TowardPositive).0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // A generic helper exercises the trait through `T: RealScalar`, which
    // is exactly how the ball kernels reach the scalar engine. If any
    // delegation recursed, these would overflow the stack.
    fn round_trip_add_sub<T: RealScalar>(a: &T, b: &T) -> (Option<Ordering>, Status) {
        let (sum, _) = a.add(b, RoundingMode::NearestEven);
        let (back, _) = sum.sub(b, RoundingMode::NearestEven);
        back.compare(a)
    }

    fn primitives<T: RealScalar>(x: &T) {
        let (up, _) = x.next_up();
        let (down, _) = x.next_down();
        assert_eq!(down.compare(x).0, Some(Ordering::Less));
        assert_eq!(x.compare(&up).0, Some(Ordering::Less));
        let (u, _) = x.ulp();
        assert!(u.magnitude_to_mag() != Mag::Zero);
        let (scaled, _) = x.scale_by_pow2(3);
        // x·8 > x for positive x.
        assert_eq!(x.compare(&scaled).0, Some(Ordering::Less));
    }

    #[test]
    fn bigfloat_impl_delegates_without_recursion() {
        let a = BigFloat::try_from_i64_exact(7, 64).unwrap();
        let b = BigFloat::try_from_i64_exact(3, 64).unwrap();
        let (ord, st) = round_trip_add_sub(&a, &b);
        assert_eq!(ord, Some(Ordering::Equal));
        assert!(st.is_ok());
        primitives(&a);

        // Mag bridges round-trip a small integer.
        let m = a.magnitude_to_mag();
        let back = <BigFloat as RealScalar>::radius_to_scalar(m);
        assert_eq!(back.compare(&a).0, Some(Ordering::Equal));
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixedfloat_impl_delegates_without_recursion() {
        use pfloat::FixedFloat;
        let a = FixedFloat::<64>::try_from_i64_exact(7).unwrap();
        let b = FixedFloat::<64>::try_from_i64_exact(3).unwrap();
        let (ord, st) = round_trip_add_sub(&a, &b);
        assert_eq!(ord, Some(Ordering::Equal));
        assert!(st.is_ok());
        primitives(&a);

        let m = a.magnitude_to_mag();
        let back = <FixedFloat<64> as RealScalar>::radius_to_scalar(m);
        assert_eq!(back.compare(&a).0, Some(Ordering::Equal));
    }
}
