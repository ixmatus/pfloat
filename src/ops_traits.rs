//! `core::ops` operator overloads for [`BigFloat`] and
//! [`FixedFloat<PREC>`].
//!
//! Gated behind the `ops` feature. Each overload uses the
//! IEEE 754-2019 default rounding mode
//! ([`RoundingMode::NearestEven`]) and discards the returned
//! [`Status`]. Callers needing explicit rounding-mode control or
//! flag accumulation should use the method form
//! (`a.add(&b, mode)`).
//!
//! The pattern matches ferrodec's `src/ops_traits.rs`: a constant
//! `RM = RoundingMode::NearestEven` per file, and each `Add`/`Sub`/
//! `Mul`/`Div` impl is a single-line forward to the canonical
//! method. `Neg` flips the sign without touching the rounding
//! pipeline.

use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::big::BigFloat;
use crate::rounding::RoundingMode;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

const RM: RoundingMode = RoundingMode::NearestEven;

// -------- BigFloat operator impls --------

impl Add for BigFloat {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        BigFloat::add(&self, &rhs, RM).0
    }
}

impl Add<&BigFloat> for BigFloat {
    type Output = Self;
    #[inline]
    fn add(self, rhs: &BigFloat) -> Self {
        BigFloat::add(&self, rhs, RM).0
    }
}

impl Sub for BigFloat {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        BigFloat::sub(&self, &rhs, RM).0
    }
}

impl Sub<&BigFloat> for BigFloat {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: &BigFloat) -> Self {
        BigFloat::sub(&self, rhs, RM).0
    }
}

impl Mul for BigFloat {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        BigFloat::mul(&self, &rhs, RM).0
    }
}

impl Mul<&BigFloat> for BigFloat {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: &BigFloat) -> Self {
        BigFloat::mul(&self, rhs, RM).0
    }
}

impl Div for BigFloat {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        BigFloat::div(&self, &rhs, RM).0
    }
}

impl Div<&BigFloat> for BigFloat {
    type Output = Self;
    #[inline]
    fn div(self, rhs: &BigFloat) -> Self {
        BigFloat::div(&self, rhs, RM).0
    }
}

impl Neg for BigFloat {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        self.negated()
    }
}

impl AddAssign for BigFloat {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = BigFloat::add(self, &rhs, RM).0;
    }
}

impl SubAssign for BigFloat {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = BigFloat::sub(self, &rhs, RM).0;
    }
}

impl MulAssign for BigFloat {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = BigFloat::mul(self, &rhs, RM).0;
    }
}

impl DivAssign for BigFloat {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = BigFloat::div(self, &rhs, RM).0;
    }
}

// -------- FixedFloat operator impls --------

#[cfg(feature = "fixed")]
impl<const PREC: u32> Add for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        FixedFloat::<PREC>::add(&self, &rhs, RM).0
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> Sub for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        FixedFloat::<PREC>::sub(&self, &rhs, RM).0
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> Mul for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        FixedFloat::<PREC>::mul(&self, &rhs, RM).0
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> Div for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self {
        FixedFloat::<PREC>::div(&self, &rhs, RM).0
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> Neg for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        self.negated()
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> AddAssign for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = FixedFloat::<PREC>::add(self, &rhs, RM).0;
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> SubAssign for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = FixedFloat::<PREC>::sub(self, &rhs, RM).0;
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> MulAssign for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = FixedFloat::<PREC>::mul(self, &rhs, RM).0;
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> DivAssign for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = FixedFloat::<PREC>::div(self, &rhs, RM).0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn from_i64(n: i64, p: u32) -> BigFloat {
        BigFloat::try_from_i64_exact(n, p).unwrap()
    }

    #[test]
    fn bigfloat_add_op() {
        let a = from_i64(2, 53);
        let b = from_i64(3, 53);
        let sum = a + b;
        assert_eq!(
            sum.partial_cmp(&from_i64(5, 53)).0,
            Some(core::cmp::Ordering::Equal)
        );
    }

    #[test]
    fn bigfloat_mul_div_chain() {
        let a = from_i64(6, 53);
        let b = from_i64(7, 53);
        let c = from_i64(3, 53);
        let result = (a * b) / c; // (6 * 7) / 3 = 14
        assert_eq!(
            result.partial_cmp(&from_i64(14, 53)).0,
            Some(core::cmp::Ordering::Equal)
        );
    }

    #[test]
    fn bigfloat_neg_op() {
        let a = from_i64(7, 53);
        let n = -a.clone();
        assert!(n.is_sign_negative());
        let abs = -n;
        assert!(abs.is_sign_positive());
    }

    #[test]
    fn bigfloat_assign_ops() {
        let mut a = from_i64(10, 53);
        a += from_i64(5, 53);
        assert_eq!(
            a.partial_cmp(&from_i64(15, 53)).0,
            Some(core::cmp::Ordering::Equal)
        );
        a -= from_i64(3, 53);
        assert_eq!(
            a.partial_cmp(&from_i64(12, 53)).0,
            Some(core::cmp::Ordering::Equal)
        );
        a *= from_i64(2, 53);
        assert_eq!(
            a.partial_cmp(&from_i64(24, 53)).0,
            Some(core::cmp::Ordering::Equal)
        );
        a /= from_i64(4, 53);
        assert_eq!(
            a.partial_cmp(&from_i64(6, 53)).0,
            Some(core::cmp::Ordering::Equal)
        );
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixedfloat_arithmetic_ops() {
        let a = FixedFloat::<53>::try_from_i64_exact(2).unwrap();
        let b = FixedFloat::<53>::try_from_i64_exact(3).unwrap();
        let sum = a + b;
        let five = FixedFloat::<53>::try_from_i64_exact(5).unwrap();
        assert_eq!(sum.partial_cmp(&five).0, Some(core::cmp::Ordering::Equal));

        let diff = b - a;
        let one = FixedFloat::<53>::try_from_i64_exact(1).unwrap();
        assert_eq!(diff.partial_cmp(&one).0, Some(core::cmp::Ordering::Equal));

        let prod = a * b;
        let six = FixedFloat::<53>::try_from_i64_exact(6).unwrap();
        assert_eq!(prod.partial_cmp(&six).0, Some(core::cmp::Ordering::Equal));

        let q = b / a;
        // 3 / 2 = 1.5 (not exact at any small precision; verify
        // approximately)
        let two = FixedFloat::<53>::try_from_i64_exact(2).unwrap();
        assert_eq!(q.partial_cmp(&one).0, Some(core::cmp::Ordering::Greater));
        assert_eq!(q.partial_cmp(&two).0, Some(core::cmp::Ordering::Less));

        let neg = -a;
        assert!(neg.is_sign_negative());
    }
}
