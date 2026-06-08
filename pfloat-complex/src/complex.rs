//! [`Complex<T>`]: a complex number as a pair of pfloat scalar components.

use pfloat::{RoundingMode, Status};

use crate::scalar::RealScalar;

/// A complex number `re + im·i`, both components a [`RealScalar`].
///
/// The fields are public: a complex number carries no validity invariant
/// (any pair of real components, including NaN and infinity parts, denotes a
/// valid value), so there is nothing for accessors to protect. The arithmetic
/// is componentwise and correctly rounded under an explicit rounding mode,
/// returning the merged [`Status`] of the two component operations.
#[derive(Clone, Debug)]
pub struct Complex<T> {
    /// The real part.
    pub re: T,
    /// The imaginary part.
    pub im: T,
}

impl<T: RealScalar> Complex<T> {
    /// A complex number from its real and imaginary parts.
    pub fn new(re: T, im: T) -> Self {
        Self { re, im }
    }

    /// The real part.
    pub fn re(&self) -> &T {
        &self.re
    }

    /// The imaginary part.
    pub fn im(&self) -> &T {
        &self.im
    }

    /// `true` if either component is NaN.
    pub fn is_nan(&self) -> bool {
        self.re.is_nan() || self.im.is_nan()
    }

    /// Negation `−(re + im·i) = −re − im·i`. Exact (a sign-bit flip per
    /// component), so it has no rounding mode and no [`Status`].
    #[must_use]
    pub fn neg(&self) -> Self {
        Self {
            re: self.re.negated(),
            im: self.im.negated(),
        }
    }

    /// Complex conjugate `conj(re + im·i) = re − im·i`. Exact (the
    /// imaginary part's sign is flipped), so it has no rounding mode and no
    /// [`Status`].
    #[must_use]
    pub fn conj(&self) -> Self {
        Self {
            re: self.re.clone(),
            im: self.im.negated(),
        }
    }

    /// `self + other`, componentwise, each part correctly rounded under
    /// `mode`. The returned [`Status`] is the OR of the two component
    /// statuses (`INEXACT` if either part rounded, and so on).
    #[must_use]
    pub fn add(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let (re, s_re) = self.re.add(&other.re, mode);
        let (im, s_im) = self.im.add(&other.im, mode);
        (Self { re, im }, s_re | s_im)
    }

    /// `self − other`, componentwise, each part correctly rounded under
    /// `mode`. The returned [`Status`] is the OR of the two component
    /// statuses.
    #[must_use]
    pub fn sub(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let (re, s_re) = self.re.sub(&other.re, mode);
        let (im, s_im) = self.im.sub(&other.im, mode);
        (Self { re, im }, s_re | s_im)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;
    use pfloat::BigFloat;

    fn bf(n: i64) -> BigFloat {
        BigFloat::try_from_i64_exact(n, 64).unwrap()
    }

    fn c(re: i64, im: i64) -> Complex<BigFloat> {
        Complex::new(bf(re), bf(im))
    }

    fn eq(a: &Complex<BigFloat>, b: &Complex<BigFloat>) -> bool {
        a.re.compare(&b.re).0 == Some(Ordering::Equal)
            && a.im.compare(&b.im).0 == Some(Ordering::Equal)
    }

    #[test]
    fn add_is_componentwise() {
        // (2 + 3i) + (5 + 7i) = 7 + 10i.
        let (r, s) = c(2, 3).add(&c(5, 7), RoundingMode::NearestEven);
        assert!(s.is_ok());
        assert!(eq(&r, &c(7, 10)));
    }

    #[test]
    fn sub_is_componentwise() {
        // (5 + 7i) − (2 + 3i) = 3 + 4i.
        let (r, _) = c(5, 7).sub(&c(2, 3), RoundingMode::NearestEven);
        assert!(eq(&r, &c(3, 4)));
    }

    #[test]
    fn neg_flips_both_parts() {
        // −(2 + 3i) = −2 − 3i.
        assert!(eq(&c(2, 3).neg(), &c(-2, -3)));
    }

    #[test]
    fn conj_flips_only_the_imaginary_part() {
        // conj(2 + 3i) = 2 − 3i; conj is its own inverse.
        let z = c(2, 3);
        assert!(eq(&z.conj(), &c(2, -3)));
        assert!(eq(&z.conj().conj(), &z));
    }

    #[test]
    fn add_then_sub_round_trips() {
        let z = c(11, -4);
        let w = c(6, 9);
        let (sum, _) = z.add(&w, RoundingMode::NearestEven);
        let (back, _) = sum.sub(&w, RoundingMode::NearestEven);
        assert!(eq(&back, &z));
    }

    #[test]
    fn status_merges_across_components() {
        // A real-part add that rounds (1 + 2^-70 loses the tiny addend at
        // precision 64, since ulp(1) = 2^-63) with an exact imaginary part:
        // the merged status must be INEXACT, proving the OR across the two
        // component operations.
        let (tiny, _) = bf(1).scale_by_pow2(-70);
        let z = Complex::new(bf(1), bf(0));
        let w = Complex::new(tiny, bf(0));
        let (_, s) = z.add(&w, RoundingMode::NearestEven);
        assert!(s.inexact());
    }

    #[test]
    fn nan_part_makes_complex_nan() {
        let nan = BigFloat::try_new_quiet_nan(pfloat::Sign::Positive, 64, &[]).unwrap();
        let z = Complex::new(nan, bf(1));
        assert!(z.is_nan());
        assert!(!c(1, 2).is_nan());
    }
}
