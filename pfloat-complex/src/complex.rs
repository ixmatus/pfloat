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

    /// `self · other`, the complex product `(a + bi)(c + di) =
    /// (a·c − b·d) + (a·d + b·c)i`. Each component is one fused two-product
    /// (the C1 primitive `mul_sub_mul` / `mul_add_mul`), correctly rounded
    /// with a single rounding (ADR-0088), so the product is componentwise
    /// correctly rounded with no Ziv loop. The returned [`Status`] is the OR
    /// of the two component statuses.
    #[must_use]
    pub fn mul(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        // a, b = self.re, self.im;  c, d = other.re, other.im.
        // re = a·c − b·d ;  im = a·d + b·c.
        let (re, s_re) = self.re.mul_sub_mul(&other.re, &self.im, &other.im, mode);
        let (im, s_im) = self.re.mul_add_mul(&other.im, &self.im, &other.re, mode);
        (Self { re, im }, s_re | s_im)
    }

    /// `self / other`, the complex quotient
    /// `(a + bi)/(c + di) = [(ac + bd) + (bc − ad)i] / (c² + d²)`,
    /// componentwise correctly rounded under `mode` (ADR-0090).
    ///
    /// Unlike multiply, the quotient is not correctly rounded by rounding
    /// each part's numerator and denominator separately. The kernel brackets
    /// the true numerators and denominator with their directed fused
    /// two-product pairs at a working precision above the output precision,
    /// forms the quotient interval, and grows the precision until both ends
    /// round to the same value (a Ziv loop, capped at five iterations). The
    /// working precision exceeds any `FixedFloat<PREC>`, so the kernel runs
    /// in `BigFloat` and bridges back through [`RealScalar::from_big`].
    #[must_use]
    pub fn div(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let p = self
            .re
            .precision()
            .max(self.im.precision())
            .max(other.re.precision())
            .max(other.im.precision());
        let (re, im, status) = crate::div::complex_div_big(
            &self.re.to_big(),
            &self.im.to_big(),
            &other.re.to_big(),
            &other.im.to_big(),
            p,
            mode,
        );
        (
            Self {
                re: T::from_big(&re),
                im: T::from_big(&im),
            },
            status,
        )
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

    #[test]
    fn mul_basic() {
        // (2 + 3i)(4 + 5i) = (8 − 15) + (10 + 12)i = −7 + 22i.
        let (r, s) = c(2, 3).mul(&c(4, 5), RoundingMode::NearestEven);
        assert!(s.is_ok());
        assert!(eq(&r, &c(-7, 22)));
    }

    #[test]
    fn i_squared_is_minus_one() {
        // (0 + 1i)^2 = −1 + 0i.
        let i = c(0, 1);
        let (r, _) = i.mul(&i, RoundingMode::NearestEven);
        assert!(eq(&r, &c(-1, 0)));
    }

    #[test]
    fn z_times_conj_is_norm_squared() {
        // (3 + 4i)(3 − 4i) = 9 + 16 = 25, with an exactly-zero imaginary
        // part (no spurious INEXACT on the cancelling component).
        let z = c(3, 4);
        let (r, s) = z.mul(&z.conj(), RoundingMode::NearestEven);
        assert!(eq(&r, &c(25, 0)));
        assert!(!s.inexact());
    }

    #[test]
    fn mul_is_commutative() {
        let z = c(2, -7);
        let w = c(-3, 5);
        let (zw, _) = z.mul(&w, RoundingMode::NearestEven);
        let (wz, _) = w.mul(&z, RoundingMode::NearestEven);
        assert!(eq(&zw, &wz));
    }

    #[test]
    fn div_by_real_is_componentwise() {
        // (6 + 8i)/(2) = 3 + 4i.
        let (r, _) = c(6, 8).div(&c(2, 0), RoundingMode::NearestEven);
        assert!(eq(&r, &c(3, 4)));
    }

    #[test]
    fn div_one_by_i_is_minus_i() {
        // 1/i = −i. re = (0+0)/1 = 0; im = (0−1)/1 = −1.
        let (r, _) = c(1, 0).div(&c(0, 1), RoundingMode::NearestEven);
        assert!(eq(&r, &c(0, -1)));
    }

    #[test]
    fn div_by_self_is_one_exactly() {
        let z = c(2, 3);
        let (r, s) = z.div(&z, RoundingMode::NearestEven);
        assert!(eq(&r, &c(1, 0)));
        assert!(!s.inexact(), "z/z = 1 exactly");
    }

    #[test]
    fn div_inverts_mul_exactly() {
        // (2 + 3i)(4 + 5i) = −7 + 22i, so (−7 + 22i)/(4 + 5i) = 2 + 3i,
        // an exact integer quotient.
        let z = c(-7, 22);
        let w = c(4, 5);
        let (r, s) = z.div(&w, RoundingMode::NearestEven);
        assert!(eq(&r, &c(2, 3)));
        assert!(!s.inexact());
    }

    #[test]
    fn div_round_trips_a_product() {
        // (z·w)/w == z, with z, w chosen so z·w is exact and z·w/w is too.
        let z = c(13, -29);
        let w = c(5, 8);
        let (zw, _) = z.mul(&w, RoundingMode::NearestEven);
        let (back, _) = zw.div(&w, RoundingMode::NearestEven);
        assert!(eq(&back, &z));
    }

    #[test]
    fn div_real_matches_scalar_div_bit_for_bit() {
        // (1 + 0i)/(3 + 0i) real part = 1/3, correctly rounded. It must equal
        // the scalar correctly-rounded 1/3 bit for bit, with INEXACT set —
        // the componentwise-CR claim reduced to the real axis.
        let (r, s) = c(1, 0).div(&c(3, 0), RoundingMode::NearestEven);
        let scalar_third = bf(1).div(&bf(3), RoundingMode::NearestEven).0;
        assert_eq!(r.re.compare(&scalar_third).0, Some(Ordering::Equal));
        assert!(r.im.is_zero());
        assert!(s.inexact());
    }

    #[test]
    fn div_real_matches_scalar_div_all_modes() {
        // The componentwise-CR claim under every rounding mode: (1)/(7) real
        // part must equal the scalar correctly-rounded 1/7 bit for bit in
        // each of the five modes (1/7 rounds differently per mode).
        use RoundingMode::{NearestAway, NearestEven, TowardNegative, TowardPositive, TowardZero};
        for mode in [
            NearestEven,
            NearestAway,
            TowardZero,
            TowardPositive,
            TowardNegative,
        ] {
            let (r, _) = c(1, 0).div(&c(7, 0), mode);
            let scalar = bf(1).div(&bf(7), mode).0;
            assert_eq!(
                r.re.compare(&scalar).0,
                Some(Ordering::Equal),
                "real part diverged from scalar 1/7 under {mode:?}"
            );
            assert!(r.im.is_zero());
        }
    }

    #[test]
    fn div_by_zero_is_nan_invalid_basic() {
        // Componentwise (C3): the numerator a·c + b·d is also zero when the
        // divisor c + di is zero, so each component is 0/0 = NaN + INVALID.
        // The C99 Annex G complex-infinity refinement is a later slice (C4).
        let (r, s) = c(1, 1).div(&c(0, 0), RoundingMode::NearestEven);
        assert!(r.is_nan());
        assert!(s.invalid());
    }
}
