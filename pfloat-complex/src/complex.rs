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

    /// The magnitude `|self| = hypot(re, im)`, a real [`RealScalar`],
    /// correctly rounded under `mode`. Delegates to the scalar `hypot`
    /// kernel (IEEE 754-2019 §9.2.1), so it is correctly rounded and
    /// inherits the infinity-dominates-NaN special cases; it is not the
    /// lossy `sqrt(self.norm_sqr())`.
    #[cfg(feature = "exp-log")]
    #[must_use]
    pub fn abs(&self, mode: RoundingMode) -> (T, Status) {
        self.re.hypot(&self.im, mode)
    }

    /// The squared magnitude `|self|² = re² + im²`, a real [`RealScalar`],
    /// correctly rounded under `mode` with a single rounding (the fused
    /// two-product `mul_add_mul`, ADR-0088). This is the squared norm, not
    /// `abs()²`; it can overflow to `+∞` where `abs` stays finite, and it is
    /// exact (no spurious `INEXACT`) when `re² + im²` is representable.
    #[must_use]
    pub fn norm_sqr(&self, mode: RoundingMode) -> (T, Status) {
        self.re.mul_add_mul(&self.re, &self.im, &self.im, mode)
    }

    /// The phase `arg(self) = atan2(im, re)` in `(−π, π]`, a real
    /// [`RealScalar`], correctly rounded under `mode`. The branch cut on the
    /// negative real axis and the signed-zero discrimination (`arg(−1 + 0i)
    /// = +π`, `arg(−1 − 0i) = −π`) follow C99 Annex G, carried by the scalar
    /// `atan2` kernel's IEEE 754-2019 §9.2.1 table.
    #[cfg(feature = "trig")]
    #[must_use]
    pub fn arg(&self, mode: RoundingMode) -> (T, Status) {
        self.im.atan2(&self.re, mode)
    }

    /// The polar form `(|self|, arg(self))`, each part correctly rounded
    /// under `mode`. The returned [`Status`] is the OR of the magnitude and
    /// phase statuses. (`trig` implies `exp-log`, so `abs` is available.)
    #[cfg(feature = "trig")]
    #[must_use]
    pub fn to_polar(&self, mode: RoundingMode) -> (T, T, Status) {
        let (r, s_r) = self.abs(mode);
        let (theta, s_t) = self.arg(mode);
        (r, theta, s_r | s_t)
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
        // Annex G §G.5.1 infinity recovery (ADR-0091): a complex infinity
        // times a value with a zero part collapses the fused cross products to
        // (NaN, NaN) (each carries a 0·∞ term). When that happens and an
        // operand is a complex infinity, recover the mandated infinity. The
        // recovery runs in BigFloat (the rare path only), bridged via
        // to_big/from_big; a genuine NaN with no infinity returns None and the
        // naive (NaN, NaN) stands.
        if re.is_nan() && im.is_nan() {
            let p = self
                .re
                .precision()
                .max(self.im.precision())
                .max(other.re.precision())
                .max(other.im.precision());
            if let Some((r, i, s)) = crate::specials::recover_mul(
                &self.re.to_big(),
                &self.im.to_big(),
                &other.re.to_big(),
                &other.im.to_big(),
                p,
            ) {
                return (
                    Self {
                        re: T::from_big(&r),
                        im: T::from_big(&i),
                    },
                    s,
                );
            }
        }
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

    #[cfg(feature = "exp-log")]
    #[test]
    fn abs_three_four_is_five() {
        // |3 + 4i| = hypot(3, 4) = 5, exact.
        let (r, s) = c(3, 4).abs(RoundingMode::NearestEven);
        assert!(s.is_ok());
        assert_eq!(r.compare(&bf(5)).0, Some(Ordering::Equal));
    }

    #[cfg(feature = "exp-log")]
    #[test]
    fn abs_matches_scalar_hypot_inexact() {
        // |1 + 1i| = hypot(1, 1) = √2, INEXACT, bit-for-bit the scalar hypot.
        let (r, s) = c(1, 1).abs(RoundingMode::NearestEven);
        let scalar = bf(1).hypot(&bf(1), RoundingMode::NearestEven).0;
        assert!(s.inexact());
        assert_eq!(r.compare(&scalar).0, Some(Ordering::Equal));
    }

    #[test]
    fn norm_sqr_is_exact_squared_norm() {
        // |3 + 4i|² = 9 + 16 = 25, exact (no spurious INEXACT), and distinct
        // from abs (it can overflow where abs would not).
        let (r, s) = c(3, 4).norm_sqr(RoundingMode::NearestEven);
        assert!(!s.inexact());
        assert_eq!(r.compare(&bf(25)).0, Some(Ordering::Equal));
    }

    #[cfg(feature = "trig")]
    #[test]
    fn arg_signed_zero_selects_the_branch() {
        // The Annex G branch cut on the negative real axis, carried by
        // atan2's signed-zero table: arg(−1 + 0i) = +π, arg(−1 − 0i) = −π.
        let neg1 = bf(-1);
        let pz = BigFloat::try_new_zero(pfloat::Sign::Positive, 64).unwrap();
        let nz = BigFloat::try_new_zero(pfloat::Sign::Negative, 64).unwrap();
        let (above, _) = Complex::new(neg1.clone(), pz).arg(RoundingMode::NearestEven);
        let (below, _) = Complex::new(neg1, nz).arg(RoundingMode::NearestEven);
        assert!(above.is_sign_positive(), "arg(−1 + 0i) = +π");
        assert!(below.is_sign_negative(), "arg(−1 − 0i) = −π");
        // Same magnitude (both π): the only difference is the branch sign.
        assert_eq!(
            above.abs().partial_cmp(&below.abs()).0,
            Some(Ordering::Equal)
        );
    }

    #[cfg(feature = "trig")]
    #[test]
    fn to_polar_of_i_is_one_and_quarter_turn() {
        // i = 0 + 1i: |i| = 1 exactly, arg(i) = π/2 (so 1 < θ < 2).
        let (r, theta, _) = c(0, 1).to_polar(RoundingMode::NearestEven);
        assert_eq!(r.compare(&bf(1)).0, Some(Ordering::Equal));
        assert!(theta.is_sign_positive());
        assert_eq!(theta.partial_cmp(&bf(1)).0, Some(Ordering::Greater));
        assert_eq!(theta.partial_cmp(&bf(2)).0, Some(Ordering::Less));
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
    fn div_finite_by_zero_is_complex_infinity() {
        // C4 (Annex G §G.5.1): a finite nonzero dividend over a complex-zero
        // divisor is a complex infinity, not the C3 componentwise NaN. For
        // (1 + 1i)/(0 + 0i) the directed infinity has sign from c = +0, so
        // both parts are +∞.
        let (r, _) = c(1, 1).div(&c(0, 0), RoundingMode::NearestEven);
        assert!(r.re.is_infinite() && r.re.is_sign_positive());
        assert!(r.im.is_infinite() && r.im.is_sign_positive());
    }

    #[test]
    fn div_zero_by_zero_is_nan() {
        // 0/0 stays NaN + INVALID (it reaches the §G.5.1 D1 branch and falls
        // out as (NaN, NaN) via ∞·0).
        let (r, s) = c(0, 0).div(&c(0, 0), RoundingMode::NearestEven);
        assert!(r.is_nan());
        assert!(s.invalid());
    }
}
