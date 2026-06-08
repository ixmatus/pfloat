//! Fused two-product sum and difference for [`BigFloat`]:
//! `mul_add_mul(a, b, c, d) = a·b + c·d` and
//! `mul_sub_mul(a, b, c, d) = a·b − c·d`, each correctly rounded with a
//! single rounding step.
//!
//! These are the cancellation-safe building blocks the complex multiply
//! and divide kernels need (`pfloat-complex`): a complex product
//! `(a + bi)(c + di)` has real part `a·c − b·d` and imaginary part
//! `a·d + b·c`, and complex division forms `a·c + b·d`, `b·c − a·d`, and
//! the denominator `c² + d²`. Every one is `a·b ± c·d`.
//!
//! # Why a single rounding is correct (no Ziv loop) — ADR-0088
//!
//! The exact product of a `p`-bit significand and a `q`-bit significand
//! has at most `p + q` bits, so [`mul_round`](BigFloat::mul_round) at
//! precision `c.precision + d.precision` performs no rounding: `c·d` is
//! computed exactly. Feeding that exact value as the addend to
//! [`fma`](BigFloat::fma) — whose contract is one rounding of the exact
//! `a·b + addend` — yields `round(a·b + c·d)`, correctly rounded. The
//! difference form negates the exact `c·d` first.
//!
//! Catastrophic cancellation (`a·b ≈ ∓c·d`) does not threaten this. The
//! exact sum is a single real value with a bounded representation (each
//! product is finite precision, so their exponent-aligned sum is exact in
//! the arbitrary-precision significand), and rounding an exact value once
//! is correctly rounded by definition. There is no accumulated error for a
//! guard band to bound, so no Ziv loop is needed. This is the structural
//! difference from [`hypot`](BigFloat::hypot), whose `sqrt` is irrational
//! and therefore genuinely requires the Ziv driver: there the sum of
//! squares is exact but its square root is not.
//!
//! Special-case behaviour (NaN, infinity, signed zero) composes from
//! `mul_round` (forming `c·d`) and `fma` (forming `a·b + addend`): each
//! follows IEEE 754-2019, so an `INVALID` from a `0 × ∞` product or a NaN
//! operand propagates into the returned `Status`.

use crate::big::{BigFloat, BuildError};
use crate::rounding::RoundingMode;
use crate::status::Status;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `a·b + c·d` (with `a = self`), correctly rounded under `mode` to a
    /// precision of `max(self, b, c, d)`.
    #[must_use]
    pub fn mul_add_mul(&self, b: &Self, c: &Self, d: &Self, mode: RoundingMode) -> (Self, Status) {
        let target = self
            .precision
            .max(b.precision)
            .max(c.precision)
            .max(d.precision);
        self.mul_add_mul_round(b, c, d, target, mode)
            .expect("max of four valid precisions is valid")
    }

    /// `a·b + c·d` with an explicit result precision.
    pub fn mul_add_mul_round(
        &self,
        b: &Self,
        c: &Self,
        d: &Self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(two_product_kernel(
            self,
            b,
            c,
            d,
            false,
            target_precision,
            mode,
        ))
    }

    /// `a·b + c·d` accumulating into a caller-supplied flag bag.
    #[must_use]
    pub fn mul_add_mul_with_flags(
        &self,
        b: &Self,
        c: &Self,
        d: &Self,
        mode: RoundingMode,
        flags: &mut Status,
    ) -> Self {
        let (value, status) = self.mul_add_mul(b, c, d, mode);
        *flags |= status;
        value
    }

    /// `a·b + c·d` with explicit precision, accumulating into a flag bag.
    pub fn mul_add_mul_round_with_flags(
        &self,
        b: &Self,
        c: &Self,
        d: &Self,
        target_precision: u32,
        mode: RoundingMode,
        flags: &mut Status,
    ) -> Result<Self, BuildError> {
        let (value, status) = self.mul_add_mul_round(b, c, d, target_precision, mode)?;
        *flags |= status;
        Ok(value)
    }

    /// `a·b − c·d` (with `a = self`), correctly rounded under `mode` to a
    /// precision of `max(self, b, c, d)`.
    #[must_use]
    pub fn mul_sub_mul(&self, b: &Self, c: &Self, d: &Self, mode: RoundingMode) -> (Self, Status) {
        let target = self
            .precision
            .max(b.precision)
            .max(c.precision)
            .max(d.precision);
        self.mul_sub_mul_round(b, c, d, target, mode)
            .expect("max of four valid precisions is valid")
    }

    /// `a·b − c·d` with an explicit result precision.
    pub fn mul_sub_mul_round(
        &self,
        b: &Self,
        c: &Self,
        d: &Self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(two_product_kernel(
            self,
            b,
            c,
            d,
            true,
            target_precision,
            mode,
        ))
    }

    /// `a·b − c·d` accumulating into a caller-supplied flag bag.
    #[must_use]
    pub fn mul_sub_mul_with_flags(
        &self,
        b: &Self,
        c: &Self,
        d: &Self,
        mode: RoundingMode,
        flags: &mut Status,
    ) -> Self {
        let (value, status) = self.mul_sub_mul(b, c, d, mode);
        *flags |= status;
        value
    }

    /// `a·b − c·d` with explicit precision, accumulating into a flag bag.
    pub fn mul_sub_mul_round_with_flags(
        &self,
        b: &Self,
        c: &Self,
        d: &Self,
        target_precision: u32,
        mode: RoundingMode,
        flags: &mut Status,
    ) -> Result<Self, BuildError> {
        let (value, status) = self.mul_sub_mul_round(b, c, d, target_precision, mode)?;
        *flags |= status;
        Ok(value)
    }
}

/// `a·b ± c·d`, one rounding. `subtract` selects the difference form.
///
/// `c·d` is formed exactly (precision `c.precision + d.precision` cannot
/// round a product of those operands), negated for the difference form,
/// and handed to `fma` as the addend, so the only rounding is `fma`'s
/// single rounding of the exact `a·b + addend`.
fn two_product_kernel(
    a: &BigFloat,
    b: &BigFloat,
    c: &BigFloat,
    d: &BigFloat,
    subtract: bool,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    // Exact `c·d`: a (p)-bit by (q)-bit product has at most p+q bits, so
    // this precision performs no rounding. `saturating_add` guards the
    // (unreachable for real precisions) u32 overflow; `.max(1)` keeps it a
    // valid precision even if both inputs were somehow precision 0 (the
    // BigFloat invariant forbids that, but the kernel stays total).
    let cd_precision = c.precision.saturating_add(d.precision).max(1);
    let (cd_exact, cd_status) = c
        .mul_round(d, cd_precision, RoundingMode::NearestEven)
        .expect("cd_precision >= 1");
    let addend = if subtract {
        cd_exact.negated()
    } else {
        cd_exact
    };
    // Single rounding: fma forms the exact `a·b`, adds the exact addend,
    // and rounds once to the target.
    let (value, fma_status) = a
        .fma_round(b, &addend, target_precision, mode)
        .expect("target_precision validated by the caller");
    // `cd_status` carries any INVALID from a `0 × ∞` product or a NaN
    // operand in `(c, d)`; the product is exact so it never adds INEXACT.
    (value, fma_status | cd_status)
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `a·b + c·d` for `FixedFloat`, delegating to [`BigFloat::mul_add_mul`].
    #[must_use]
    pub fn mul_add_mul(&self, b: &Self, c: &Self, d: &Self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self
            .to_big()
            .mul_add_mul(&b.to_big(), &c.to_big(), &d.to_big(), mode);
        (
            Self::try_from_big_exact(big).expect("precision matches PREC"),
            status,
        )
    }

    /// `a·b − c·d` for `FixedFloat`, delegating to [`BigFloat::mul_sub_mul`].
    #[must_use]
    pub fn mul_sub_mul(&self, b: &Self, c: &Self, d: &Self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self
            .to_big()
            .mul_sub_mul(&b.to_big(), &c.to_big(), &d.to_big(), mode);
        (
            Self::try_from_big_exact(big).expect("precision matches PREC"),
            status,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sign::Sign;
    use core::cmp::Ordering;

    fn bf(n: i64, p: u32) -> BigFloat {
        BigFloat::try_from_i64_exact(n, p).unwrap()
    }

    fn eq(a: &BigFloat, b: &BigFloat) -> bool {
        matches!(a.partial_cmp(b).0, Some(Ordering::Equal)) && a.precision() == b.precision()
    }

    #[test]
    fn add_form_basic() {
        // 2·3 + 4·5 = 26.
        let (r, s) = bf(2, 53).mul_add_mul(
            &bf(3, 53),
            &bf(4, 53),
            &bf(5, 53),
            RoundingMode::NearestEven,
        );
        assert!(s.is_ok());
        assert!(eq(&r, &bf(26, 53)));
    }

    #[test]
    fn sub_form_basic() {
        // 7·3 − 4·5 = 1.
        let (r, s) = bf(7, 53).mul_sub_mul(
            &bf(3, 53),
            &bf(4, 53),
            &bf(5, 53),
            RoundingMode::NearestEven,
        );
        assert!(s.is_ok());
        assert!(eq(&r, &bf(1, 53)));
    }

    #[test]
    fn sub_form_exact_zero_on_equal_products() {
        // a·b − c·d with a·b == c·d exactly is exactly 0 (no spurious
        // INEXACT): 6·4 − 8·3 = 24 − 24 = 0.
        let (r, s) = bf(6, 53).mul_sub_mul(
            &bf(4, 53),
            &bf(8, 53),
            &bf(3, 53),
            RoundingMode::NearestEven,
        );
        assert!(r.is_zero());
        assert!(!s.inexact(), "an exact cancellation must not raise INEXACT");
    }

    #[test]
    fn catastrophic_cancellation_is_exact_not_lossy() {
        // Two products that agree in their high bits but differ in a low
        // bit: the exact difference is tiny but exact, and one rounding
        // returns it without a Ziv loop. (2^40 + 1)·1 − 2^40·1 = 1.
        let big = bf(1i64 << 40, 80);
        let big_plus_one = bf((1i64 << 40) + 1, 80);
        let one = bf(1, 80);
        let (r, s) = big_plus_one.mul_sub_mul(&one, &big, &one, RoundingMode::NearestEven);
        assert!(eq(&r, &bf(1, 80)), "got {r:?}");
        assert!(!s.inexact());
    }

    #[test]
    fn single_rounding_beats_round_then_subtract() {
        // A target precision low enough that rounding each product first
        // would differ from rounding the exact difference once.
        // a·b = 33 (0b100001), c·d = 1; exact a·b − c·d = 32 = 0b100000,
        // exact at 3 bits. Round-then-subtract would round 33 -> 32 (3
        // bits), then 32 − 1 = 31 -> rounds to 32 as well here, so use a
        // case where the difference is what rounds: target 3 bits.
        let (one_round, _) = bf(33, 53)
            .mul_sub_mul_round(
                &bf(1, 53),
                &bf(1, 53),
                &bf(1, 53),
                3,
                RoundingMode::NearestEven,
            )
            .unwrap();
        // exact 33 − 1 = 32 = 0b100000 -> exact at 3 bits.
        assert!(eq(
            &one_round,
            &bf(32, 53)
                .round_to_precision(3, RoundingMode::NearestEven)
                .unwrap()
                .0
        ));
    }

    #[test]
    fn matches_separate_exact_then_round() {
        // For inputs whose products are exact and whose sum is exact at the
        // target, the fused result equals the obvious two-step computation.
        let a = bf(123, 200);
        let b = bf(456, 200);
        let c = bf(7, 200);
        let d = bf(89, 200);
        let (fused, _) = a.mul_add_mul(&b, &c, &d, RoundingMode::NearestEven);
        // 123·456 + 7·89 = 56088 + 623 = 56711.
        assert!(eq(&fused, &bf(56711, 200)));
    }

    #[test]
    fn nan_operand_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = bf(2, 53).mul_add_mul(&bf(3, 53), &q, &bf(5, 53), RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn zero_times_inf_product_is_invalid() {
        // c·d = 0 × ∞ is invalid; the flag must reach the result.
        let zero = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let inf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, s) = bf(2, 53).mul_add_mul(&bf(3, 53), &zero, &inf, RoundingMode::NearestEven);
        assert!(r.is_nan());
        assert!(s.invalid());
    }

    #[test]
    fn rejects_zero_precision() {
        assert_eq!(
            bf(1, 53).mul_add_mul_round(
                &bf(1, 53),
                &bf(1, 53),
                &bf(1, 53),
                0,
                RoundingMode::NearestEven
            ),
            Err(BuildError::PrecisionZero)
        );
    }
}
