//! IEEE 754-2019 §6.5 multiplication for [`BigFloat`].
//!
//! Multiplication is the cleanest of the four arithmetic operations
//! in algorithmic shape: the mantissa product is an exact integer
//! product of `m_a × m_b`, and the result exponent is determined by
//! the product's leading bit. The kernel routes the product through
//! [`crate::rounding::round_finite_to_precision`] to round to the
//! target precision.
//!
//! Special cases per IEEE 754-2019 §6.2 and §7.2:
//!
//! - NaN operand: propagated via [`super::propagate_nan2`]; sNaN
//!   raises `INVALID`.
//! - `±0 × ±∞` (either order): invalid, returns qNaN + `INVALID`.
//! - `±∞ × finite` (finite non-zero): returns `±∞` with the
//!   XOR-combined sign.
//! - `±∞ × ±∞`: returns `±∞` with XOR-combined sign.
//! - `±0 × finite`: returns signed zero with XOR-combined sign.
//! - Two finite non-zero values: schoolbook or Karatsuba
//!   multiplication of the mantissas, then rounding.
//!
//! The multi-precision multiplication dispatcher lives in
//! [`super::limbs::multiply_limbs`]: schoolbook for either operand
//! ≤ [`super::limbs::KARATSUBA_THRESHOLD`] limbs, Karatsuba above.
//! FFT (Schönhage-Strassen) is deferred to 1.x per ADR-0010.

use alloc::vec;
use alloc::vec::Vec;

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::mantissa::limbs_for;
use crate::rounding::{round_finite_to_precision, RoundingMode};
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::limbs::{extract_as_integer, multiply_limbs, or_left_shifted_into, top_set_bit};
use super::propagate_nan2;

impl BigFloat {
    /// IEEE 754-2019 `multiplication(self, other)`.
    ///
    /// Returns the product rounded under `mode` to a precision of
    /// `max(self.precision, other.precision)`.
    #[must_use]
    pub fn mul(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let target = self.precision.max(other.precision);
        self.mul_round(other, target, mode)
            .expect("max of two valid precisions is valid")
    }

    /// IEEE 754-2019 `multiplication(self, other)` with explicit
    /// result precision.
    ///
    /// Returns [`BuildError::PrecisionZero`] when
    /// `target_precision == 0`.
    pub fn mul_round(
        &self,
        other: &Self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(mul_kernel(self, other, target_precision, mode))
    }

    /// `mul` accumulating into a caller-supplied flag bag.
    #[must_use]
    pub fn mul_with_flags(&self, other: &Self, mode: RoundingMode, flags: &mut Status) -> Self {
        let (value, status) = self.mul(other, mode);
        *flags |= status;
        value
    }

    /// `mul_round` accumulating into a caller-supplied flag bag.
    pub fn mul_round_with_flags(
        &self,
        other: &Self,
        target_precision: u32,
        mode: RoundingMode,
        flags: &mut Status,
    ) -> Result<Self, BuildError> {
        let (value, status) = self.mul_round(other, target_precision, mode)?;
        *flags |= status;
        Ok(value)
    }
}

fn mul_kernel(
    a: &BigFloat,
    b: &BigFloat,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    if let Some(propagated) = propagate_nan2(a, b, target_precision) {
        return propagated;
    }

    let sign_a = a.sign();
    let sign_b = b.sign();
    let result_sign = sign_a.xor(sign_b);

    match (&a.class, &b.class) {
        // Inf × 0 (either order): invalid.
        (Class::Infinity { .. }, Class::Zero { .. })
        | (Class::Zero { .. }, Class::Infinity { .. }) => {
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("BigFloat invariant: precision >= 1");
            auto_raise(Status::INVALID);
            (nan, Status::INVALID)
        }
        // Inf × Inf, Inf × finite, finite × Inf: signed infinity.
        (Class::Infinity { .. }, _) | (_, Class::Infinity { .. }) => {
            let inf = BigFloat::try_new_infinity(result_sign, target_precision)
                .expect("BigFloat invariant: precision >= 1");
            (inf, Status::OK)
        }
        // 0 × anything (excluding Inf and NaN, both handled above):
        // signed zero with XOR-combined sign.
        (Class::Zero { .. }, _) | (_, Class::Zero { .. }) => {
            let z = BigFloat::try_new_zero(result_sign, target_precision)
                .expect("BigFloat invariant: precision >= 1");
            (z, Status::OK)
        }
        (Class::Normal { .. }, Class::Normal { .. }) => {
            mul_finite_finite(a, b, result_sign, target_precision, mode)
        }
        _ => unreachable!("NaN already handled by propagate_nan2"),
    }
}

/// Multiply two finite non-zero values via integer-mantissa product
/// and route through the rounding pipeline.
fn mul_finite_finite(
    a: &BigFloat,
    b: &BigFloat,
    result_sign: Sign,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    let (e_a, m_a, p_a) = match &a.class {
        Class::Normal {
            exponent, mantissa, ..
        } => (*exponent, mantissa.as_slice(), a.precision),
        _ => unreachable!("mul_finite_finite called with non-Normal a"),
    };
    let (e_b, m_b, p_b) = match &b.class {
        Class::Normal {
            exponent, mantissa, ..
        } => (*exponent, mantissa.as_slice(), b.precision),
        _ => unreachable!("mul_finite_finite called with non-Normal b"),
    };

    // Convert top-aligned mantissas to bottom-aligned integers, then
    // multiply.
    let a_int = extract_as_integer(m_a, p_a);
    let b_int = extract_as_integer(m_b, p_b);
    let product = multiply_limbs(&a_int, &b_int);

    // Find the top bit of the product. It must be set because both
    // operands have their top bit set, so the product is at least
    // 2^((p_a - 1) + (p_b - 1)) = 2^(p_a + p_b - 2), well above zero.
    let top_bit = top_set_bit(&product).expect("non-zero product");

    // Convert to the rounding-pipeline's intermediate frame:
    //   - intermediate_precision = top_bit + 1 (bits in the product's
    //     mantissa-as-integer)
    //   - intermediate storage is top-aligned, so we place the
    //     product's bits [0, top_bit] at the top of the storage.
    // The product of two operands near the documented u32::MAX
    // precision ceiling spans more bits than a u32 can name; the raw
    // cast wrapped silently (pf-9wb2, ADR-0107). Reaching this needs
    // ~2^31-bit operands (a half-gigabyte mantissa each) — at the
    // ceiling's edge the saturation below keeps the arithmetic
    // self-consistent, and the debug assertion documents the
    // envelope.
    debug_assert!(
        top_bit < u32::MAX as usize,
        "operand-precision sum exceeds the u32 ceiling (ADR-0002 edge)"
    );
    let intermediate_precision = u32::try_from(top_bit.saturating_add(1)).unwrap_or(u32::MAX);
    let intermediate_limbs = limbs_for(intermediate_precision);
    let mut intermediate: Vec<u64> = vec![0u64; intermediate_limbs];
    let dst_low_zero =
        ((intermediate_limbs as u64) * 64 - u64::from(intermediate_precision)) as u32;
    or_left_shifted_into(
        &mut intermediate,
        &product,
        intermediate_precision,
        dst_low_zero,
    );

    // Result exponent:
    //   v_a × v_b = m_a_int × m_b_int × 2^((e_a - p_a + 1) + (e_b - p_b + 1))
    //            = M_prod × 2^(e_a + e_b - p_a - p_b + 2)
    // where M_prod has its top bit at integer position `top_bit`.
    // The pfloat exponent of the result is the position of its
    // MSB, which is `top_bit + (e_a + e_b - p_a - p_b + 2)`.
    // `e_a` and `e_b` can each approach `i64::MAX` (operands
    // produced by, e.g., `exp` of a large argument), so the true
    // result exponent can fall outside the `i64` range pfloat uses
    // for exponents. The bare `i64` additions would panic
    // (debug-overflow) or wrap (release) on such operands — a
    // caller-reachable defect (fuzz-found via Airy `bi_prime`,
    // pre-existing since slice 1d). Compute in `i128` (the sum of a
    // small `top_bit` and a few `i64`s cannot overflow `i128`) and
    // saturate to the `i64` range, flagging `OVERFLOW`/`UNDERFLOW`.
    // This is the same saturating contract `round_finite_to_precision`
    // already applies when a round-up pushes the exponent past
    // `i64::MAX` (see `rounding.rs`): pfloat has no `emax`, so an
    // exponent of `i64::MAX`/`i64::MIN` is a saturated finite value,
    // not `±∞`.
    let result_exp_wide = i128::from(top_bit as i64) + i128::from(e_a) + i128::from(e_b)
        - i128::from(p_a)
        - i128::from(p_b)
        + 2;
    let mut exp_saturation = Status::OK;
    let result_exp = if result_exp_wide > i128::from(i64::MAX) {
        exp_saturation = Status::OVERFLOW;
        i64::MAX
    } else if result_exp_wide < i128::from(i64::MIN) {
        exp_saturation = Status::UNDERFLOW;
        i64::MIN
    } else {
        result_exp_wide as i64
    };

    let (value, round_status) = round_finite_to_precision(
        result_sign,
        result_exp,
        &intermediate,
        intermediate_precision,
        false, // no upstream sticky from alignment (mul is exact before rounding)
        target_precision,
        mode,
    );
    let status = round_status | exp_saturation;
    auto_raise(status);
    (value, status)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn from_i64(n: i64, p: u32) -> BigFloat {
        BigFloat::try_from_i64_exact(n, p).unwrap()
    }

    fn assert_eq_bf(a: &BigFloat, b: &BigFloat) {
        use core::cmp::Ordering;
        assert_eq!(
            a.partial_cmp(b).0,
            Some(Ordering::Equal),
            "expected equal: {a:?} vs {b:?}"
        );
        assert_eq!(a.precision(), b.precision());
    }

    #[test]
    fn mul_one_by_one() {
        let one = from_i64(1, 53);
        let (p, s) = one.mul(&one, RoundingMode::NearestEven);
        assert!(s.is_ok());
        assert_eq_bf(&p, &one);
    }

    #[test]
    fn mul_two_by_three() {
        let two = from_i64(2, 53);
        let three = from_i64(3, 53);
        let (p, s) = two.mul(&three, RoundingMode::NearestEven);
        assert!(s.is_ok());
        let six = from_i64(6, 53);
        assert_eq_bf(&p, &six);
    }

    #[test]
    fn mul_sign_rule() {
        let pa = from_i64(3, 53);
        let nb = from_i64(-4, 53);
        let (p, _) = pa.mul(&nb, RoundingMode::NearestEven);
        let expected = from_i64(-12, 53);
        assert_eq_bf(&p, &expected);
        let (p2, _) = nb.mul(&nb, RoundingMode::NearestEven);
        let sixteen = from_i64(16, 53);
        assert_eq_bf(&p2, &sixteen);
    }

    #[test]
    fn mul_by_zero_is_zero() {
        let a = from_i64(42, 53);
        let zero = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (p, s) = a.mul(&zero, RoundingMode::NearestEven);
        assert!(s.is_ok());
        assert!(p.is_zero());
        assert!(p.is_sign_positive());
    }

    #[test]
    fn mul_by_negative_zero() {
        // 5 × -0 = -0
        let five = from_i64(5, 53);
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (p, _) = five.mul(&nz, RoundingMode::NearestEven);
        assert!(p.is_zero());
        assert!(p.is_sign_negative());
    }

    #[test]
    fn mul_inf_by_zero_is_invalid() {
        let inf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let zero = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (p, s) = inf.mul(&zero, RoundingMode::NearestEven);
        assert!(p.is_quiet_nan());
        assert!(s.invalid());
        let (p2, s2) = zero.mul(&inf, RoundingMode::NearestEven);
        assert!(p2.is_quiet_nan());
        assert!(s2.invalid());
    }

    #[test]
    fn mul_inf_by_finite() {
        let inf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let neg_two = from_i64(-2, 53);
        let (p, _) = inf.mul(&neg_two, RoundingMode::NearestEven);
        assert!(p.is_infinite());
        assert!(p.is_sign_negative());
    }

    #[test]
    fn mul_inf_by_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (p1, _) = pi.mul(&pi, RoundingMode::NearestEven);
        assert!(p1.is_infinite());
        assert!(p1.is_sign_positive());
        let (p2, _) = pi.mul(&ni, RoundingMode::NearestEven);
        assert!(p2.is_infinite());
        assert!(p2.is_sign_negative());
    }

    #[test]
    fn mul_qnan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let two = from_i64(2, 53);
        let (p, s) = q.mul(&two, RoundingMode::NearestEven);
        assert!(p.is_quiet_nan());
        assert!(!s.invalid());
    }

    #[test]
    fn mul_snan_raises_invalid() {
        let s = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let two = from_i64(2, 53);
        let (p, st) = s.mul(&two, RoundingMode::NearestEven);
        assert!(p.is_quiet_nan());
        assert!(st.invalid());
    }

    #[test]
    fn mul_large_values_exact_at_high_precision() {
        // 2^30 × 2^30 = 2^60. Both fit comfortably in 53-bit
        // precision; product 2^60 also fits in 1 bit (top bit only).
        let a = from_i64(1 << 30, 53);
        let b = from_i64(1 << 30, 53);
        let (p, s) = a.mul(&b, RoundingMode::NearestEven);
        assert!(s.is_ok());
        let expected = from_i64(1i64 << 60, 53);
        assert_eq_bf(&p, &expected);
    }

    #[test]
    fn mul_inexact_when_product_exceeds_precision() {
        // 3 × 5 = 15 = 0b1111, 4 significant bits. At precision 3,
        // rounding required.
        let three = from_i64(3, 53);
        let five = from_i64(5, 53);
        let (_p, status) = three
            .mul_round(&five, 3, RoundingMode::NearestEven)
            .unwrap();
        assert!(status.inexact());
    }

    #[test]
    fn mul_cross_precision_promotes() {
        let a = BigFloat::try_from_i64_exact(7, 53).unwrap();
        let b = BigFloat::try_from_i64_exact(11, 113).unwrap();
        let (p, _) = a.mul(&b, RoundingMode::NearestEven);
        assert_eq!(p.precision(), 113);
        let expected = BigFloat::try_from_i64_exact(77, 113).unwrap();
        assert_eq_bf(&p, &expected);
    }

    #[test]
    fn mul_round_rejects_zero_precision() {
        let one = from_i64(1, 53);
        assert_eq!(
            one.mul_round(&one, 0, RoundingMode::NearestEven),
            Err(BuildError::PrecisionZero)
        );
    }

    #[test]
    fn mul_with_flags_accumulates() {
        let three = from_i64(3, 53);
        let five = from_i64(5, 53);
        let mut flags = Status::OK;
        let _ = three
            .mul_round_with_flags(&five, 3, RoundingMode::NearestEven, &mut flags)
            .unwrap();
        assert!(flags.inexact());
    }

    #[test]
    fn mul_extreme_exponent_saturates_not_panics() {
        // Regression (pf-rnc, fuzz-found via Airy bi_prime): the
        // result exponent `top_bit + e_a + e_b - p_a - p_b + 2` was
        // computed in `i64` and overflowed (debug-panic / release-
        // wrap) once `e_a + e_b` passed `i64::MAX`. Squaring 2
        // repeatedly doubles the exponent; within ~63 steps the true
        // product exponent exceeds the `i64` range and must saturate
        // with `OVERFLOW` to a finite value (pfloat has no `emax`,
        // so the exponent clamps to `i64::MAX`, never `±∞`), and it
        // must never panic or yield `NaN`.
        let mut x = from_i64(2, 53);
        let mut saw_overflow = false;
        for _ in 0..70 {
            let (sq, st) = x.mul(&x, RoundingMode::NearestEven);
            assert!(!sq.is_nan(), "exponent saturation must not produce NaN");
            if st.overflow() {
                saw_overflow = true;
                assert!(!sq.is_nan());
                break;
            }
            x = sq;
        }
        assert!(
            saw_overflow,
            "repeated squaring must hit exponent saturation (OVERFLOW), not panic"
        );
    }
}
