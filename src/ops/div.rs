//! IEEE 754-2019 §6.5 division for [`BigFloat`].
//!
//! Division produces an infinite binary expansion in the general
//! case (e.g., `1/3 = 0.0101...`), so the kernel computes a
//! finite-precision quotient with guard bits and a remainder, then
//! routes through [`crate::rounding::round_finite_to_precision`]
//! with `pre_sticky = (remainder != 0)`.
//!
//! Special cases per IEEE 754-2019 §6.2, §7.2, and §7.3:
//!
//! - NaN operand: propagated via [`super::propagate_nan2`]; sNaN
//!   raises `INVALID`.
//! - `±0 / ±0`: invalid, returns qNaN + `INVALID`.
//! - `±∞ / ±∞`: invalid, returns qNaN + `INVALID`.
//! - `finite_nonzero / ±0`: returns ±∞ with the XOR-combined sign
//!   and raises `DIV_BY_ZERO` (§7.3).
//! - `±0 / finite_nonzero` (denominator non-zero non-NaN): returns
//!   signed zero with XOR-combined sign.
//! - `±∞ / finite`: returns ±∞ with XOR-combined sign.
//! - `finite / ±∞`: returns signed zero with XOR-combined sign.
//! - Two finite non-zero values: integer-mantissa long division
//!   via [`super::limbs::divmod_limbs`], then rounding.
//!
//! The mantissa-quotient algorithm is bit-by-bit long division
//! (the "schoolbook" path per the plan). Phase 7 may replace with
//! Knuth Algorithm D or Newton iteration for better asymptotic
//! complexity at large precisions.

use alloc::vec;
use alloc::vec::Vec;

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::mantissa::limbs_for;
use crate::rounding::{round_finite_to_precision, RoundingMode};
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::limbs::{divmod_limbs, extract_as_integer, or_left_shifted_into, top_set_bit};
use super::propagate_nan2;

impl BigFloat {
    /// IEEE 754-2019 `division(self, other)`.
    ///
    /// Returns `self / other` rounded under `mode` to a precision
    /// of `max(self.precision, other.precision)`.
    #[must_use]
    pub fn div(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let target = self.precision.max(other.precision);
        self.div_round(other, target, mode)
            .expect("max of two valid precisions is valid")
    }

    /// IEEE 754-2019 `division(self, other)` with explicit result
    /// precision.
    ///
    /// Returns [`BuildError::PrecisionZero`] when
    /// `target_precision == 0`.
    pub fn div_round(
        &self,
        other: &Self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(div_kernel(self, other, target_precision, mode))
    }

    /// `div` accumulating into a caller-supplied flag bag.
    #[must_use]
    pub fn div_with_flags(&self, other: &Self, mode: RoundingMode, flags: &mut Status) -> Self {
        let (value, status) = self.div(other, mode);
        *flags |= status;
        value
    }

    /// `div_round` accumulating into a caller-supplied flag bag.
    pub fn div_round_with_flags(
        &self,
        other: &Self,
        target_precision: u32,
        mode: RoundingMode,
        flags: &mut Status,
    ) -> Result<Self, BuildError> {
        let (value, status) = self.div_round(other, target_precision, mode)?;
        *flags |= status;
        Ok(value)
    }
}

fn div_kernel(
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
        // 0 / 0 and Inf / Inf: invalid.
        (Class::Zero { .. }, Class::Zero { .. })
        | (Class::Infinity { .. }, Class::Infinity { .. }) => {
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("BigFloat invariant: precision >= 1");
            auto_raise(Status::INVALID);
            (nan, Status::INVALID)
        }
        // ±Inf / anything finite (including 0): signed Inf, no flag.
        // IEEE 754-2019 §7.3 signals divideByZero only for a FINITE
        // dividend, so this must precede the (_, Zero) arm: Inf/0 is an
        // exact infinite result, not a §7.3 exception. (Inf/Inf and
        // Inf/NaN are already handled above.) Review 2026-05-29.
        (Class::Infinity { .. }, _) => {
            let inf = BigFloat::try_new_infinity(result_sign, target_precision)
                .expect("BigFloat invariant: precision >= 1");
            (inf, Status::OK)
        }
        // finite_nonzero / 0: ±Inf + DIV_BY_ZERO. Numerator non-NaN,
        // non-Inf, non-zero implies it is Normal.
        (_, Class::Zero { .. }) => {
            let inf = BigFloat::try_new_infinity(result_sign, target_precision)
                .expect("BigFloat invariant: precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            (inf, Status::DIV_BY_ZERO)
        }
        // 0 / finite_nonzero: signed zero with XOR-combined sign.
        (Class::Zero { .. }, _) => {
            let z = BigFloat::try_new_zero(result_sign, target_precision)
                .expect("BigFloat invariant: precision >= 1");
            (z, Status::OK)
        }
        // finite / Inf: signed zero with XOR-combined sign.
        (_, Class::Infinity { .. }) => {
            let z = BigFloat::try_new_zero(result_sign, target_precision)
                .expect("BigFloat invariant: precision >= 1");
            (z, Status::OK)
        }
        (Class::Normal { .. }, Class::Normal { .. }) => {
            div_finite_finite(a, b, result_sign, target_precision, mode)
        }
        _ => unreachable!("NaN already handled by propagate_nan2"),
    }
}

fn div_finite_finite(
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
        _ => unreachable!("div_finite_finite called with non-Normal a"),
    };
    let (e_b, m_b, p_b) = match &b.class {
        Class::Normal {
            exponent, mantissa, ..
        } => (*exponent, mantissa.as_slice(), b.precision),
        _ => unreachable!("div_finite_finite called with non-Normal b"),
    };

    // Extract m_a and m_b as bottom-aligned integers.
    let m_a_int = extract_as_integer(m_a, p_a);
    let m_b_int = extract_as_integer(m_b, p_b);

    // Shift the dividend left so the quotient lands at a useful
    // precision. We want the quotient to have at least
    // `target_precision + GUARD` bits so the rounding pipeline can
    // dispose of guard/sticky bits correctly.
    //
    // After left-shifting m_a by L bits and dividing by m_b, the
    // quotient has approximately (p_a + L - p_b) bits. To guarantee
    // at least `target_precision + GUARD` bits, set
    // L >= target_precision + GUARD + p_b - p_a + 1 (slack 1 for the
    // floor/ceil bounce).
    //
    // Use a generous L for safety; the divmod cost is dominated by
    // the bit-count of the dividend.
    let guard: u32 = 8;
    // Pick L large enough that the quotient is always at least
    // `target_precision + guard` bits. Worst case: m_a is at its
    // minimum (top bit only) and m_b is at its maximum (all bits
    // set), so the quotient is at the low end of its range and we
    // lose two bits relative to the simple estimate.
    let l = target_precision
        .saturating_add(guard)
        .saturating_add(p_b.saturating_sub(p_a).saturating_add(2));

    let total_bits = p_a + l;
    let shifted_limbs = limbs_for(total_bits);
    let mut shifted_dividend: Vec<u64> = vec![0u64; shifted_limbs];
    or_left_shifted_into(&mut shifted_dividend, &m_a_int, p_a, l);

    let (quotient, remainder) = divmod_limbs(&shifted_dividend, &m_b_int);

    // The quotient is non-zero because m_a is non-zero and L was
    // chosen to make the dividend bigger than m_b in magnitude.
    let top_bit = top_set_bit(&quotient).expect("non-zero quotient");

    // Build top-aligned intermediate at (top_bit + 1) bits.
    let intermediate_precision = (top_bit + 1) as u32;
    let intermediate_limbs = limbs_for(intermediate_precision);
    let mut intermediate: Vec<u64> = vec![0u64; intermediate_limbs];
    let dst_low_zero = (intermediate_limbs as u32) * 64 - intermediate_precision;
    or_left_shifted_into(
        &mut intermediate,
        &quotient,
        intermediate_precision,
        dst_low_zero,
    );

    // Sticky bit: any non-zero limb in the remainder means we lost
    // bits below the quotient's LSB.
    let pre_sticky = remainder.iter().any(|&l| l != 0);

    // Result exponent derivation:
    //   v_a / v_b = (m_a × 2^L / m_b) × 2^(scale_a - scale_b - L)
    //            = Q × 2^(scale_a - scale_b - L)
    // where Q has its MSB at integer position `top_bit`.
    // pfloat's exponent is the position of the result MSB:
    //   result_exp = top_bit + (scale_a - scale_b - L)
    // with scale_a = e_a - p_a + 1, scale_b = e_b - p_b + 1.
    // `e_a` and `e_b` can each approach the `i64` limits (operands
    // produced by, e.g., `exp` of a large argument), so the quotient
    // exponent `e_a − e_b ± …` can fall outside the `i64` range
    // pfloat uses for exponents. The bare `i64` arithmetic would
    // panic (debug-overflow) or wrap (release) — the same
    // caller-reachable defect fixed in `mul` (pf-rnc, fuzz-found via
    // Airy `bi_prime`, which composes both). Compute in `i128` and
    // saturate to the `i64` range, flagging `OVERFLOW`/`UNDERFLOW`
    // (the saturating contract `round_finite_to_precision` already
    // applies to a round-up past `i64::MAX`; pfloat has no `emax`,
    // so a saturated exponent is a finite value, not `±∞`).
    let scale_diff =
        (i128::from(e_a) - i128::from(p_a) + 1) - (i128::from(e_b) - i128::from(p_b) + 1);
    let result_exp_wide = i128::from(top_bit as i64) + scale_diff - i128::from(l);
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
        pre_sticky,
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
    fn div_one_by_one() {
        let one = from_i64(1, 53);
        let (q, s) = one.div(&one, RoundingMode::NearestEven);
        assert!(s.is_ok());
        assert_eq_bf(&q, &one);
    }

    #[test]
    fn div_six_by_three() {
        let six = from_i64(6, 53);
        let three = from_i64(3, 53);
        let (q, s) = six.div(&three, RoundingMode::NearestEven);
        assert!(s.is_ok());
        let two = from_i64(2, 53);
        assert_eq_bf(&q, &two);
    }

    #[test]
    fn div_one_by_two_is_half() {
        let one = from_i64(1, 53);
        let two = from_i64(2, 53);
        let (q, s) = one.div(&two, RoundingMode::NearestEven);
        assert!(s.is_ok());
        // 0.5 has mantissa with top bit only set, exponent -1.
        // We can check this via comparison: 0.5 + 0.5 = 1.
        let (sum, _) = q.add(&q, RoundingMode::NearestEven);
        assert_eq_bf(&sum, &one);
    }

    #[test]
    fn div_sign_rule() {
        let six = from_i64(6, 53);
        let neg_three = from_i64(-3, 53);
        let (q, _) = six.div(&neg_three, RoundingMode::NearestEven);
        let neg_two = from_i64(-2, 53);
        assert_eq_bf(&q, &neg_two);
        let (q2, _) = neg_three.div(&neg_three, RoundingMode::NearestEven);
        let one = from_i64(1, 53);
        assert_eq_bf(&q2, &one);
    }

    #[test]
    fn div_inexact_one_third() {
        // 1/3 cannot be represented exactly in binary; INEXACT must
        // be raised.
        let one = from_i64(1, 53);
        let three = from_i64(3, 53);
        let (_q, s) = one.div(&three, RoundingMode::NearestEven);
        assert!(s.inexact());
    }

    #[test]
    fn div_by_zero_raises_flag_and_signs_correctly() {
        let one = from_i64(1, 53);
        let pos_zero = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let neg_zero = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (q1, s1) = one.div(&pos_zero, RoundingMode::NearestEven);
        assert!(q1.is_infinite());
        assert!(q1.is_sign_positive());
        assert!(s1.div_by_zero());

        let (q2, s2) = one.div(&neg_zero, RoundingMode::NearestEven);
        assert!(q2.is_infinite());
        assert!(q2.is_sign_negative());
        assert!(s2.div_by_zero());

        let neg_one = from_i64(-1, 53);
        let (q3, _) = neg_one.div(&pos_zero, RoundingMode::NearestEven);
        assert!(q3.is_infinite());
        assert!(q3.is_sign_negative());
    }

    #[test]
    fn zero_div_zero_is_qnan_invalid() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (q, s) = pz.div(&pz, RoundingMode::NearestEven);
        assert!(q.is_quiet_nan());
        assert!(s.invalid());
    }

    #[test]
    fn inf_div_inf_is_qnan_invalid() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (q, s) = pi.div(&pi, RoundingMode::NearestEven);
        assert!(q.is_quiet_nan());
        assert!(s.invalid());
    }

    #[test]
    fn zero_div_finite_is_signed_zero() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let two = from_i64(2, 53);
        let (q, s) = pz.div(&two, RoundingMode::NearestEven);
        assert!(q.is_zero());
        assert!(q.is_sign_positive());
        assert!(s.is_ok());

        let neg_two = from_i64(-2, 53);
        let (q2, _) = pz.div(&neg_two, RoundingMode::NearestEven);
        assert!(q2.is_zero());
        assert!(q2.is_sign_negative());
    }

    #[test]
    fn inf_div_finite_is_signed_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let neg_two = from_i64(-2, 53);
        let (q, _) = pi.div(&neg_two, RoundingMode::NearestEven);
        assert!(q.is_infinite());
        assert!(q.is_sign_negative());
    }

    #[test]
    fn finite_div_inf_is_signed_zero() {
        let two = from_i64(2, 53);
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (q, _) = two.div(&pi, RoundingMode::NearestEven);
        assert!(q.is_zero());
        assert!(q.is_sign_positive());

        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (q2, _) = two.div(&ni, RoundingMode::NearestEven);
        assert!(q2.is_zero());
        assert!(q2.is_sign_negative());
    }

    #[test]
    fn div_qnan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let two = from_i64(2, 53);
        let (r, s) = q.div(&two, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(!s.invalid());
    }

    #[test]
    fn div_snan_raises_invalid() {
        let s = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let two = from_i64(2, 53);
        let (r, st) = s.div(&two, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(st.invalid());
    }

    #[test]
    fn div_round_trip_with_mul() {
        // (a * b) / b = a when there is no rounding loss.
        let a = from_i64(7, 53);
        let b = from_i64(11, 53);
        let (p, _) = a.mul(&b, RoundingMode::NearestEven);
        let (back, _) = p.div(&b, RoundingMode::NearestEven);
        assert_eq_bf(&back, &a);
    }

    #[test]
    fn div_round_rejects_zero_precision() {
        let one = from_i64(1, 53);
        assert_eq!(
            one.div_round(&one, 0, RoundingMode::NearestEven),
            Err(BuildError::PrecisionZero)
        );
    }

    #[test]
    fn div_cross_precision_promotes() {
        let a = BigFloat::try_from_i64_exact(15, 53).unwrap();
        let b = BigFloat::try_from_i64_exact(3, 113).unwrap();
        let (q, _) = a.div(&b, RoundingMode::NearestEven);
        assert_eq!(q.precision(), 113);
        let five = BigFloat::try_from_i64_exact(5, 113).unwrap();
        assert_eq_bf(&q, &five);
    }

    #[test]
    fn div_with_flags_accumulates() {
        let one = from_i64(1, 53);
        let three = from_i64(3, 53);
        let mut flags = Status::OK;
        let _ = one.div_with_flags(&three, RoundingMode::NearestEven, &mut flags);
        assert!(flags.inexact());
    }

    #[test]
    fn div_extreme_exponent_saturates_not_panics() {
        // Regression (pf-rnc, fuzz-found via Airy bi_prime): the
        // quotient exponent (scale_diff) was computed in i64 and
        // overflowed once e_a − e_b passed i64::MAX. Square 2 until
        // the next square would saturate (so `big`'s exponent
        // exceeds i64::MAX/2), take its reciprocal (exponent ≈ −that),
        // then big / tiny has exponent ≈ 2·that > i64::MAX and must
        // saturate with OVERFLOW — never panic or yield NaN.
        let mut big = from_i64(2, 53);
        loop {
            let (sq, st) = big.mul(&big, RoundingMode::NearestEven);
            if st.overflow() {
                break;
            }
            big = sq;
        }
        let one = from_i64(1, 53);
        let (tiny, _) = one.div(&big, RoundingMode::NearestEven);
        let (q, st) = big.div(&tiny, RoundingMode::NearestEven);
        assert!(!q.is_nan(), "exponent saturation must not produce NaN");
        assert!(st.overflow(), "huge / tiny must flag exponent OVERFLOW");
    }
}
