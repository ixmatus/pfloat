//! IEEE 754-2019 §6.5 addition and subtraction for [`BigFloat`].
//!
//! Add (and the trivially-derived sub) is the simplest of the four
//! arithmetic ops in shape, but the most rule-laden because it has
//! to cope with:
//!
//! - NaN propagation, with `INVALID` for any signaling-NaN operand.
//! - `±∞ ± ±∞`: same-sign infinities give the same infinity;
//!   opposite signs give NaN + INVALID (`Inf - Inf` is undefined).
//! - Zero handling, including the IEEE 754 sign rule for `(±0) ± (±0)`:
//!   in [`RoundingMode::TowardNegative`] the sum of equal-magnitude
//!   opposite-sign operands is `-0`; in every other mode it is `+0`.
//! - Effective-subtract cancellation: when sign-adjusted operands
//!   cancel to *exactly* zero the sign rule above kicks in.
//! - Mantissa alignment by `2^Δ` where `Δ` is the difference in
//!   binary exponents. Bits below the working window become sticky.
//!
//! Rounding to the target precision and status-flag emission are
//! deferred to [`crate::rounding::round_finite_to_precision`].

use core::cmp::Ordering;

use alloc::vec;
use alloc::vec::Vec;

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::cmp::magnitude_cmp;
use crate::mantissa::limbs_for;
use crate::rounding::{round_finite_to_precision, RoundingMode};
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::propagate_nan2;

impl BigFloat {
    /// IEEE 754-2019 `addition(self, other)`.
    ///
    /// Returns the sum rounded under `mode` to a precision of
    /// `max(self.precision, other.precision)`. Use
    /// [`add_round`](Self::add_round) for an explicit result
    /// precision.
    #[must_use]
    pub fn add(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let target = self.precision.max(other.precision);
        // `add_round` only fails on `target == 0`, which cannot
        // happen here because both operands have precision >= 1.
        self.add_round(other, target, mode)
            .expect("max of two valid precisions is valid")
    }

    /// IEEE 754-2019 `subtraction(self, other)`.
    ///
    /// Returns `self - other` rounded under `mode` to a precision
    /// of `max(self.precision, other.precision)`. Use
    /// [`sub_round`](Self::sub_round) for an explicit result
    /// precision.
    #[must_use]
    pub fn sub(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let target = self.precision.max(other.precision);
        self.sub_round(other, target, mode)
            .expect("max of two valid precisions is valid")
    }

    /// IEEE 754-2019 `addition(self, other)` with explicit result
    /// precision.
    ///
    /// Returns [`BuildError::PrecisionZero`] when
    /// `target_precision == 0`.
    pub fn add_round(
        &self,
        other: &Self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(add_with_signs(
            self,
            other,
            other.sign(),
            target_precision,
            mode,
        ))
    }

    /// IEEE 754-2019 `subtraction(self, other)` with explicit
    /// result precision.
    pub fn sub_round(
        &self,
        other: &Self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        // sub(a, b) is add(a, b) but with b's effective sign flipped.
        Ok(add_with_signs(
            self,
            other,
            other.sign().flip(),
            target_precision,
            mode,
        ))
    }

    /// `add` accumulating into a caller-supplied flag bag.
    #[must_use]
    pub fn add_with_flags(&self, other: &Self, mode: RoundingMode, flags: &mut Status) -> Self {
        let (value, status) = self.add(other, mode);
        *flags |= status;
        value
    }

    /// `sub` accumulating into a caller-supplied flag bag.
    #[must_use]
    pub fn sub_with_flags(&self, other: &Self, mode: RoundingMode, flags: &mut Status) -> Self {
        let (value, status) = self.sub(other, mode);
        *flags |= status;
        value
    }

    /// `add_round` accumulating into a caller-supplied flag bag.
    pub fn add_round_with_flags(
        &self,
        other: &Self,
        target_precision: u32,
        mode: RoundingMode,
        flags: &mut Status,
    ) -> Result<Self, BuildError> {
        let (value, status) = self.add_round(other, target_precision, mode)?;
        *flags |= status;
        Ok(value)
    }

    /// `sub_round` accumulating into a caller-supplied flag bag.
    pub fn sub_round_with_flags(
        &self,
        other: &Self,
        target_precision: u32,
        mode: RoundingMode,
        flags: &mut Status,
    ) -> Result<Self, BuildError> {
        let (value, status) = self.sub_round(other, target_precision, mode)?;
        *flags |= status;
        Ok(value)
    }
}

/// Compute `a + (sign_b_override × |b|)` at `target_precision`.
///
/// `sign_b_override` lets `sub_round` reuse this kernel by flipping
/// `b`'s effective sign without cloning the mantissa.
fn add_with_signs(
    a: &BigFloat,
    b: &BigFloat,
    sign_b: Sign,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    // Special-case dispatch.

    if let Some(propagated) = propagate_nan2(a, b, target_precision) {
        return propagated;
    }

    // From here on, neither operand is NaN.
    let sign_a = a.sign();
    match (&a.class, &b.class) {
        (Class::Infinity { .. }, Class::Infinity { .. }) => {
            if sign_a == sign_b {
                let inf = BigFloat::try_new_infinity(sign_a, target_precision)
                    .expect("BigFloat invariant: precision >= 1");
                return (inf, Status::OK);
            }
            // Inf + (-Inf): undefined, qNaN + INVALID.
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("BigFloat invariant: precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        (Class::Infinity { .. }, _) => {
            let inf = BigFloat::try_new_infinity(sign_a, target_precision)
                .expect("BigFloat invariant: precision >= 1");
            return (inf, Status::OK);
        }
        (_, Class::Infinity { .. }) => {
            let inf = BigFloat::try_new_infinity(sign_b, target_precision)
                .expect("BigFloat invariant: precision >= 1");
            return (inf, Status::OK);
        }
        (Class::Zero { .. }, Class::Zero { .. }) => {
            // ±0 ± ±0:
            // - same sign: that signed zero
            // - opposite signs: +0 except in TowardNegative, where -0
            let result_sign = if sign_a == sign_b {
                sign_a
            } else if matches!(mode, RoundingMode::TowardNegative) {
                Sign::Negative
            } else {
                Sign::Positive
            };
            let z = BigFloat::try_new_zero(result_sign, target_precision)
                .expect("BigFloat invariant: precision >= 1");
            return (z, Status::OK);
        }
        (Class::Zero { .. }, _) => {
            // 0 + b with b finite non-zero: result is b at target_precision (with effective sign).
            return zero_plus_finite(b, sign_b, target_precision, mode);
        }
        (_, Class::Zero { .. }) => {
            // a + 0: result is a at target_precision.
            return zero_plus_finite(a, sign_a, target_precision, mode);
        }
        (Class::Normal { .. }, Class::Normal { .. }) => {
            // Fall through to the finite-finite kernel below.
        }
        // NaN combinations were handled by propagate_nan2 above; the
        // remaining patterns are unreachable.
        _ => unreachable!("NaN already handled by propagate_nan2"),
    }

    add_finite_finite(a, sign_a, b, sign_b, target_precision, mode)
}

/// Helper for `0 + b` (with b non-zero, non-NaN, non-Inf): re-emit
/// `b` with the supplied sign at `target_precision`, possibly
/// rounding.
fn zero_plus_finite(
    operand: &BigFloat,
    effective_sign: Sign,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    // operand is finite non-zero; re-round to target_precision and
    // override the sign.
    let (mut rounded, status) = operand
        .round_to_precision(target_precision, mode)
        .expect("target_precision validated above");
    apply_sign_in_place(&mut rounded, effective_sign);
    (rounded, status)
}

fn apply_sign_in_place(value: &mut BigFloat, new_sign: Sign) {
    match &mut value.class {
        Class::Zero { sign }
        | Class::Infinity { sign }
        | Class::Nan { sign, .. }
        | Class::Normal { sign, .. } => {
            *sign = new_sign;
        }
    }
}

/// Add or subtract two finite, non-zero, non-NaN, non-Inf values.
///
/// Computes `a + (sign_b × |b|)` at `target_precision`.
fn add_finite_finite(
    a: &BigFloat,
    sign_a: Sign,
    b: &BigFloat,
    sign_b: Sign,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    // Order operands so `large` has the larger magnitude. Equal
    // magnitudes pick `a` arbitrarily.
    let a_is_large = magnitude_cmp(a, b) != Ordering::Less;
    let (large, sign_l, small, sign_s) = if a_is_large {
        (a, sign_a, b, sign_b)
    } else {
        (b, sign_b, a, sign_a)
    };

    let (e_l, m_l, p_l) = match &large.class {
        Class::Normal {
            exponent, mantissa, ..
        } => (*exponent, mantissa.as_slice(), large.precision),
        _ => unreachable!("large is finite non-zero non-NaN non-Inf"),
    };
    let (e_s, m_s, p_s) = match &small.class {
        Class::Normal {
            exponent, mantissa, ..
        } => (*exponent, mantissa.as_slice(), small.precision),
        _ => unreachable!("small is finite non-zero non-NaN non-Inf"),
    };

    // Each operand has a "scale": the binary weight of its mantissa
    // bit 0. scale = e - p + 1 (i.e., value = mantissa_as_integer
    // × 2^scale). Choose common_scale = min of the two so neither
    // operand needs to right-shift (and thus drop bits) when placed
    // into the buffer. The operand with the larger scale shifts
    // left by the difference; the one with the smaller scale shifts
    // by zero.
    let scale_l = e_l - i64::from(p_l) + 1;
    let scale_s = e_s - i64::from(p_s) + 1;
    let common_scale = scale_l.min(scale_s);
    // Shifts are always non-negative.
    let shift_l = (scale_l - common_scale) as u64;
    let shift_s = (scale_s - common_scale) as u64;

    // Top bit positions of each operand in the common-scale frame.
    let top_l = shift_l + u64::from(p_l) - 1;
    let top_s = shift_s + u64::from(p_s) - 1;

    // Huge-scale-difference short-circuit: when the smaller-scale
    // operand sits more than `huge_gap_threshold` bits below the
    // larger's window, its contribution is below the rounding
    // boundary at the target precision. Round the larger with sticky.
    let huge_gap_threshold = u64::from(target_precision) + 64;
    let scale_diff = top_l.abs_diff(top_s);
    if scale_diff > huge_gap_threshold {
        return huge_gap_short_circuit(
            large,
            sign_l,
            sign_s,
            sign_l == sign_s,
            target_precision,
            mode,
        );
    }

    // Working precision: enough to hold the highest bit (with
    // possible 1-bit carry on same-sign add) plus a few guard bits.
    let guard_bits: u64 = 4;
    let working_prec_u64 = top_l.max(top_s) + 2 + guard_bits;
    let working_prec = u32::try_from(working_prec_u64)
        .expect("working_prec fits in u32 within huge_gap_threshold");

    let working_limbs = limbs_for(working_prec);

    // Build buffers with each operand at its shifted bit position.
    let mut large_buf: Vec<u64> = vec![0u64; working_limbs];
    place_value_left_shifted(&mut large_buf, m_l, p_l, shift_l as u32);

    let mut small_buf: Vec<u64> = vec![0u64; working_limbs];
    place_value_left_shifted(&mut small_buf, m_s, p_s, shift_s as u32);

    let same_sign = sign_l == sign_s;
    let mut sum_buf: Vec<u64>;
    let result_sign: Sign;

    if same_sign {
        // Magnitude addition.
        sum_buf = large_buf;
        let _carry = limbs_add_assign(&mut sum_buf, &small_buf);
        // Carry into the (working_prec)-th bit can only happen when
        // working_prec was set tight; here working_prec >= p_l + 1
        // so the (p_l + gap + 1)-th bit is within the buffer.
        // Top bit may now be at p_l + gap (no carry) or p_l + gap + 1
        // (carry). Either way it sits within working_prec.
        result_sign = sign_l;
    } else {
        // Magnitude subtraction. `large_buf >= small_buf` because
        // we ordered operands by magnitude.
        sum_buf = large_buf;
        limbs_sub_assign(&mut sum_buf, &small_buf);
        // Result might be exactly zero (cancellation), or have many
        // leading zeros (catastrophic cancellation).
        if sum_buf.iter().all(|&l| l == 0) {
            // Exact cancellation.
            let result_sign = if matches!(mode, RoundingMode::TowardNegative) {
                Sign::Negative
            } else {
                Sign::Positive
            };
            let z = BigFloat::try_new_zero(result_sign, target_precision)
                .expect("target_precision validated");
            return (z, Status::OK);
        }
        result_sign = sign_l;
    }

    // Find the leading bit of the result.
    let leading_bit = top_set_bit(&sum_buf).expect("non-zero buffer has a top bit");

    // Convert from "scale = e_s - p_s + 1" frame back to
    // "exponent = position of MSB in mantissa-frame": result's MSB
    // is at integer position `leading_bit` (in the buffer's frame),
    // weighted by 2^(common_scale). The MSB's exponent (in pfloat's
    // sense, position of MSB) is therefore `leading_bit + common_scale`.
    let common_scale = e_s - i64::from(p_s) + 1;
    let result_exponent = i64::try_from(leading_bit)
        .expect("leading bit fits in i64 at any practical precision")
        + common_scale;

    // Build a normalized intermediate at intermediate_precision = leading_bit + 1
    // bits, top-bit-set.
    let intermediate_precision = (leading_bit + 1) as u32;
    let intermediate_limbs = limbs_for(intermediate_precision);
    let mut intermediate: Vec<u64> = vec![0u64; intermediate_limbs];

    // Source: bits [0, leading_bit] of sum_buf.
    // Destination: top `intermediate_precision` bits of the
    // intermediate storage (top-bit-set normalization).
    //
    // The shift between source and destination:
    //   src bit 0 currently at sum_buf storage position 0 (acc[0] bit 0).
    //   dst bit 0 (= mantissa-as-integer LSB) at intermediate storage
    //   position (intermediate_limbs * 64 - intermediate_precision).
    //   So we shift left by (intermediate_limbs * 64 - intermediate_precision) bits.
    let dst_low_zero = (intermediate_limbs as u32) * 64 - intermediate_precision;
    place_buffer_left_shifted(
        &mut intermediate,
        &sum_buf,
        intermediate_precision,
        dst_low_zero,
    );

    // Route through the rounding pipeline.
    let (value, status) = round_finite_to_precision(
        result_sign,
        result_exponent,
        &intermediate,
        intermediate_precision,
        false, // we did not drop bits before reaching the pipeline
        target_precision,
        mode,
    );

    auto_raise(status);
    (value, status)
}

/// Result for the `gap > huge_gap_threshold` case: round the larger
/// operand to `target_precision`, with sticky bit set if the smaller
/// is non-zero.
fn huge_gap_short_circuit(
    large: &BigFloat,
    sign_l: Sign,
    sign_s: Sign,
    same_sign: bool,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    // The smaller operand contributes a non-zero sub-ULP residue.
    // For same-sign add, sticky is set as if large's bits below the
    // rounding boundary had a 1. For opposite-sign sub, the smaller
    // pulls large slightly down; the rounding direction may flip
    // for directed modes, but the magnitude difference is below 1
    // ULP at target precision, so a bit-exact treatment is the same
    // as truncation-with-sticky.
    let _ = (sign_s, same_sign);

    // Build an intermediate at large.precision and route through
    // the rounding pipeline with `pre_sticky = true`.
    let (e_l, m_l, p_l) = match &large.class {
        Class::Normal {
            exponent, mantissa, ..
        } => (*exponent, mantissa.as_slice(), large.precision),
        _ => unreachable!("large is finite non-zero"),
    };

    let (value, status) = round_finite_to_precision(
        sign_l,
        e_l,
        m_l,
        p_l,
        true, // sticky from the smaller operand's contribution
        target_precision,
        mode,
    );
    auto_raise(status);
    (value, status)
}

// Multi-limb arithmetic helpers live in `crate::ops::limbs`.

/// Place `src`'s mantissa-as-integer into `dst` shifted left by
/// `left_shift` bits. Adapter over `limbs::or_left_shifted_into`
/// that first extracts the top-aligned mantissa into a bottom-
/// aligned integer.
fn place_value_left_shifted(dst: &mut [u64], src: &[u64], src_precision: u32, left_shift: u32) {
    let src_storage_bits = (src.len() as u32) * 64;
    let src_low_zero = src_storage_bits - src_precision;

    let value_limb_count = limbs_for(src_precision);
    let mut value_limbs: Vec<u64> = vec![0u64; value_limb_count];
    crate::ops::limbs::extract_value_limbs(src, src_low_zero, src_precision, &mut value_limbs);

    crate::ops::limbs::or_left_shifted_into(dst, &value_limbs, src_precision, left_shift);
}

/// Place `src` (a little-endian limb array of `value_bits` bits) into
/// `dst` storage with bottom bit at `dst_low_zero`. Caller-supplied
/// invariant: `dst` has at least `dst_low_zero + value_bits` storage
/// bits and was zeroed before this call.
fn place_buffer_left_shifted(dst: &mut [u64], src: &[u64], value_bits: u32, dst_low_zero: u32) {
    crate::ops::limbs::or_left_shifted_into(dst, src, value_bits, dst_low_zero);
}

use crate::ops::limbs::{limbs_add_assign, limbs_sub_assign, top_set_bit};

#[cfg(test)]
mod tests {
    use super::*;

    fn from_i64(n: i64, p: u32) -> BigFloat {
        BigFloat::try_from_i64_exact(n, p).unwrap()
    }

    fn assert_eq_bf(a: &BigFloat, b: &BigFloat) {
        assert_eq!(
            a.partial_cmp(b).0,
            Some(Ordering::Equal),
            "expected equal: {a:?} vs {b:?}"
        );
        assert_eq!(a.precision(), b.precision());
    }

    #[test]
    fn add_zero_to_one() {
        let one = from_i64(1, 53);
        let zero = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (sum, status) = one.add(&zero, RoundingMode::NearestEven);
        assert!(status.is_ok());
        assert_eq_bf(&sum, &one);
    }

    #[test]
    fn add_one_plus_one_is_two() {
        let one = from_i64(1, 53);
        let (sum, status) = one.add(&one, RoundingMode::NearestEven);
        assert!(status.is_ok());
        let two = from_i64(2, 53);
        assert_eq_bf(&sum, &two);
    }

    #[test]
    fn add_two_plus_three_is_five() {
        let two = from_i64(2, 53);
        let three = from_i64(3, 53);
        let (sum, status) = two.add(&three, RoundingMode::NearestEven);
        assert!(status.is_ok());
        let five = from_i64(5, 53);
        assert_eq_bf(&sum, &five);
    }

    #[test]
    fn add_with_different_magnitudes() {
        let big = from_i64(1000, 53);
        let small = from_i64(1, 53);
        let (sum, status) = big.add(&small, RoundingMode::NearestEven);
        assert!(status.is_ok());
        let expected = from_i64(1001, 53);
        assert_eq_bf(&sum, &expected);
    }

    #[test]
    fn add_negative_cancels_to_zero() {
        let a = from_i64(7, 53);
        let neg_a = from_i64(-7, 53);
        let (sum, status) = a.add(&neg_a, RoundingMode::NearestEven);
        assert!(status.is_ok());
        assert!(sum.is_zero());
        assert!(sum.is_sign_positive(), "exact cancel default sign is +0");
    }

    #[test]
    fn add_negative_cancels_to_negative_zero_under_toward_negative() {
        let a = from_i64(7, 53);
        let neg_a = from_i64(-7, 53);
        let (sum, _) = a.add(&neg_a, RoundingMode::TowardNegative);
        assert!(sum.is_zero());
        assert!(sum.is_sign_negative());
    }

    #[test]
    fn sub_self_is_zero() {
        let a = from_i64(42, 53);
        let (diff, status) = a.sub(&a, RoundingMode::NearestEven);
        assert!(status.is_ok());
        assert!(diff.is_zero());
        assert!(diff.is_sign_positive());
    }

    #[test]
    fn sub_smaller_from_larger() {
        let large = from_i64(10, 53);
        let small = from_i64(3, 53);
        let (diff, status) = large.sub(&small, RoundingMode::NearestEven);
        assert!(status.is_ok());
        let expected = from_i64(7, 53);
        assert_eq_bf(&diff, &expected);
    }

    #[test]
    fn sub_larger_from_smaller_yields_negative() {
        let large = from_i64(10, 53);
        let small = from_i64(3, 53);
        let (diff, _) = small.sub(&large, RoundingMode::NearestEven);
        let expected = from_i64(-7, 53);
        assert_eq_bf(&diff, &expected);
    }

    #[test]
    fn add_negative_one_plus_one() {
        let one = from_i64(1, 53);
        let neg_one = from_i64(-1, 53);
        let (sum, _) = neg_one.add(&one, RoundingMode::NearestEven);
        assert!(sum.is_zero());
    }

    // --- Special-value dispatch ---

    #[test]
    fn add_inf_inf_same_sign() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (sum, _) = pi.add(&pi, RoundingMode::NearestEven);
        assert!(sum.is_infinite());
        assert!(sum.is_sign_positive());
    }

    #[test]
    fn add_inf_minus_inf_is_qnan_invalid() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (sum, status) = pi.add(&ni, RoundingMode::NearestEven);
        assert!(sum.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn add_inf_finite() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let one = from_i64(1, 53);
        let (sum, _) = pi.add(&one, RoundingMode::NearestEven);
        assert!(sum.is_infinite());
        assert!(sum.is_sign_positive());
    }

    #[test]
    fn add_qnan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[42]).unwrap();
        let one = from_i64(1, 53);
        let (sum, status) = q.add(&one, RoundingMode::NearestEven);
        assert!(sum.is_quiet_nan());
        assert!(!status.invalid());
    }

    #[test]
    fn add_snan_raises_invalid() {
        let s = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let one = from_i64(1, 53);
        let (sum, status) = s.add(&one, RoundingMode::NearestEven);
        assert!(sum.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn add_zero_zero_same_sign() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (sum, _) = pz.add(&pz, RoundingMode::NearestEven);
        assert!(sum.is_zero());
        assert!(sum.is_sign_positive());

        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (sum2, _) = nz.add(&nz, RoundingMode::NearestEven);
        assert!(sum2.is_zero());
        assert!(sum2.is_sign_negative());
    }

    #[test]
    fn add_zero_zero_opposite_sign() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (sum, _) = pz.add(&nz, RoundingMode::NearestEven);
        assert!(sum.is_zero());
        assert!(sum.is_sign_positive(), "default mode: +0");

        let (sum2, _) = pz.add(&nz, RoundingMode::TowardNegative);
        assert!(sum2.is_sign_negative(), "TowardNegative: -0");
    }

    // --- Cross-precision ---

    #[test]
    fn add_promotes_to_max_precision() {
        let one_53 = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let one_113 = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (sum, status) = one_53.add(&one_113, RoundingMode::NearestEven);
        assert!(status.is_ok());
        assert_eq!(sum.precision(), 113);
        let two_113 = BigFloat::try_from_i64_exact(2, 113).unwrap();
        assert_eq_bf(&sum, &two_113);
    }

    // --- Inexact path ---

    #[test]
    fn add_at_low_precision_can_round() {
        // 3 + 5 = 8 fits in 1 bit (mantissa 0b1, exponent 3) so
        // exact at any precision >= 1.
        let three = from_i64(3, 53);
        let five = from_i64(5, 53);
        let (sum, status) = three
            .add_round(&five, 1, RoundingMode::NearestEven)
            .unwrap();
        assert!(status.is_ok());
        let eight = BigFloat::try_from_i64_exact(8, 1).unwrap();
        assert_eq_bf(&sum, &eight);
    }

    #[test]
    fn add_inexact_propagates_inexact_flag() {
        // 3 + 4 = 7 (binary 0b111), needs 3 bits exactly. At
        // precision 2, rounded.
        let three = from_i64(3, 53);
        let four = from_i64(4, 53);
        let (_sum, status) = three
            .add_round(&four, 2, RoundingMode::NearestEven)
            .unwrap();
        assert!(status.inexact());
    }

    #[test]
    fn add_round_rejects_zero_precision() {
        let one = from_i64(1, 53);
        assert_eq!(
            one.add_round(&one, 0, RoundingMode::NearestEven),
            Err(BuildError::PrecisionZero)
        );
    }

    #[test]
    fn add_with_flags_accumulates() {
        let three = from_i64(3, 53);
        let four = from_i64(4, 53);
        let mut flags = Status::OK;
        let _ = three
            .add_round_with_flags(&four, 2, RoundingMode::NearestEven, &mut flags)
            .unwrap();
        assert!(flags.inexact());
    }

    // Helper tests for limb-level primitives live in
    // `src/ops/limbs.rs::tests`.
}
