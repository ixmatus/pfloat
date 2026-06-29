//! IEEE 754-2019 §5.4.1 square root for [`BigFloat`].
//!
//! Algorithm: extract the mantissa as a bottom-aligned integer,
//! left-shift it by enough bits so that the resulting integer square
//! root has at least `target_precision + GUARD` bits, then take the
//! integer sqrt and route through the rounding pipeline with
//! `pre_sticky = (remainder != 0)`. The shift amount is chosen to
//! make the result-exponent split clean (same parity as the
//! operand's scale).
//!
//! Sqrt is rare among arithmetic operations in that its tie cases
//! never arise: if the operand's value is a perfect square at the
//! shifted precision, the sqrt is exact and rounding is trivial;
//! otherwise the true sqrt is irrational and cannot lie exactly on
//! a half-way rounding boundary. The rounding pipeline therefore
//! handles all modes via the same `(intermediate_precision,
//! pre_sticky)` interface used by add/sub/mul/div.
//!
//! Special cases per IEEE 754-2019 §6.2 and §7.2:
//!
//! - NaN: propagated; sNaN raises `INVALID`.
//! - `sqrt(±0)`: returns the input (preserves the sign of zero).
//! - `sqrt(+∞)`: returns `+∞`.
//! - `sqrt(-∞)` and `sqrt(negative finite)`: returns qNaN +
//!   `INVALID` (real square root of a negative is not defined).

use alloc::vec;
use alloc::vec::Vec;

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::mantissa::limbs_for;
use crate::rounding::{round_finite_to_precision, RoundingMode};
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::limbs::{extract_as_integer, isqrt_limbs, or_left_shifted_into, top_set_bit};

impl BigFloat {
    /// IEEE 754-2019 `squareRoot(self)`.
    ///
    /// Returns `sqrt(self)` rounded under `mode` to a precision of
    /// `self.precision`.
    #[must_use]
    pub fn sqrt(&self, mode: RoundingMode) -> (Self, Status) {
        self.sqrt_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// IEEE 754-2019 `squareRoot(self)` with explicit result
    /// precision.
    pub fn sqrt_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(sqrt_kernel(self, target_precision, mode))
    }

    /// `sqrt` accumulating into a caller-supplied flag bag.
    #[must_use]
    pub fn sqrt_with_flags(&self, mode: RoundingMode, flags: &mut Status) -> Self {
        let (value, status) = self.sqrt(mode);
        *flags |= status;
        value
    }

    /// `sqrt_round` accumulating into a caller-supplied flag bag.
    pub fn sqrt_round_with_flags(
        &self,
        target_precision: u32,
        mode: RoundingMode,
        flags: &mut Status,
    ) -> Result<Self, BuildError> {
        let (value, status) = self.sqrt_round(target_precision, mode)?;
        *flags |= status;
        Ok(value)
    }
}

fn sqrt_kernel(a: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    match &a.class {
        Class::Nan {
            quiet,
            sign,
            payload,
        } => {
            if *quiet {
                let nan = BigFloat::try_new_quiet_nan(*sign, target_precision, payload)
                    .expect("BigFloat invariant: precision >= 1");
                (nan, Status::OK)
            } else {
                // sNaN raises INVALID and propagates qNaN.
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("BigFloat invariant: precision >= 1");
                auto_raise(Status::INVALID);
                (nan, Status::INVALID)
            }
        }
        Class::Zero { sign } => {
            // sqrt(±0) = ±0 per IEEE 754-2019 §6.3.
            let z = BigFloat::try_new_zero(*sign, target_precision)
                .expect("BigFloat invariant: precision >= 1");
            (z, Status::OK)
        }
        Class::Infinity { sign } => match sign {
            Sign::Positive => {
                let inf = BigFloat::try_new_infinity(Sign::Positive, target_precision)
                    .expect("BigFloat invariant: precision >= 1");
                (inf, Status::OK)
            }
            Sign::Negative => {
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("BigFloat invariant: precision >= 1");
                auto_raise(Status::INVALID);
                (nan, Status::INVALID)
            }
        },
        Class::Normal {
            sign: Sign::Negative,
            ..
        } => {
            // sqrt of a negative real is not defined.
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("BigFloat invariant: precision >= 1");
            auto_raise(Status::INVALID);
            (nan, Status::INVALID)
        }
        Class::Normal {
            sign: Sign::Positive,
            exponent,
            mantissa,
        } => sqrt_finite_positive(*exponent, mantissa, a.precision, target_precision, mode),
    }
}

fn sqrt_finite_positive(
    e_a: i64,
    m_a: &[u64],
    p_a: u32,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    let m_a_int = extract_as_integer(m_a, p_a);
    // i128 so an exponent near i64::MIN (reachable by dividing into a
    // saturated-exponent value under the no-emax design) does not
    // overflow the scale; the final result exponent ~ e_a/2 is well
    // within i64. Review 2026-05-29.
    let scale_a = i128::from(e_a) - i128::from(p_a) + 1;

    // Pick a shift L so:
    //   1. p_a + L >= 2 × (target_precision + guard) so the integer
    //      sqrt has at least target_precision + guard bits.
    //   2. (scale_a - L) is even (so the result exponent splits
    //      cleanly: sqrt(v_a) = isqrt(m_a × 2^L) × 2^((scale_a-L)/2)).
    let guard: u32 = 16;
    let mut l = (2u32.saturating_mul(target_precision.saturating_add(guard))).saturating_sub(p_a);
    let parity = ((scale_a - i128::from(l)) & 1) != 0;
    if parity {
        l += 1;
    }
    // Source clamp for uniformity with div/cbrt (pf-9wb2, ADR-0108);
    // sqrt's 1-bit parity bump happens to fit limbs_for's rounding
    // slack, but the invariant should not rest on that coincidence.
    l = l.min(u32::MAX - p_a);

    // Saturating (pf-9wb2, ADR-0107): wraps only with operands near
    // the documented u32::MAX precision ceiling (the allocation wall
    // itself); the saturated width keeps the buffer self-consistent.
    let total_bits = p_a.saturating_add(l);
    let shifted_limbs = limbs_for(total_bits);
    let mut shifted: Vec<u64> = vec![0u64; shifted_limbs];
    or_left_shifted_into(&mut shifted, &m_a_int, p_a, l);

    let (s, r) = isqrt_limbs(&shifted);

    let top_bit_s = top_set_bit(&s).expect("non-zero sqrt of non-zero N");

    let intermediate_precision = (top_bit_s + 1) as u32;
    let intermediate_limbs = limbs_for(intermediate_precision);
    let mut intermediate: Vec<u64> = vec![0u64; intermediate_limbs];
    let dst_low_zero =
        ((intermediate_limbs as u64) * 64 - u64::from(intermediate_precision)) as u32;
    or_left_shifted_into(&mut intermediate, &s, intermediate_precision, dst_low_zero);

    let pre_sticky = r.iter().any(|&v| v != 0);

    // Result exponent:
    //   v_a = m_a_int × 2^scale_a = (m_a_int × 2^L) × 2^(scale_a - L) = N × 2^(scale_a - L)
    //   sqrt(v_a) = sqrt(N) × 2^((scale_a - L) / 2)
    //   In pfloat's exponent (position of MSB):
    //     result_exp = top_bit_s + (scale_a - L) / 2
    let scale_diff = scale_a - i128::from(l);
    debug_assert!(
        scale_diff & 1 == 0,
        "scale_diff must be even by construction"
    );
    let result_exp_i = i128::from(top_bit_s as i64) + scale_diff / 2;
    let result_exp =
        i64::try_from(result_exp_i).unwrap_or(if result_exp_i < 0 { i64::MIN } else { i64::MAX });

    let (value, status) = round_finite_to_precision(
        Sign::Positive,
        result_exp,
        &intermediate,
        intermediate_precision,
        pre_sticky,
        target_precision,
        mode,
    );
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
    fn sqrt_zero() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, s) = pz.sqrt(RoundingMode::NearestEven);
        assert!(r.is_zero());
        assert!(r.is_sign_positive());
        assert!(s.is_ok());

        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r2, _) = nz.sqrt(RoundingMode::NearestEven);
        // sqrt(-0) = -0 per IEEE 754-2019 §6.3.
        assert!(r2.is_zero());
        assert!(r2.is_sign_negative());
    }

    #[test]
    fn sqrt_one() {
        let one = from_i64(1, 53);
        let (r, s) = one.sqrt(RoundingMode::NearestEven);
        assert!(s.is_ok());
        assert_eq_bf(&r, &one);
    }

    #[test]
    fn sqrt_four_is_two() {
        let four = from_i64(4, 53);
        let (r, s) = four.sqrt(RoundingMode::NearestEven);
        assert!(s.is_ok());
        let two = from_i64(2, 53);
        assert_eq_bf(&r, &two);
    }

    #[test]
    fn sqrt_nine_is_three() {
        let nine = from_i64(9, 53);
        let (r, s) = nine.sqrt(RoundingMode::NearestEven);
        assert!(s.is_ok());
        let three = from_i64(3, 53);
        assert_eq_bf(&r, &three);
    }

    #[test]
    fn sqrt_sixteen_is_four() {
        let sixteen = from_i64(16, 53);
        let (r, s) = sixteen.sqrt(RoundingMode::NearestEven);
        assert!(s.is_ok());
        let four = from_i64(4, 53);
        assert_eq_bf(&r, &four);
    }

    #[test]
    fn sqrt_two_is_inexact() {
        let two = from_i64(2, 53);
        let (r, s) = two.sqrt(RoundingMode::NearestEven);
        assert!(s.inexact());
        // sqrt(2) ≈ 1.414... so result is in (1, 2).
        let one = from_i64(1, 53);
        let two_val = from_i64(2, 53);
        assert_eq!(r.partial_cmp(&one).0, Some(core::cmp::Ordering::Greater));
        assert_eq!(r.partial_cmp(&two_val).0, Some(core::cmp::Ordering::Less));
    }

    #[test]
    fn sqrt_squared_round_trip() {
        // sqrt(a)² == a when a is a perfect square; otherwise close.
        let nine = from_i64(9, 53);
        let (r, _) = nine.sqrt(RoundingMode::NearestEven);
        let (back, _) = r.mul(&r, RoundingMode::NearestEven);
        assert_eq_bf(&back, &nine);
    }

    #[test]
    fn sqrt_inf_pos() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, s) = pi.sqrt(RoundingMode::NearestEven);
        assert!(r.is_infinite());
        assert!(r.is_sign_positive());
        assert!(s.is_ok());
    }

    #[test]
    fn sqrt_inf_neg_is_invalid() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, s) = ni.sqrt(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(s.invalid());
    }

    #[test]
    fn sqrt_negative_finite_is_invalid() {
        let neg = from_i64(-4, 53);
        let (r, s) = neg.sqrt(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(s.invalid());
    }

    #[test]
    fn sqrt_qnan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Negative, 53, &[7]).unwrap();
        let (r, s) = q.sqrt(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(r.is_sign_negative());
        assert!(!s.invalid());
    }

    #[test]
    fn sqrt_snan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, s) = sn.sqrt(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(s.invalid());
    }

    #[test]
    fn sqrt_round_rejects_zero_precision() {
        let one = from_i64(1, 53);
        assert_eq!(
            one.sqrt_round(0, RoundingMode::NearestEven),
            Err(BuildError::PrecisionZero)
        );
    }

    #[test]
    fn sqrt_with_flags_accumulates() {
        let two = from_i64(2, 53);
        let mut flags = Status::OK;
        let _ = two.sqrt_with_flags(RoundingMode::NearestEven, &mut flags);
        assert!(flags.inexact());
    }

    #[test]
    fn sqrt_large_perfect_square() {
        // 100² = 10000.
        let ten_k = from_i64(10000, 53);
        let (r, _) = ten_k.sqrt(RoundingMode::NearestEven);
        let hundred = from_i64(100, 53);
        assert_eq_bf(&r, &hundred);
    }

    #[test]
    fn sqrt_fraction_half() {
        // sqrt(0.25) = 0.5. Compute 1/4 via division then sqrt.
        let one = from_i64(1, 53);
        let four = from_i64(4, 53);
        let (quarter, _) = one.div(&four, RoundingMode::NearestEven);
        let (half_from_sqrt, _) = quarter.sqrt(RoundingMode::NearestEven);
        // 0.5 = 1/2.
        let two = from_i64(2, 53);
        let (half, _) = one.div(&two, RoundingMode::NearestEven);
        assert_eq_bf(&half_from_sqrt, &half);
    }
}
