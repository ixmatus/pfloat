//! IEEE 754-2019 §9.2.1 cube root for [`BigFloat`].
//!
//! Algorithm: like [`super::sqrt`], `cbrt` is an exact-integer root, not
//! a Ziv-driven transcendental (ADR-0056). Extract `|mantissa|` as a
//! bottom-aligned integer, left-shift it by enough bits that the integer
//! cube root has at least `target_precision + GUARD` bits, take the
//! integer cube root via [`super::limbs::iroot_limbs`] with `k = 3`, and
//! route through the rounding pipeline with `pre_sticky = (remainder !=
//! 0)`. The shift amount makes the result-exponent split clean:
//! `(scale − L)` is a multiple of 3 — the cube-root generalization of
//! sqrt's even-parity shift.
//!
//! `cbrt` is an odd function, so the real sign is carried through:
//! `cbrt(−x) = −cbrt(|x|)`. A negative operand is in-domain (unlike
//! sqrt) and raises nothing.
//!
//! Like sqrt, cbrt's tie cases never arise: a perfect cube at the
//! shifted precision is exact (`remainder == 0`); otherwise the true
//! cube root is irrational and cannot lie on a half-way rounding
//! boundary. The rounding pipeline therefore handles all five modes via
//! the same `(intermediate_precision, pre_sticky)` interface used by
//! add/sub/mul/div/sqrt, with no Ziv loop and no calibration entry.
//!
//! Special cases per IEEE 754-2019 §6.2, §7.2, and §9.2.1:
//!
//! - NaN: propagated; sNaN raises `INVALID`.
//! - `cbrt(±0) = ±0` (sign preserved).
//! - `cbrt(±∞) = ±∞` (odd root preserves the sign of infinity).
//! - `cbrt(negative finite) = −cbrt(|x|)`: the real cube root.

use alloc::vec;
use alloc::vec::Vec;

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::mantissa::limbs_for;
use crate::rounding::{round_finite_to_precision, RoundingMode};
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::limbs::{extract_as_integer, iroot_limbs, or_left_shifted_into, top_set_bit};

impl BigFloat {
    /// IEEE 754-2019 `rootn(self, 3)`: the real cube root.
    ///
    /// Returns `cbrt(self)` rounded under `mode` to a precision of
    /// `self.precision`.
    #[must_use]
    pub fn cbrt(&self, mode: RoundingMode) -> (Self, Status) {
        self.cbrt_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `cbrt(self)` with an explicit result precision.
    pub fn cbrt_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(cbrt_kernel(self, target_precision, mode))
    }

    /// `cbrt` accumulating into a caller-supplied flag bag.
    #[must_use]
    pub fn cbrt_with_flags(&self, mode: RoundingMode, flags: &mut Status) -> Self {
        let (value, status) = self.cbrt(mode);
        *flags |= status;
        value
    }

    /// `cbrt_round` accumulating into a caller-supplied flag bag.
    pub fn cbrt_round_with_flags(
        &self,
        target_precision: u32,
        mode: RoundingMode,
        flags: &mut Status,
    ) -> Result<Self, BuildError> {
        let (value, status) = self.cbrt_round(target_precision, mode)?;
        *flags |= status;
        Ok(value)
    }
}

fn cbrt_kernel(a: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            // cbrt(±0) = ±0 per IEEE 754-2019 §6.3.
            let z = BigFloat::try_new_zero(*sign, target_precision)
                .expect("BigFloat invariant: precision >= 1");
            (z, Status::OK)
        }
        Class::Infinity { sign } => {
            // cbrt(±∞) = ±∞: the odd root preserves the sign of infinity.
            let inf = BigFloat::try_new_infinity(*sign, target_precision)
                .expect("BigFloat invariant: precision >= 1");
            (inf, Status::OK)
        }
        Class::Normal {
            sign,
            exponent,
            mantissa,
        } => cbrt_finite(
            *sign,
            *exponent,
            mantissa,
            a.precision,
            target_precision,
            mode,
        ),
    }
}

fn cbrt_finite(
    sign: Sign,
    e_a: i64,
    m_a: &[u64],
    p_a: u32,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    let m_a_int = extract_as_integer(m_a, p_a);
    // i128 so an exponent near i64::MIN/MAX (reachable under the no-emax
    // design) does not overflow the scale; the final result exponent
    // ~ e_a/3 is well within i64.
    let scale_a = i128::from(e_a) - i128::from(p_a) + 1;

    // Pick a shift L so:
    //   1. p_a + L >= 3 × (target_precision + guard) so the integer cube
    //      root has at least target_precision + guard bits.
    //   2. (scale_a - L) ≡ 0 (mod 3) so the result exponent splits
    //      cleanly: cbrt(v_a) = iroot(m_a × 2^L, 3) × 2^((scale_a-L)/3).
    let guard: u32 = 16;
    let mut l = (3u32.saturating_mul(target_precision.saturating_add(guard))).saturating_sub(p_a);
    // Drive (scale_a - L) to a multiple of 3. Increasing L only adds root
    // bits, so requirement 1 is preserved.
    let rem3 = (((scale_a - i128::from(l)) % 3 + 3) % 3) as u32;
    l += rem3;
    // Clamp at the SOURCE post-bump (pf-9wb2, ADR-0108; the div
    // rationale): the rem3 bump can push `p_a + l` one past the u32
    // domain at the ceiling, dropping the mantissa's top bits and
    // panicking the non-zero-root expect. The clamp may break the
    // mod-3 alignment only at the ceiling's edge, where the result
    // still rounds from a self-consistent buffer.
    l = l.min(u32::MAX - p_a);

    // Saturating (pf-9wb2, ADR-0107): wraps only with operands near
    // the documented u32::MAX precision ceiling (the allocation wall
    // itself); the saturated width keeps the buffer self-consistent.
    let total_bits = p_a.saturating_add(l);
    let shifted_limbs = limbs_for(total_bits);
    let mut shifted: Vec<u64> = vec![0u64; shifted_limbs];
    or_left_shifted_into(&mut shifted, &m_a_int, p_a, l);

    let (s, r) = iroot_limbs(&shifted, 3);

    let top_bit_s = top_set_bit(&s).expect("non-zero cube root of non-zero N");

    let intermediate_precision = (top_bit_s + 1) as u32;
    let intermediate_limbs = limbs_for(intermediate_precision);
    let mut intermediate: Vec<u64> = vec![0u64; intermediate_limbs];
    let dst_low_zero =
        ((intermediate_limbs as u64) * 64 - u64::from(intermediate_precision)) as u32;
    or_left_shifted_into(&mut intermediate, &s, intermediate_precision, dst_low_zero);

    let pre_sticky = r.iter().any(|&v| v != 0);

    // Result exponent:
    //   v_a = m_a_int × 2^scale_a = N × 2^(scale_a - L), N = m_a_int × 2^L
    //   cbrt(v_a) = cbrt(N) × 2^((scale_a - L) / 3)
    //   result_exp = top_bit_s + (scale_a - L) / 3
    let scale_diff = scale_a - i128::from(l);
    debug_assert!(
        scale_diff % 3 == 0,
        "scale_diff must be a multiple of 3 by construction"
    );
    let result_exp_i = i128::from(top_bit_s as i64) + scale_diff / 3;
    let result_exp =
        i64::try_from(result_exp_i).unwrap_or(if result_exp_i < 0 { i64::MIN } else { i64::MAX });

    let (value, status) = round_finite_to_precision(
        sign,
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
    fn cbrt_zero_preserves_sign() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, s) = pz.cbrt(RoundingMode::NearestEven);
        assert!(r.is_zero() && r.is_sign_positive() && s.is_ok());

        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r2, _) = nz.cbrt(RoundingMode::NearestEven);
        assert!(r2.is_zero() && r2.is_sign_negative());
    }

    #[test]
    fn cbrt_one_is_one() {
        let one = from_i64(1, 53);
        let (r, s) = one.cbrt(RoundingMode::NearestEven);
        assert!(s.is_ok());
        assert_eq_bf(&r, &one);
    }

    #[test]
    fn cbrt_eight_is_two() {
        let (r, s) = from_i64(8, 53).cbrt(RoundingMode::NearestEven);
        assert!(s.is_ok());
        assert_eq_bf(&r, &from_i64(2, 53));
    }

    #[test]
    fn cbrt_twentyseven_is_three() {
        let (r, s) = from_i64(27, 53).cbrt(RoundingMode::NearestEven);
        assert!(s.is_ok());
        assert_eq_bf(&r, &from_i64(3, 53));
    }

    #[test]
    fn cbrt_thousand_is_ten() {
        let (r, s) = from_i64(1000, 53).cbrt(RoundingMode::NearestEven);
        assert!(s.is_ok());
        assert_eq_bf(&r, &from_i64(10, 53));
    }

    #[test]
    fn cbrt_neg_eight_is_neg_two() {
        // The defining divergence from sqrt: the real cube root of a
        // negative is the negative real root, with no INVALID.
        let (r, s) = from_i64(-8, 53).cbrt(RoundingMode::NearestEven);
        assert!(s.is_ok());
        assert!(!s.invalid());
        assert_eq_bf(&r, &from_i64(-2, 53));
    }

    #[test]
    fn cbrt_neg_twentyseven_is_neg_three() {
        let (r, s) = from_i64(-27, 53).cbrt(RoundingMode::NearestEven);
        assert!(s.is_ok());
        assert_eq_bf(&r, &from_i64(-3, 53));
    }

    #[test]
    fn cbrt_two_is_inexact_and_between_one_and_two() {
        use core::cmp::Ordering;
        let (r, s) = from_i64(2, 53).cbrt(RoundingMode::NearestEven);
        assert!(s.inexact());
        // 2^(1/3) ≈ 1.2599
        assert_eq!(r.partial_cmp(&from_i64(1, 53)).0, Some(Ordering::Greater));
        assert_eq!(r.partial_cmp(&from_i64(2, 53)).0, Some(Ordering::Less));
    }

    #[test]
    fn cbrt_cubed_round_trip() {
        // cbrt(a)³ == a for a perfect cube; the inverse of mul·mul.
        let n = from_i64(27, 53);
        let (c, _) = n.cbrt(RoundingMode::NearestEven);
        let (sq, _) = c.mul(&c, RoundingMode::NearestEven);
        let (cube, _) = sq.mul(&c, RoundingMode::NearestEven);
        assert_eq_bf(&cube, &n);
    }

    #[test]
    fn cbrt_fraction_one_eighth_is_one_half() {
        // cbrt(1/8) = 1/2.
        let one = from_i64(1, 53);
        let eight = from_i64(8, 53);
        let (eighth, _) = one.div(&eight, RoundingMode::NearestEven);
        let (root, _) = eighth.cbrt(RoundingMode::NearestEven);
        let two = from_i64(2, 53);
        let (half, _) = one.div(&two, RoundingMode::NearestEven);
        assert_eq_bf(&root, &half);
    }

    #[test]
    fn cbrt_inf_preserves_sign() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, s) = pi.cbrt(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive() && s.is_ok());

        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r2, s2) = ni.cbrt(RoundingMode::NearestEven);
        assert!(r2.is_infinite() && r2.is_sign_negative() && s2.is_ok());
        assert!(!s2.invalid());
    }

    #[test]
    fn cbrt_qnan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Negative, 53, &[7]).unwrap();
        let (r, s) = q.cbrt(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan() && r.is_sign_negative() && !s.invalid());
    }

    #[test]
    fn cbrt_snan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, s) = sn.cbrt(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan() && s.invalid());
    }

    #[test]
    fn cbrt_round_rejects_zero_precision() {
        assert_eq!(
            from_i64(1, 53).cbrt_round(0, RoundingMode::NearestEven),
            Err(BuildError::PrecisionZero)
        );
    }

    #[test]
    fn cbrt_with_flags_accumulates() {
        let mut flags = Status::OK;
        let _ = from_i64(2, 53).cbrt_with_flags(RoundingMode::NearestEven, &mut flags);
        assert!(flags.inexact());
    }

    #[test]
    fn cbrt_higher_precision_perfect_cube() {
        // 5³ = 125 at 200 bits.
        let (r, _) = from_i64(125, 200).cbrt(RoundingMode::NearestEven);
        assert_eq!(r.precision(), 200);
        assert_eq_bf(&r, &from_i64(5, 200));
    }

    #[test]
    fn cbrt_is_odd() {
        // cbrt(-x) = -cbrt(x) for a non-perfect cube.
        let (a, _) = from_i64(5, 80).cbrt(RoundingMode::NearestEven);
        let (b, _) = from_i64(-5, 80).cbrt(RoundingMode::NearestEven);
        assert_eq_bf(&b, &a.negated());
    }

    #[test]
    fn cbrt_explicit_round_precision() {
        let (r, _) = from_i64(8, 113)
            .cbrt_round(53, RoundingMode::NearestEven)
            .unwrap();
        assert_eq!(r.precision(), 53);
        assert_eq_bf(&r, &from_i64(2, 53));
    }
}
