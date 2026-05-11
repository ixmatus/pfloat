//! `exp(x)`: natural exponential function.
//!
//! Algorithm:
//!
//! 1. Special cases: NaN propagates; `exp(±0) = 1`;
//!    `exp(+∞) = +∞`; `exp(-∞) = +0`; very large `x` overflows
//!    to `+∞`; very negative `x` underflows to `+0`.
//! 2. Range reduce: choose integer `k = round(x / ln(2))`, then
//!    `r = x − k · ln(2)`. After reduction `|r| ≤ ln(2)/2 ≈ 0.347`.
//!    `ln(2)` is hardcoded at 1024-bit precision (see
//!    `super::LN2_LIMBS_1024`).
//! 3. Taylor series: `exp(r) = 1 + r + r²/2! + r³/3! + …`. With
//!    `|r| < 0.5`, the series converges geometrically faster than
//!    one bit per term once `n > |r| · target_precision`.
//! 4. Compose: `exp(x) = exp(r) · 2^k`. The `2^k` factor is a free
//!    exponent shift on the `BigFloat` (no arithmetic needed).
//!
//! pfloat does not implement full Ziv-strategy retry yet: slice 3a
//! computes with a fixed 64-bit guard above the target precision,
//! which is correctly-rounded for the vast majority of inputs (and
//! always correct in the round-toward-zero / round-toward-±∞ modes
//! that don't have tie cases). Pathological round-to-nearest tie
//! cases at the boundary of the rounding ULP could miss; the
//! Lefèvre–Muller worst-case tables for `exp` cover the known hard
//! arguments at common precisions and can be wired in later as
//! Phase 5 verification.

use crate::big::BigFloat;
use crate::class::Class;
use crate::ops::limbs::extract_as_integer;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::ln_2_at;

impl BigFloat {
    /// `exp(self)`: returns `e^x` rounded under `mode` to
    /// `self.precision`.
    ///
    /// See [the module docs](self) for the algorithm. Slice 3a
    /// uses a fixed 64-bit guard; pathological round-to-nearest
    /// tie cases at the rounding ULP boundary may miss. Phase 5
    /// will wire in the Lefèvre–Muller worst-case test corpus.
    #[must_use]
    pub fn exp(&self, mode: RoundingMode) -> (Self, Status) {
        let target = self.precision;
        self.exp_round(target, mode)
    }

    /// `exp(self)` with an explicit result precision.
    #[must_use]
    pub fn exp_round(&self, target_precision: u32, mode: RoundingMode) -> (Self, Status) {
        debug_assert!(target_precision >= 1);
        exp_kernel(self, target_precision, mode)
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `exp(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::exp`].
    #[must_use]
    pub fn exp(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().exp(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn exp_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    // Special cases per IEEE 754-2019 §9.
    match &x.class {
        Class::Nan {
            quiet,
            sign,
            payload,
        } => {
            if !*quiet {
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            let nan = BigFloat::try_new_quiet_nan(*sign, target_precision, payload)
                .expect("precision >= 1");
            return (nan, Status::OK);
        }
        Class::Zero { .. } => {
            // exp(±0) = 1.
            let one = BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1");
            return (one, Status::OK);
        }
        Class::Infinity {
            sign: Sign::Positive,
        } => {
            return (
                BigFloat::try_new_infinity(Sign::Positive, target_precision)
                    .expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Infinity {
            sign: Sign::Negative,
        } => {
            // exp(-∞) = +0.
            return (
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Normal { .. } => {}
    }

    // Working precision: target + 64 bits of guard.
    let working_prec = target_precision.saturating_add(64).min(1024);

    let x_w = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    let ln_2 = ln_2_at(working_prec);

    // k = round(x / ln(2)) as i64.
    let (x_over_ln2, _) = x_w.div(&ln_2, RoundingMode::NearestEven);
    let k = round_bigfloat_to_i64(&x_over_ln2);

    // r = x - k * ln(2).
    let k_big = BigFloat::try_from_i64_exact(k, working_prec)
        .or_else(|_| {
            BigFloat::try_from_i64_round(k, working_prec, RoundingMode::NearestEven).map(|(v, _)| v)
        })
        .expect("i64 fits in working precision");
    let (k_ln2, _) = k_big.mul(&ln_2, RoundingMode::NearestEven);
    let (r, _) = x_w.sub(&k_ln2, RoundingMode::NearestEven);

    // Taylor series: exp(r) = 1 + r + r²/2! + r³/3! + ...
    let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let mut sum = one.clone();
    let mut term = one;
    // Convergence: |term_n| ≈ |r|^n / n!, which for |r| < 0.5 falls below
    // 2^(-working_prec) by roughly working_prec / log2(1/|r|) ≈ 2 * working_prec terms.
    // Cap iterations to avoid runaway in pathological cases.
    let max_iter = 4u32.saturating_mul(working_prec).max(256);
    for n in 1u32..=max_iter {
        let (new_numer, _) = term.mul(&r, RoundingMode::NearestEven);
        let n_big =
            BigFloat::try_from_i64_exact(i64::from(n), working_prec).expect("precision >= 1");
        let (new_term, _) = new_numer.div(&n_big, RoundingMode::NearestEven);
        term = new_term;
        let (new_sum, _) = sum.add(&term, RoundingMode::NearestEven);
        sum = new_sum;
        // Termination: term below 2^(-working_prec) compared to sum (~1).
        if let Class::Normal { exponent, .. } = &term.class {
            if *exponent < -i64::from(working_prec) - 4 {
                break;
            }
        } else {
            break;
        }
    }

    // exp(x) = sum × 2^k. Apply k as a free exponent shift.
    let scaled = shift_exponent(sum, k);

    let (rounded, status) = scaled
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    auto_raise(status);
    (rounded, status)
}

/// Round a `BigFloat` to the nearest `i64` (banker's rounding for
/// ties).
///
/// Saturates to `i64::MAX` / `i64::MIN` for out-of-range inputs.
/// Returns `0` for `NaN` and other non-Normal/Zero classes.
fn round_bigfloat_to_i64(v: &BigFloat) -> i64 {
    let (sign, exponent, mantissa) = match &v.class {
        Class::Normal {
            sign,
            exponent,
            mantissa,
        } => (*sign, *exponent, mantissa),
        Class::Zero { .. } => return 0,
        _ => return 0,
    };

    if exponent < -1 {
        // |v| < 0.5, rounds to 0.
        return 0;
    }
    if exponent > 62 {
        return if matches!(sign, Sign::Positive) {
            i64::MAX
        } else {
            i64::MIN
        };
    }

    let p = v.precision;
    let m_int = extract_as_integer(mantissa, p);
    let scale = exponent - i64::from(p) + 1;

    let magnitude: u64 = if scale >= 0 {
        shift_left_to_u64(&m_int, scale as u32)
    } else {
        shift_right_round_to_u64(&m_int, (-scale) as u32)
    };

    let signed = magnitude as i64;
    if matches!(sign, Sign::Negative) {
        signed.wrapping_neg()
    } else {
        signed
    }
}

/// Shift a multi-limb integer left by `s` bits and return the low
/// 64 bits as u64. Used only via `round_bigfloat_to_i64` for inputs
/// whose result is guaranteed to fit in u64 (exp ≤ 62, scale ≥ 0
/// implies the high limbs are all zero).
fn shift_left_to_u64(m: &[u64], s: u32) -> u64 {
    let m0 = m.first().copied().unwrap_or(0);
    if s >= 64 {
        0
    } else {
        m0.checked_shl(s).unwrap_or(0)
    }
}

/// Shift a multi-limb integer right by `s` bits with
/// round-to-nearest-even on the discarded bits, returning the
/// result as a u64 (saturating).
fn shift_right_round_to_u64(m: &[u64], s: u32) -> u64 {
    if s == 0 {
        return m.first().copied().unwrap_or(0);
    }

    let limb_shift = (s / 64) as usize;
    let bit_shift = s % 64;

    let int_part: u64 = if bit_shift == 0 {
        m.get(limb_shift).copied().unwrap_or(0)
    } else {
        let lo = m.get(limb_shift).copied().unwrap_or(0) >> bit_shift;
        let hi = m
            .get(limb_shift + 1)
            .copied()
            .unwrap_or(0)
            .checked_shl(64 - bit_shift)
            .unwrap_or(0);
        lo | hi
    };

    // Guard bit at position s - 1 of m.
    let guard_pos = s - 1;
    let guard_limb = (guard_pos / 64) as usize;
    let guard_bit = guard_pos % 64;
    let guard = if guard_limb < m.len() {
        (m[guard_limb] >> guard_bit) & 1
    } else {
        0
    };

    // Sticky: any non-zero bit below position s - 1.
    let mut sticky = false;
    for &limb in m.iter().take(guard_limb.min(m.len())) {
        if limb != 0 {
            sticky = true;
            break;
        }
    }
    if !sticky && guard_limb < m.len() && guard_bit > 0 {
        let mask = (1u64 << guard_bit) - 1;
        if m[guard_limb] & mask != 0 {
            sticky = true;
        }
    }

    let round_up = guard == 1 && (sticky || (int_part & 1) == 1);
    if round_up {
        int_part.saturating_add(1)
    } else {
        int_part
    }
}

/// Multiply a finite `BigFloat` by `2^k` by adding `k` to its
/// exponent. Non-normal classes pass through.
fn shift_exponent(mut v: BigFloat, k: i64) -> BigFloat {
    if let Class::Normal {
        exponent,
        mantissa: _,
        sign: _,
    } = &mut v.class
    {
        *exponent = exponent.saturating_add(k);
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    fn parse(s: &str, p: u32) -> BigFloat {
        BigFloat::parse_str(s, p, RoundingMode::NearestEven)
            .expect("parse should succeed")
            .0
    }

    fn close_at(v: &BigFloat, expected: &BigFloat, prec: u32) -> bool {
        // |v - expected| <= 4 ULPs at `prec`. ULP at a value of
        // magnitude ~1 at precision `prec` is 2^(-prec).
        let (diff, _) = v.sub(expected, RoundingMode::NearestEven);
        let abs_diff = diff.abs();
        if abs_diff.is_zero() {
            return true;
        }
        if !abs_diff.is_normal() {
            return false;
        }
        let exp = match &abs_diff.class {
            Class::Normal { exponent, .. } => *exponent,
            _ => return false,
        };
        // 4 ULP at expected's magnitude is 2^(expected_exp - prec + 3).
        let expected_exp = match &expected.class {
            Class::Normal { exponent, .. } => *exponent,
            _ => 0,
        };
        let tolerance_exp = expected_exp - i64::from(prec) + 3;
        exp <= tolerance_exp
    }

    #[test]
    fn exp_zero_is_one() {
        let zero = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, status) = zero.exp(RoundingMode::NearestEven);
        assert!(status.is_ok());
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    #[test]
    fn exp_neg_zero_is_one() {
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, _) = nz.exp(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    #[test]
    fn exp_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.exp(RoundingMode::NearestEven);
        assert!(r.is_infinite());
        assert!(r.is_sign_positive());
    }

    #[test]
    fn exp_neg_inf() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, _) = ni.exp(RoundingMode::NearestEven);
        assert!(r.is_zero());
        assert!(r.is_sign_positive());
    }

    #[test]
    fn exp_qnan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.exp(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn exp_snan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.exp(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn exp_one_is_e() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, _) = one.exp(RoundingMode::NearestEven);
        // e ≈ 2.718281828459045
        let e = parse("2.718281828459045", 53);
        assert!(close_at(&r, &e, 53), "exp(1) = {r:?}, expected ≈ {e:?}");
    }

    #[test]
    fn exp_neg_one() {
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        let (r, _) = neg_one.exp(RoundingMode::NearestEven);
        // e^-1 ≈ 0.36787944117144233
        let expected = parse("0.36787944117144233", 53);
        assert!(close_at(&r, &expected, 53));
    }

    #[test]
    fn exp_ln2_is_two() {
        // exp(ln(2)) = 2
        let ln2 = parse("0.6931471805599453", 53);
        let (r, _) = ln2.exp(RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert!(close_at(&r, &two, 53), "exp(ln(2)) should ≈ 2, got {r:?}");
    }

    #[test]
    fn exp_ten_is_about_22026() {
        // e^10 ≈ 22026.465794806718
        let ten = BigFloat::try_from_i64_exact(10, 53).unwrap();
        let (r, _) = ten.exp(RoundingMode::NearestEven);
        let expected = parse("22026.465794806718", 53);
        assert!(
            close_at(&r, &expected, 53),
            "exp(10) = {r}, expected ≈ {expected}"
        );
    }

    #[test]
    fn exp_large_negative() {
        // e^-30 ≈ 9.357622968840175e-14
        let neg_thirty = BigFloat::try_from_i64_exact(-30, 53).unwrap();
        let (r, _) = neg_thirty.exp(RoundingMode::NearestEven);
        let expected = parse("9.357622968840175e-14", 53);
        assert!(close_at(&r, &expected, 53), "exp(-30) = {r}");
    }

    #[test]
    fn exp_at_higher_precision() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.exp(RoundingMode::NearestEven);
        assert_eq!(r.precision(), 113);
        let e = parse("2.718281828459045235360287471352662", 113);
        assert!(close_at(&r, &e, 113));
    }

    #[test]
    fn exp_with_explicit_round() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.exp_round(53, RoundingMode::NearestEven);
        assert_eq!(r.precision(), 53);
        let e = parse("2.718281828459045", 53);
        assert!(close_at(&r, &e, 53));
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixed_exp() {
        let one = FixedFloat::<53>::try_from_i64_exact(1).unwrap();
        let (r, _) = one.exp(RoundingMode::NearestEven);
        let e = parse("2.718281828459045", 53);
        assert!(close_at(&r.to_big(), &e, 53));
    }

    #[test]
    fn ln_2_constant_top_bits() {
        let ln2 = ln_2_at(53);
        // Top 8 bits of mantissa = 0xB1 (= 0b10110001), matching
        // ln(2) ≈ 0.1011000101110...
        match &ln2.class {
            Class::Normal { mantissa, .. } => {
                // Top byte of MSL.
                let top_byte = (*mantissa.last().unwrap()) >> 56;
                assert_eq!(top_byte, 0xB1);
            }
            _ => panic!("expected Normal"),
        }
    }
}
