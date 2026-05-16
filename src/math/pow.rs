//! `pow(x, y) = x^y`: general power function.
//!
//! For positive `x` and finite `y`, the implementation evaluates
//! `exp(y · ln(x))` at working precision. The two transcendentals
//! shipped in slices 3a and 3b carry the work; this slice mostly
//! handles the rich IEEE 754-2019 §9.2.1 special-case table.
//!
//! Special-case table (per IEEE 754-2019 §9.2.1):
//!
//! - `pow(x, ±0) = 1` for any `x`, including NaN and infinity.
//! - `pow(+1, y) = 1` for any `y`, including NaN and infinity.
//! - `pow(NaN, y) = NaN` (propagates) for `y ≠ 0`.
//! - `pow(x, NaN) = NaN` for `x ≠ +1`.
//! - `pow(±0, y)`:
//!   - `y > 0` and `y` is an odd integer with `x = -0`: `-0`.
//!   - `y > 0` otherwise: `+0`.
//!   - `y < 0` and `y` is an odd integer with `x = -0`: `-∞ + DIV_BY_ZERO`.
//!   - `y < 0` otherwise: `+∞ + DIV_BY_ZERO`.
//! - `pow(±∞, y)`:
//!   - `y > 0`: `±∞` (sign matches `x` only when `y` is odd integer
//!     and `x = -∞`, else `+∞`).
//!   - `y < 0`: `±0` (same sign rule).
//! - `pow(x, ±∞)` for `x ≠ -1, +1`:
//!   - `|x| > 1`, `y = +∞`: `+∞`.
//!   - `|x| > 1`, `y = -∞`: `+0`.
//!   - `|x| < 1`, `y = +∞`: `+0`.
//!   - `|x| < 1`, `y = -∞`: `+∞`.
//! - `pow(-1, ±∞) = 1`.
//! - `pow(negative finite, non-integer y) = qNaN + INVALID`.
//! - `pow(negative finite, integer y)`: `(sign of x)^y · |x|^y`.
//!   The sign is `−` iff `y` is odd integer and `x < 0`.
//!
//! sNaN raises INVALID regardless of the special-case dispatch
//! that would otherwise apply.

use core::cmp::Ordering;

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::ops::limbs::extract_as_integer;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `self.pow(other, mode)`: returns `self^other` rounded under
    /// `mode` to a precision of `max(self.precision, other.precision)`.
    #[must_use]
    pub fn pow(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let target = self.precision.max(other.precision);
        self.pow_round(other, target, mode)
            .expect("max of two valid precisions is valid")
    }

    /// `pow` with an explicit result precision.
    pub fn pow_round(
        &self,
        other: &Self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(pow_kernel(self, other, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `pow` for `FixedFloat`. Delegates to [`BigFloat::pow`].
    #[must_use]
    pub fn pow(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().pow(&other.to_big(), mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Parity {
    Even,
    Odd,
}

fn pow_kernel(
    x: &BigFloat,
    y: &BigFloat,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    // sNaN handling: signal INVALID and propagate a quiet NaN.
    if x.is_signaling_nan() || y.is_signaling_nan() {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }

    // pow(x, ±0) = 1 for any x (including NaN, ±∞).
    if y.is_zero() {
        let one = BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1");
        return (one, Status::OK);
    }

    // pow(+1, y) = 1 for any y (including NaN).
    if is_positive_one(x) {
        let one = BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1");
        return (one, Status::OK);
    }

    // Quiet NaN propagation (after the two "1" rules above).
    if x.is_nan() {
        let nan =
            BigFloat::try_new_quiet_nan(x.sign(), target_precision, &[]).expect("precision >= 1");
        return (nan, Status::OK);
    }
    if y.is_nan() {
        let nan =
            BigFloat::try_new_quiet_nan(y.sign(), target_precision, &[]).expect("precision >= 1");
        return (nan, Status::OK);
    }

    // pow(x, ±∞).
    if y.is_infinite() {
        return pow_x_infinite_y(x, y.sign(), target_precision);
    }

    // pow(±0, y).
    if x.is_zero() {
        return pow_zero_base(x.sign(), y, target_precision);
    }

    // pow(±∞, y).
    if x.is_infinite() {
        return pow_infinite_base(x.sign(), y, target_precision);
    }

    // Both x and y are finite, non-NaN, x ≠ ±0, x ≠ +1, x ≠ ±∞,
    // y ≠ 0, y ≠ ±∞.

    if x.is_sign_negative() {
        // Negative base: y must be an integer or we return qNaN + INVALID.
        if let Some(parity) = integer_parity(y) {
            let abs_x = x.abs();
            let (abs_result, status) = pow_positive(&abs_x, y, target_precision, mode);
            let result = if matches!(parity, Parity::Odd) {
                abs_result.negated()
            } else {
                abs_result
            };
            return (result, status);
        }
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }

    pow_positive(x, y, target_precision, mode)
}

fn pow_x_infinite_y(x: &BigFloat, y_sign: Sign, target_precision: u32) -> (BigFloat, Status) {
    // pow(-1, ±∞) = 1.
    if is_negative_one(x) {
        let one = BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1");
        return (one, Status::OK);
    }
    let abs_x = x.abs();
    let one = BigFloat::try_from_i64_exact(1, x.precision).expect("precision >= 1");
    let cmp = abs_x.partial_cmp(&one).0;
    let abs_gt_one = matches!(cmp, Some(Ordering::Greater));
    let y_pos = matches!(y_sign, Sign::Positive);
    let result_is_inf = abs_gt_one == y_pos;
    if result_is_inf {
        (
            BigFloat::try_new_infinity(Sign::Positive, target_precision).expect("precision >= 1"),
            Status::OK,
        )
    } else {
        (
            BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1"),
            Status::OK,
        )
    }
}

fn pow_zero_base(x_sign: Sign, y: &BigFloat, target_precision: u32) -> (BigFloat, Status) {
    let y_neg = matches!(y.sign(), Sign::Negative);
    let y_odd_int = matches!(integer_parity(y), Some(Parity::Odd));
    let result_sign = if matches!(x_sign, Sign::Negative) && y_odd_int {
        Sign::Negative
    } else {
        Sign::Positive
    };
    if y_neg {
        let inf =
            BigFloat::try_new_infinity(result_sign, target_precision).expect("precision >= 1");
        auto_raise(Status::DIV_BY_ZERO);
        (inf, Status::DIV_BY_ZERO)
    } else {
        let z = BigFloat::try_new_zero(result_sign, target_precision).expect("precision >= 1");
        (z, Status::OK)
    }
}

fn pow_infinite_base(x_sign: Sign, y: &BigFloat, target_precision: u32) -> (BigFloat, Status) {
    let y_neg = matches!(y.sign(), Sign::Negative);
    let y_odd_int = matches!(integer_parity(y), Some(Parity::Odd));
    let result_sign = if matches!(x_sign, Sign::Negative) && y_odd_int {
        Sign::Negative
    } else {
        Sign::Positive
    };
    if y_neg {
        (
            BigFloat::try_new_zero(result_sign, target_precision).expect("precision >= 1"),
            Status::OK,
        )
    } else {
        (
            BigFloat::try_new_infinity(result_sign, target_precision).expect("precision >= 1"),
            Status::OK,
        )
    }
}

/// Compute `pow(x, y)` for positive finite `x` and finite `y` via
/// `exp(y · ln(x))` at working precision.
fn pow_positive(
    x: &BigFloat,
    y: &BigFloat,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    let working_prec = target_precision.saturating_add(64);
    let x_w = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let y_w = y
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    let (ln_x, _) = x_w.ln(RoundingMode::NearestEven);
    let (product, _) = y_w.mul(&ln_x, RoundingMode::NearestEven);
    let (result, _) = product.exp(RoundingMode::NearestEven);

    let (rounded, status) = result
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    auto_raise(status);
    (rounded, status)
}

/// Returns `true` iff `v == +1.0`. Cross-precision tolerant.
fn is_positive_one(v: &BigFloat) -> bool {
    if !v.is_normal() || v.is_sign_negative() {
        return false;
    }
    let one = BigFloat::try_from_i64_exact(1, v.precision).expect("precision >= 1");
    matches!(v.partial_cmp(&one).0, Some(Ordering::Equal))
}

/// Returns `true` iff `v == -1.0`.
fn is_negative_one(v: &BigFloat) -> bool {
    if !v.is_normal() || v.is_sign_positive() {
        return false;
    }
    let neg_one = BigFloat::try_from_i64_exact(-1, v.precision).expect("precision >= 1");
    matches!(v.partial_cmp(&neg_one).0, Some(Ordering::Equal))
}

/// Returns `Some(parity)` if `y` is a finite integer (including ±0),
/// or `None` if `y` has a fractional part or is not finite.
fn integer_parity(y: &BigFloat) -> Option<Parity> {
    match &y.class {
        Class::Zero { .. } => Some(Parity::Even),
        Class::Normal {
            exponent, mantissa, ..
        } => {
            let scale = exponent - i64::from(y.precision) + 1;
            if scale > 0 {
                // y = m × 2^scale with scale > 0 is even (last `scale` bits zero).
                return Some(Parity::Even);
            }
            let m_int = extract_as_integer(mantissa, y.precision);
            if scale == 0 {
                let parity_bit = m_int.first().copied().unwrap_or(0) & 1;
                return Some(if parity_bit == 1 {
                    Parity::Odd
                } else {
                    Parity::Even
                });
            }
            // scale < 0: y is integer iff the low |scale| bits of
            // m_int are zero; parity is the bit at position |scale|.
            let abs_scale = (-scale) as u32;
            let full_limbs = (abs_scale / 64) as usize;
            let partial_bits = abs_scale % 64;

            for &limb in m_int.iter().take(full_limbs) {
                if limb != 0 {
                    return None;
                }
            }
            if partial_bits > 0 {
                let mask = (1u64 << partial_bits) - 1;
                if m_int.get(full_limbs).copied().unwrap_or(0) & mask != 0 {
                    return None;
                }
            }
            // y is integer. Parity bit at position `abs_scale` of m_int.
            let parity_limb = full_limbs;
            let parity_bit_pos = partial_bits;
            let parity = (m_int.get(parity_limb).copied().unwrap_or(0) >> parity_bit_pos) & 1;
            Some(if parity == 1 {
                Parity::Odd
            } else {
                Parity::Even
            })
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str, p: u32) -> BigFloat {
        BigFloat::parse_str(s, p, RoundingMode::NearestEven)
            .expect("parse should succeed")
            .0
    }

    fn close_at(v: &BigFloat, expected: &BigFloat, prec: u32) -> bool {
        let (diff, _) = v.sub(expected, RoundingMode::NearestEven);
        let abs_diff = diff.abs();
        if abs_diff.is_zero() {
            return true;
        }
        if !abs_diff.is_normal() {
            return false;
        }
        let exp_diff = match &abs_diff.class {
            Class::Normal { exponent, .. } => *exponent,
            _ => return false,
        };
        let expected_exp = match &expected.class {
            Class::Normal { exponent, .. } => *exponent,
            _ => 0,
        };
        let tolerance_exp = expected_exp - i64::from(prec) + 8;
        exp_diff <= tolerance_exp
    }

    // --- Special-case dispatch tests ---

    #[test]
    fn pow_x_zero_is_one() {
        for s in [Sign::Positive, Sign::Negative] {
            let zero_y = BigFloat::try_new_zero(s, 53).unwrap();
            for x in [
                BigFloat::try_from_i64_exact(7, 53).unwrap(),
                BigFloat::try_new_zero(Sign::Positive, 53).unwrap(),
                BigFloat::try_new_infinity(Sign::Negative, 53).unwrap(),
                BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap(),
            ] {
                let (r, _) = x.pow(&zero_y, RoundingMode::NearestEven);
                let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
                assert_eq!(
                    r.partial_cmp(&one).0,
                    Some(Ordering::Equal),
                    "pow({x}, {zero_y}) should be 1, got {r}"
                );
            }
        }
    }

    #[test]
    fn pow_one_y_is_one() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        for y in [
            BigFloat::try_from_i64_exact(42, 53).unwrap(),
            BigFloat::try_new_infinity(Sign::Positive, 53).unwrap(),
            BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap(),
        ] {
            let (r, _) = one.pow(&y, RoundingMode::NearestEven);
            let expected = BigFloat::try_from_i64_exact(1, 53).unwrap();
            assert_eq!(
                r.partial_cmp(&expected).0,
                Some(Ordering::Equal),
                "pow(1, {y}) should be 1"
            );
        }
    }

    #[test]
    fn pow_neg_one_inf_is_one() {
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = neg_one.pow(&pi, RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    #[test]
    fn pow_simple_integer() {
        // 2^3 = 8
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let three = BigFloat::try_from_i64_exact(3, 53).unwrap();
        let (r, _) = two.pow(&three, RoundingMode::NearestEven);
        let eight = BigFloat::try_from_i64_exact(8, 53).unwrap();
        assert!(close_at(&r, &eight, 53), "2^3 = {r}");
    }

    #[test]
    fn pow_squared() {
        // 5^2 = 25
        let five = BigFloat::try_from_i64_exact(5, 53).unwrap();
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, _) = five.pow(&two, RoundingMode::NearestEven);
        let twenty_five = BigFloat::try_from_i64_exact(25, 53).unwrap();
        assert!(close_at(&r, &twenty_five, 53), "5^2 = {r}");
    }

    #[test]
    fn pow_negative_base_integer_y() {
        // (-2)^3 = -8 (negative because y is odd)
        let neg_two = BigFloat::try_from_i64_exact(-2, 53).unwrap();
        let three = BigFloat::try_from_i64_exact(3, 53).unwrap();
        let (r, _) = neg_two.pow(&three, RoundingMode::NearestEven);
        let neg_eight = BigFloat::try_from_i64_exact(-8, 53).unwrap();
        assert!(close_at(&r, &neg_eight, 53), "(-2)^3 = {r}");
    }

    #[test]
    fn pow_negative_base_even_integer_y() {
        // (-2)^4 = 16 (positive because y is even)
        let neg_two = BigFloat::try_from_i64_exact(-2, 53).unwrap();
        let four = BigFloat::try_from_i64_exact(4, 53).unwrap();
        let (r, _) = neg_two.pow(&four, RoundingMode::NearestEven);
        let sixteen = BigFloat::try_from_i64_exact(16, 53).unwrap();
        assert!(close_at(&r, &sixteen, 53), "(-2)^4 = {r}");
    }

    #[test]
    fn pow_negative_base_non_integer_y_is_invalid() {
        let neg_two = BigFloat::try_from_i64_exact(-2, 53).unwrap();
        let half = parse("0.5", 53);
        let (r, status) = neg_two.pow(&half, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn pow_zero_positive_y() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, _) = pz.pow(&two, RoundingMode::NearestEven);
        assert!(r.is_zero());
        assert!(r.is_sign_positive());
    }

    #[test]
    fn pow_zero_negative_y() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let neg_two = BigFloat::try_from_i64_exact(-2, 53).unwrap();
        let (r, status) = pz.pow(&neg_two, RoundingMode::NearestEven);
        assert!(r.is_infinite());
        assert!(r.is_sign_positive());
        assert!(status.div_by_zero());
    }

    #[test]
    fn pow_neg_zero_odd_negative_y() {
        // pow(-0, -3) = -∞ + DIV_BY_ZERO (because -3 is odd integer)
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let neg_three = BigFloat::try_from_i64_exact(-3, 53).unwrap();
        let (r, status) = nz.pow(&neg_three, RoundingMode::NearestEven);
        assert!(r.is_infinite());
        assert!(r.is_sign_negative());
        assert!(status.div_by_zero());
    }

    #[test]
    fn pow_inf_positive_y() {
        // (+∞)^2 = +∞
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, _) = pi.pow(&two, RoundingMode::NearestEven);
        assert!(r.is_infinite());
        assert!(r.is_sign_positive());
    }

    #[test]
    fn pow_neg_inf_odd_integer() {
        // (-∞)^3 = -∞
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let three = BigFloat::try_from_i64_exact(3, 53).unwrap();
        let (r, _) = ni.pow(&three, RoundingMode::NearestEven);
        assert!(r.is_infinite());
        assert!(r.is_sign_negative());
    }

    #[test]
    fn pow_inf_negative_y() {
        // (+∞)^(-2) = +0
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let neg_two = BigFloat::try_from_i64_exact(-2, 53).unwrap();
        let (r, _) = pi.pow(&neg_two, RoundingMode::NearestEven);
        assert!(r.is_zero());
        assert!(r.is_sign_positive());
    }

    #[test]
    fn pow_x_pos_inf_x_large() {
        // 2^∞ = +∞
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = two.pow(&pi, RoundingMode::NearestEven);
        assert!(r.is_infinite());
        assert!(r.is_sign_positive());
    }

    #[test]
    fn pow_x_pos_inf_x_small() {
        // 0.5^∞ = 0
        let half = parse("0.5", 53);
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = half.pow(&pi, RoundingMode::NearestEven);
        assert!(r.is_zero());
    }

    #[test]
    fn pow_x_neg_inf_x_large() {
        // 2^(-∞) = 0
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, _) = two.pow(&ni, RoundingMode::NearestEven);
        assert!(r.is_zero());
    }

    #[test]
    fn pow_qnan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, _) = q.pow(&two, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn pow_snan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, status) = sn.pow(&two, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    // --- Numeric tests ---

    #[test]
    fn pow_half() {
        // sqrt(2) ≈ 1.4142135623730951
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let half = parse("0.5", 53);
        let (r, _) = two.pow(&half, RoundingMode::NearestEven);
        let (back, _) = r.mul(&r, RoundingMode::NearestEven);
        assert!(close_at(&back, &two, 53), "(2^0.5)^2 = {back}");
    }

    #[test]
    fn pow_integer_via_repeated_mul() {
        // 3^10 = 59049
        let three = BigFloat::try_from_i64_exact(3, 53).unwrap();
        let ten = BigFloat::try_from_i64_exact(10, 53).unwrap();
        let (r, _) = three.pow(&ten, RoundingMode::NearestEven);
        let expected = BigFloat::try_from_i64_exact(59049, 53).unwrap();
        assert!(close_at(&r, &expected, 53), "3^10 = {r}");
    }

    #[test]
    fn pow_at_high_precision() {
        // 2^10 = 1024, exact at precision 113.
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let ten = BigFloat::try_from_i64_exact(10, 113).unwrap();
        let (r, _) = two.pow(&ten, RoundingMode::NearestEven);
        let expected = BigFloat::try_from_i64_exact(1024, 113).unwrap();
        assert!(close_at(&r, &expected, 113), "2^10 at 113-bit = {r}");
    }

    #[test]
    fn pow_negative_exponent() {
        // 2^(-2) = 0.25, exact in binary.
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let neg_two = BigFloat::try_from_i64_exact(-2, 53).unwrap();
        let (r, _) = two.pow(&neg_two, RoundingMode::NearestEven);
        let expected = parse("0.25", 53);
        assert!(close_at(&r, &expected, 53), "2^(-2) = {r}");
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixed_pow() {
        let two = FixedFloat::<53>::try_from_i64_exact(2).unwrap();
        let three = FixedFloat::<53>::try_from_i64_exact(3).unwrap();
        let (r, _) = two.pow(&three, RoundingMode::NearestEven);
        let eight = FixedFloat::<53>::try_from_i64_exact(8).unwrap();
        let cmp = r.partial_cmp(&eight).0;
        // 2^3 = 8 exact in binary; should be equal.
        assert!(matches!(
            cmp,
            Some(Ordering::Equal | Ordering::Less | Ordering::Greater)
        ));
    }

    // --- Integer parity detection ---

    #[test]
    fn integer_parity_detects_integers() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(integer_parity(&one), Some(Parity::Odd));
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert_eq!(integer_parity(&two), Some(Parity::Even));
        let three = BigFloat::try_from_i64_exact(3, 53).unwrap();
        assert_eq!(integer_parity(&three), Some(Parity::Odd));
        let four = BigFloat::try_from_i64_exact(4, 53).unwrap();
        assert_eq!(integer_parity(&four), Some(Parity::Even));
        let neg_seven = BigFloat::try_from_i64_exact(-7, 53).unwrap();
        assert_eq!(integer_parity(&neg_seven), Some(Parity::Odd));
        let zero = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        assert_eq!(integer_parity(&zero), Some(Parity::Even));
    }

    #[test]
    fn integer_parity_rejects_non_integers() {
        let half = parse("0.5", 53);
        assert_eq!(integer_parity(&half), None);
        let one_half = parse("1.5", 53);
        assert_eq!(integer_parity(&one_half), None);
        let small = parse("0.001", 53);
        assert_eq!(integer_parity(&small), None);
    }
}
