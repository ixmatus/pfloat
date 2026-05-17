//! `pow(x, y) = x^y`: general power function.
//!
//! For positive `x` and finite `y` the result is correctly rounded
//! under every IEEE rounding mode (slice 7c, ADR-0022). Two paths
//! feed a shared Ziv driver [`pow_ziv`]:
//!
//! - **Integer exponent**: when `y` is an exact integer in range,
//!   `x^|n|` is formed by square-and-multiply ([`pow_int`]) and
//!   reciprocated for `n < 0`. Exact cases (`2^10`) round bit-exactly
//!   at the first guard, matching MPFR's integer fast path.
//! - **General exponent**: `exp(y · ln(x))` evaluated at working
//!   precision (the slice-3a `exp` and slice-3b `ln` carry the work).
//!
//! [`pow_ziv`] realizes DESIGN.md §"Ziv's strategy" by the
//! recompute-and-compare test: it evaluates at a guard, rounds to
//! the target under the caller's mode, then re-evaluates at a larger
//! guard; two independent higher-precision evaluations rounding to
//! the same target value means the result is not boundary-ambiguous.
//! On disagreement the guard doubles and the loop retries, capped at
//! [`ZIV_MAX_ITERS`] (the honest pathological-input caveat MPFR also
//! carries). pfloat is the first transcendental off the
//! NearestEven-only differential tier as a result.
//!
//! This module also handles the rich IEEE 754-2019 §9.2.1
//! special-case table, which short-circuits before either path.
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
use crate::ops::limbs::{extract_as_integer, top_set_bit};
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
    ///
    /// For positive base and finite exponent the result is correctly
    /// rounded under `mode`, subject to the Ziv iteration cap
    /// [`ZIV_MAX_ITERS`]: on the measure-zero exact-tie inputs that
    /// exhaust the cap the result may be 1 ULP off in directed modes,
    /// the same caveat MPFR documents (DESIGN.md §"Ziv's strategy").
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

/// First Ziv guard: the initial evaluation uses
/// `target + ZIV_BASE_GUARD` extra bits.
const ZIV_BASE_GUARD: u32 = 64;

/// Maximum extra guard bits above the target precision. The doubling
/// schedule (64, 128, 256, 512, 1024) reaches this at the last
/// iteration.
const ZIV_GUARD_CAP: u32 = 1024;

/// Maximum recompute-and-compare iterations. On the measure-zero
/// exact-tie inputs that exhaust this many iterations the result may
/// be 1 ULP off in directed modes — the honest caveat MPFR also
/// documents (DESIGN.md §"Ziv's strategy", lines 287-299).
pub(super) const ZIV_MAX_ITERS: u32 = 5;

/// Bit-exact (binary-radix canonical) equality of two values already
/// rounded to the same target precision.
fn same_rounded(a: &BigFloat, b: &BigFloat) -> bool {
    matches!(a.partial_cmp(b).0, Some(Ordering::Equal))
}

/// Correctly round `eval`'s value to `target` precision under `mode`
/// by the Ziv recompute-and-compare test (DESIGN.md §"Ziv's
/// strategy").
///
/// `eval(working)` returns the function value computed at the given
/// working precision with `NearestEven` internal rounding; the
/// directed final round to `target` is applied here. Two consecutive
/// guards whose `target`-rounded values agree settle the result: a
/// boundary case would round differently as the guard grows. On
/// disagreement the guard doubles (capped at [`ZIV_GUARD_CAP`]) and
/// the loop retries, bounded by [`ZIV_MAX_ITERS`]; if the cap is
/// reached the last iteration's value is returned best-effort.
fn pow_ziv(eval: impl Fn(u32) -> BigFloat, target: u32, mode: RoundingMode) -> (BigFloat, Status) {
    let mut guard = ZIV_BASE_GUARD;
    let mut prev: Option<BigFloat> = None;
    let mut last: Option<(BigFloat, Status)> = None;
    for _ in 0..ZIV_MAX_ITERS {
        let working = target.saturating_add(guard);
        let hi = eval(working);
        let (cand, status) = hi
            .round_to_precision(target, mode)
            .expect("target precision >= 1");
        if let Some(p) = &prev {
            if same_rounded(&cand, p) {
                auto_raise(status);
                return (cand, status);
            }
        }
        prev = Some(cand.clone());
        last = Some((cand, status));
        guard = guard.saturating_mul(2).min(ZIV_GUARD_CAP);
    }
    let (cand, status) = last.expect("ZIV_MAX_ITERS >= 1");
    auto_raise(status);
    (cand, status)
}

/// Conservative upper bound on the integer fast path's result binary
/// exponent. Beyond this the `exp·ln` path runs instead: it carries
/// the correct OVERFLOW/UNDERFLOW status and avoids a pathological
/// square-and-multiply on an astronomically out-of-range result.
/// Generous: any realistic correctly-rounded use stays far below it.
const POW_INT_RESULT_EXPONENT_CAP: i128 = 1 << 24;

/// Compute `pow(x, y)` for positive finite `x` and finite `y`,
/// correctly rounded under `mode` via the shared [`pow_ziv`] driver.
///
/// An exact integer exponent in range takes the square-and-multiply
/// path ([`pow_int`]); everything else evaluates `exp(y · ln(x))`.
/// Both feed [`pow_ziv`] for the directed final round.
fn pow_positive(
    x: &BigFloat,
    y: &BigFloat,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    if let Some(n) = integer_exponent(y) {
        if let Some(result) = pow_int_path(x, n, target_precision, mode) {
            return result;
        }
        // Out of the feasible result-exponent range: fall through to
        // `exp·ln`, which produces the correct ±∞/±0 + OVERFLOW /
        // UNDERFLOW status for the extreme case.
    }

    pow_ziv(
        |w| {
            let x_w = x
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("w >= 1")
                .0;
            let y_w = y
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("w >= 1")
                .0;
            let (ln_x, _) = x_w.ln(RoundingMode::NearestEven);
            let (product, _) = y_w.mul(&ln_x, RoundingMode::NearestEven);
            let (result, _) = product.exp(RoundingMode::NearestEven);
            result
        },
        target_precision,
        mode,
    )
}

/// Integer-exponent fast path: `x^n` by square-and-multiply through
/// [`pow_ziv`]. Returns `None` (deferring to `exp·ln`) when the
/// predicted result exponent is past [`POW_INT_RESULT_EXPONENT_CAP`]
/// or the computed value overflowed/underflowed (so the `exp·ln`
/// path can raise the correct status).
fn pow_int_path(
    x: &BigFloat,
    n: i64,
    target_precision: u32,
    mode: RoundingMode,
) -> Option<(BigFloat, Status)> {
    let e_x = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => return None,
    };
    // |x|^n has binary exponent ≈ n·log₂|x|; bound it above by
    // |n|·(|eₓ|+1) before doing any work.
    let est = i128::from(n).saturating_mul(i128::from(e_x.unsigned_abs()) + 1);
    if est.unsigned_abs() > POW_INT_RESULT_EXPONENT_CAP.unsigned_abs() {
        return None;
    }
    let n_abs = n.unsigned_abs();
    let (bf, status) = pow_ziv(
        |w| {
            let p = pow_int(x, n_abs, w);
            if n < 0 {
                let one = BigFloat::try_from_i64_exact(1, w).expect("w >= 1");
                one.div(&p, RoundingMode::NearestEven).0
            } else {
                p
            }
        },
        target_precision,
        mode,
    );
    if bf.is_infinite() || bf.is_zero() {
        // Over/underflowed; let `exp·ln` produce the IEEE status.
        None
    } else {
        Some((bf, status))
    }
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

/// Returns `Some(n)` if `y` is an exact finite integer whose value
/// fits in `i64`, else `None`.
///
/// Sibling of [`integer_parity`]: same `scale` decomposition
/// (`exponent − precision + 1`), but it reconstructs the signed
/// value rather than just the parity bit. `None` covers a fractional
/// part, a non-finite `y`, and a magnitude past `i64::MAX`. The
/// result-exponent (`n·eₓ`) feasibility guard lives at the dispatch
/// site (it needs the base); this extractor only decides integrality
/// and `i64` range, so an out-of-range integer falls back to the
/// `exp·ln` path (which handles overflow/underflow via `exp`).
fn integer_exponent(y: &BigFloat) -> Option<i64> {
    match &y.class {
        Class::Zero { .. } => Some(0),
        Class::Normal {
            exponent, mantissa, ..
        } => {
            let scale = exponent - i64::from(y.precision) + 1;
            let m_int = extract_as_integer(mantissa, y.precision);
            // `m_int` is normalized (top bit set) so `top_set_bit`
            // is `Some`; `lsb` is the lowest set bit.
            let msb = top_set_bit(&m_int)? as i64;
            let lsb = lowest_set_bit(&m_int)? as i64;
            // y = m_int · 2^scale. Integer iff no set bit lands below
            // the binary point, i.e. lsb + scale ≥ 0.
            if lsb + scale < 0 {
                return None;
            }
            // Exponent of the integer's most significant bit. Keep it
            // below 63 so the reconstructed magnitude fits i64.
            let top = msb + scale;
            if top > 62 {
                return None;
            }
            let mut mag: i128 = 0;
            for i in lsb..=msb {
                let bit = (m_int[(i / 64) as usize] >> (i % 64)) & 1;
                if bit == 1 {
                    mag |= 1i128 << (i + scale);
                }
            }
            let signed = i64::try_from(mag).ok()?;
            Some(if y.is_sign_negative() {
                -signed
            } else {
                signed
            })
        }
        _ => None,
    }
}

/// Lowest set bit index of a little-endian limb buffer, or `None`
/// if every limb is zero.
fn lowest_set_bit(buffer: &[u64]) -> Option<usize> {
    for (i, &limb) in buffer.iter().enumerate() {
        if limb != 0 {
            return Some(i * 64 + limb.trailing_zeros() as usize);
        }
    }
    None
}

/// `x^n_abs` for a non-negative integer exponent by
/// square-and-multiply at working precision `w` (`NearestEven`
/// internal rounding). The caller reciprocates for negative `n` and
/// feeds the result through [`pow_ziv`] for the directed final
/// round, so exact powers settle bit-exactly at the first guard.
fn pow_int(x: &BigFloat, n_abs: u64, w: u32) -> BigFloat {
    let mut result = BigFloat::try_from_i64_exact(1, w).expect("w >= 1");
    let mut base = x
        .round_to_precision(w, RoundingMode::NearestEven)
        .expect("w >= 1")
        .0;
    let mut e = n_abs;
    while e > 0 {
        if e & 1 == 1 {
            result = result.mul(&base, RoundingMode::NearestEven).0;
        }
        e >>= 1;
        if e > 0 {
            base = base.mul(&base, RoundingMode::NearestEven).0;
        }
    }
    result
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

    // --- Ziv driver (slice 7c.1) ---

    #[test]
    fn pow_ziv_exact_value_is_bit_exact() {
        // An exact value (8) is identical at every working precision,
        // so the first two guards agree immediately and the result is
        // returned bit-exactly — the MPFR integer-fast-path parity
        // the integer branch relies on.
        for &mode in &[
            RoundingMode::NearestEven,
            RoundingMode::TowardPositive,
            RoundingMode::TowardNegative,
            RoundingMode::TowardZero,
            RoundingMode::NearestAway,
        ] {
            let (r, _) = pow_ziv(
                |w| BigFloat::try_from_i64_exact(8, w).expect("w >= 1"),
                53,
                mode,
            );
            let eight = BigFloat::try_from_i64_exact(8, 53).unwrap();
            assert_eq!(
                r.partial_cmp(&eight).0,
                Some(Ordering::Equal),
                "pow_ziv exact 8 under {mode:?} = {r}"
            );
            assert_eq!(r.precision, 53);
        }
    }

    #[test]
    fn pow_ziv_converges_to_correctly_rounded() {
        // A transcendental-shaped constant: as the guard grows the
        // target-rounding stabilizes, and pow_ziv must land on the
        // value obtained by rounding the constant directly to the
        // target precision under the same mode (correct rounding).
        let digits = "1.41421356237309504880168872420969807856967187537694\
                       8073176679737990732478462107038850387534327641572735";
        for &mode in &[
            RoundingMode::NearestEven,
            RoundingMode::TowardPositive,
            RoundingMode::TowardNegative,
            RoundingMode::TowardZero,
            RoundingMode::NearestAway,
        ] {
            let target = 113u32;
            let (r, _) = pow_ziv(
                |w| {
                    BigFloat::parse_str(digits, w, RoundingMode::NearestEven)
                        .expect("parse")
                        .0
                },
                target,
                mode,
            );
            let direct = BigFloat::parse_str(digits, target, mode).expect("parse").0;
            assert_eq!(
                r.partial_cmp(&direct).0,
                Some(Ordering::Equal),
                "pow_ziv vs direct round under {mode:?}: ziv={r}, direct={direct}"
            );
        }
    }

    #[test]
    fn pow_ziv_is_idempotent_under_recompute() {
        // Stable input → stable output across an independent re-run:
        // the driver does not introduce nondeterminism.
        let eval = |w: u32| {
            BigFloat::parse_str(
                "2.7182818284590452353602874713527",
                w,
                RoundingMode::NearestEven,
            )
            .expect("parse")
            .0
        };
        let (a, _) = pow_ziv(eval, 80, RoundingMode::NearestEven);
        let (b, _) = pow_ziv(eval, 80, RoundingMode::NearestEven);
        assert_eq!(a.partial_cmp(&b).0, Some(Ordering::Equal));
    }

    // --- Integer exponent extractor + pow_int (slice 7c.2) ---

    #[test]
    fn integer_exponent_detects_across_precisions() {
        // The trailing-zero mantissa of a small integer grows with
        // precision; the extractor must still recover the value at
        // every precision the differential sweep exercises.
        for &p in &[53u32, 113, 256, 1024] {
            for n in [0i64, 1, -1, 7, -7, 42, 1024, -1023] {
                let v = BigFloat::try_from_i64_exact(n, p).unwrap();
                assert_eq!(integer_exponent(&v), Some(n), "n={n} p={p}");
            }
        }
    }

    #[test]
    fn integer_exponent_rejects_non_integers() {
        for s in ["0.5", "1.5", "0.001", "-3.25"] {
            assert_eq!(integer_exponent(&parse(s, 53)), None, "{s}");
        }
    }

    #[test]
    fn integer_exponent_rejects_out_of_i64_range() {
        // 2^62 fits i64; 2^63 does not (top bit exponent 63 > 62).
        let two_62 = parse("4611686018427387904", 80);
        assert_eq!(integer_exponent(&two_62), Some(1i64 << 62));
        let two_63 = parse("9223372036854775808", 80);
        assert_eq!(integer_exponent(&two_63), None);
    }

    #[test]
    fn pow_int_equals_repeated_mul() {
        // Exactly-representable integer powers: square-and-multiply
        // must equal the closed-form integer bit-for-bit.
        for (x, n, expect) in [(3i64, 10u64, 59049i64), (2, 20, 1_048_576), (7, 3, 343)] {
            let xb = BigFloat::try_from_i64_exact(x, 113).unwrap();
            let r = pow_int(&xb, n, 113 + 64);
            let (rr, _) = r
                .round_to_precision(113, RoundingMode::NearestEven)
                .unwrap();
            let want = BigFloat::try_from_i64_exact(expect, 113).unwrap();
            assert_eq!(
                rr.partial_cmp(&want).0,
                Some(Ordering::Equal),
                "{x}^{n} = {rr}, want {expect}"
            );
        }
    }

    #[test]
    fn pow_int_zero_exponent_is_one() {
        let x = BigFloat::try_from_i64_exact(5, 53).unwrap();
        let r = pow_int(&x, 0, 53);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal));
    }

    // --- Integer-path correct rounding under directed modes ---

    #[test]
    fn pow_directed_modes_integer_exponent() {
        // ADR-0014's witness: pow(63, 9) at p=53 differed from MPFR
        // by 1 ULP under the old single-shot path. 63^9 ≈ 2^53.8 so
        // it is not representable at p=53 and must round. Build the
        // exact value by an independent naive repeated multiply at
        // ample precision, then check the public pow equals that
        // value rounded to 53 bits under each mode (correct
        // rounding, no MPFR needed).
        let sixty_three = BigFloat::try_from_i64_exact(63, 128).unwrap();
        let mut exact = BigFloat::try_from_i64_exact(1, 128).unwrap();
        for _ in 0..9 {
            exact = exact.mul(&sixty_three, RoundingMode::NearestEven).0;
        }
        for &mode in &[
            RoundingMode::NearestEven,
            RoundingMode::NearestAway,
            RoundingMode::TowardZero,
            RoundingMode::TowardPositive,
            RoundingMode::TowardNegative,
        ] {
            let want = exact.round_to_precision(53, mode).unwrap().0;
            let base = BigFloat::try_from_i64_exact(63, 53).unwrap();
            let expn = BigFloat::try_from_i64_exact(9, 53).unwrap();
            let (got, _) = base.pow(&expn, mode);
            assert_eq!(
                got.partial_cmp(&want).0,
                Some(Ordering::Equal),
                "pow(63,9) under {mode:?}: got {got}, want {want}"
            );
        }
    }

    #[test]
    fn pow_negative_base_integer_uses_fast_path() {
        // (-63)^9 routes through pow_positive(|x|,…) then negates;
        // it must equal -(63^9 rounded to p), odd exponent.
        let s = BigFloat::try_from_i64_exact(63, 128).unwrap();
        let mut exact = BigFloat::try_from_i64_exact(1, 128).unwrap();
        for _ in 0..9 {
            exact = exact.mul(&s, RoundingMode::NearestEven).0;
        }
        let want = exact
            .round_to_precision(53, RoundingMode::NearestEven)
            .unwrap()
            .0
            .negated();
        let base = BigFloat::try_from_i64_exact(-63, 53).unwrap();
        let expn = BigFloat::try_from_i64_exact(9, 53).unwrap();
        let (got, _) = base.pow(&expn, RoundingMode::NearestEven);
        assert_eq!(got.partial_cmp(&want).0, Some(Ordering::Equal));
    }
}
