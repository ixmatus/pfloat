//! Argument reduction for the forward trig kernels.
//!
//! Given a finite `x`, returns a quadrant `q ∈ {0, 1, 2, 3}` and a
//! reduced argument `r ∈ [−π/4, π/4]` such that
//!
//! ```text
//! x = (q + (r / (π/2))) · (π/2)
//! ```
//!
//! up to the working precision passed in. The implementation
//! multiplies `x` by the hardcoded 4096-bit `2/π` table from
//! [`super::TWO_OVER_PI_LIMBS_4096`], rounds the product to the
//! nearest integer to get `q` (modulo 4 for the quadrant index),
//! and scales the fractional remainder by `π/2`.
//!
//! Quadrant convention:
//!
//! | q | sin(x)   | cos(x)   |
//! |---|----------|----------|
//! | 0 | +sin(r)  | +cos(r)  |
//! | 1 | +cos(r)  | −sin(r)  |
//! | 2 | −sin(r)  | −cos(r)  |
//! | 3 | −cos(r)  | +sin(r)  |
//!
//! Range cap: the 4096-bit `2/π` table supports inputs with binary
//! exponent up to roughly `4096 − target_precision − slack`. For
//! `|x|` past that band the reduction would lose bits, so
//! [`reduce`] returns `None`; the caller is expected to translate
//! that into `INVALID` with a quiet NaN.

use core::cmp::Ordering;

use crate::big::BigFloat;
use crate::class::Class;
use crate::ops::limbs::extract_as_integer;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

use super::{pi_over_2_at, two_over_pi_at};

/// The fractional/integer split produced by the reduction.
pub(super) struct Reduction {
    /// Quadrant index, `0..=3`.
    pub quadrant: u8,
    /// Reduced argument, `|r| ≤ π/4`, at the requested working
    /// precision.
    pub r: BigFloat,
}

/// Reduces a finite normal `x` to `(quadrant, r)`. Returns `None`
/// when `|x|` exceeds the table's reduction budget.
pub(super) fn reduce(x: &BigFloat, working_prec: u32) -> Option<Reduction> {
    // Determine the input's binary exponent. For Zero/special, the
    // caller has already dispatched.
    let e_x = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        Class::Zero { .. } => {
            // x = ±0 → quadrant 0, r = 0.
            return Some(Reduction {
                quadrant: 0,
                r: BigFloat::try_new_zero(x.sign(), working_prec).expect("precision >= 1"),
            });
        }
        _ => return None,
    };

    // Range check: table is 4096 bits. Multiplying x · (2/π)
    // produces a value whose bits relevant to the reduction span
    // positions roughly [e_x − 1, e_x − 1 − working_prec − slack].
    // For the lowest of these to fall inside the table we need
    // `e_x + working_prec + 64 < 4096`.
    let slack = i64::from(working_prec) + 64;
    // Saturating: an exponent near i64::MAX (reachable by repeated
    // squaring, which saturates the exponent under the no-emax design)
    // must route to the out-of-range `None` path, not overflow i64 and
    // wrap below 4096. Review 2026-05-29.
    if e_x.saturating_add(slack) >= 4096 {
        return None;
    }

    // For |x| ≤ π/4 (exponent ≤ −1), no reduction needed.
    // Strictly: π/4 ≈ 0.785, so |x| < 1 has exponent ≤ −1. Pre-cut
    // to avoid the round-trip through the multiplier.
    if e_x <= -1 {
        // We still need to verify |x| ≤ π/4; for x with exponent
        // exactly −1, x ∈ [−1, −0.5] ∪ [0.5, 1] and π/4 ≈ 0.785
        // sits inside the latter interval. For safety we fall
        // through to the general path when exponent == 0 may push
        // us across π/4. Skip the shortcut for `exponent == −1` and
        // do the multiplier path uniformly.
        if e_x <= -2 {
            // |x| < 0.5 < π/4, definitely in range.
            let r = x
                .round_to_precision(working_prec, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            return Some(Reduction { quadrant: 0, r });
        }
    }

    // The product table multiplication picks max(x.precision, 4096)
    // bits of precision; passing it through explicit
    // `mul_round` with a chosen working width is cleaner. Floor at
    // 2048 so even tight-input reductions get enough headroom, and
    // cap at 4096 to match the hardcoded table width.
    // The product x·(2/π) has magnitude ~2^e_x; to retain
    // `working_prec + 64` accurate bits BELOW the binary point (the
    // reduced argument is that fractional remainder) the product must
    // carry `e_x + working_prec + 64` significant bits. The range check
    // above guarantees this is < 4096, so the 4096-bit table covers it.
    // The earlier formula omitted e_x and silently truncated the
    // reduction to the 2048-bit floor for large |x|, collapsing
    // sin/cos/tan to 0/±1 with no flag (review 2026-05-29, root cause 1).
    let needed = e_x
        .max(0)
        .saturating_add(i64::from(working_prec))
        .saturating_add(64);
    let mul_prec = u32::try_from(needed)
        .unwrap_or(4096)
        .max(x.precision)
        .max(2048)
        .clamp(2048, 4096);
    let two_over_pi = two_over_pi_at(mul_prec);
    let (y, _) = x
        .mul_round(&two_over_pi, mul_prec, RoundingMode::NearestEven)
        .expect("precision >= 1");

    // Round y to nearest integer with banker's tie-breaking. The
    // helper handles |y| < 0.5 (rounds to 0), |y| in [0.5, 1)
    // (round-to-even rounds 0.5 to 0, anything above to ±1), and
    // |y| ≥ 1 (round-to-precision with precision = exponent + 1).
    let q_int = round_to_nearest_int(&y);

    // r_unscaled = y − q. Exact in real math; rounded to working
    // precision so we don't carry the full mul-width afterward.
    let (r_unscaled, _) = y.sub(&q_int, RoundingMode::NearestEven);

    // quadrant = q mod 4 (signed mod, taken into [0, 4)).
    let quadrant = mod_4(&q_int);

    let pi_over_2 = pi_over_2_at(working_prec.saturating_add(64));
    let (r_scaled, _) = r_unscaled.mul(&pi_over_2, RoundingMode::NearestEven);
    let (r, _) = r_scaled
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1");

    Some(Reduction { quadrant, r })
}

/// Round `y` to the nearest integer (banker's tie-breaking),
/// returning the integer as a `BigFloat`.
fn round_to_nearest_int(y: &BigFloat) -> BigFloat {
    match &y.class {
        Class::Zero { .. } => y.clone(),
        Class::Normal {
            exponent,
            sign,
            mantissa,
        } => match (*exponent).cmp(&-1) {
            Ordering::Less => {
                // |y| < 0.5 → 0 (keep y's sign for downstream r computation).
                BigFloat::try_new_zero(*sign, y.precision).expect("precision >= 1")
            }
            Ordering::Equal => {
                // |y| ∈ [0.5, 1). Banker's: y = ±0.5 → 0; else → ±1.
                let m_int = extract_as_integer(mantissa, y.precision);
                if only_top_bit_set(&m_int, y.precision) {
                    BigFloat::try_new_zero(*sign, y.precision).expect("precision >= 1")
                } else {
                    let value = if matches!(sign, Sign::Negative) {
                        -1
                    } else {
                        1
                    };
                    BigFloat::try_from_i64_exact(value, y.precision).expect("precision >= 1")
                }
            }
            Ordering::Greater => {
                // exponent ≥ 0. Round mantissa to (exponent + 1) bits.
                let int_prec = u32::try_from(*exponent + 1).unwrap_or(u32::MAX);
                y.round_to_precision(int_prec, RoundingMode::NearestEven)
                    .expect("precision >= 1")
                    .0
            }
        },
        _ => panic!("round_to_nearest_int requires finite y"),
    }
}

/// Check whether a mantissa-as-integer (extracted via
/// [`extract_as_integer`]) has only its top bit set, i.e., the
/// mantissa-integer equals `2^(precision − 1)`.
fn only_top_bit_set(m_int: &[u64], precision: u32) -> bool {
    let top_pos = (precision - 1) as usize;
    let top_limb = top_pos / 64;
    let top_bit = top_pos % 64;
    let expected = 1u64 << top_bit;
    if m_int.get(top_limb).copied().unwrap_or(0) != expected {
        return false;
    }
    for (i, &limb) in m_int.iter().enumerate() {
        if i != top_limb && limb != 0 {
            return false;
        }
    }
    true
}

/// `(value as signed integer) mod 4`, mapped into `[0, 4)`.
///
/// `q` must represent a finite integer (`round_to_nearest_int`
/// guarantees this). The mantissa-as-integer is extracted via
/// [`extract_as_integer`] so positions are tight; the value is
/// `m_int · 2^(exponent − precision + 1)`.
fn mod_4(q: &BigFloat) -> u8 {
    match &q.class {
        Class::Zero { .. } => 0,
        Class::Normal {
            sign,
            exponent,
            mantissa,
        } => {
            if *exponent < 0 {
                return 0;
            }
            let scale = *exponent - i64::from(q.precision) + 1;
            let m_int = extract_as_integer(mantissa, q.precision);
            // Bit position in m_int that corresponds to value-bit 0:
            // bit_0(value) = bit_(−scale)(m_int).
            let abs_low = if scale >= 2 {
                0
            } else if scale == 1 {
                (m_int.first().copied().unwrap_or(0) & 1) << 1
            } else if scale == 0 {
                m_int.first().copied().unwrap_or(0) & 3
            } else {
                let pos = (-scale) as usize;
                let limb0 = pos / 64;
                let bit0 = pos % 64;
                let lo_limb = m_int.get(limb0).copied().unwrap_or(0);
                let bit_0_val = (lo_limb >> bit0) & 1;
                let bit_1_val = if bit0 == 63 {
                    m_int.get(limb0 + 1).copied().unwrap_or(0) & 1
                } else {
                    (lo_limb >> (bit0 + 1)) & 1
                };
                (bit_1_val << 1) | bit_0_val
            };
            let abs_mod4 = (abs_low & 3) as u8;
            if matches!(sign, Sign::Negative) {
                (4 - abs_mod4) & 3
            } else {
                abs_mod4
            }
        }
        _ => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close_at(v: &BigFloat, expected: &BigFloat, bits: u32) -> bool {
        let (diff, _) = v.sub(expected, RoundingMode::NearestEven);
        let abs_diff = diff.abs();
        if abs_diff.is_zero() {
            return true;
        }
        let p = v.precision().max(expected.precision());
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let abs_b = expected.abs();
        let mut bound = if abs_b.is_zero() { one } else { abs_b };
        for _ in 0..bits {
            bound = bound.div(&two, RoundingMode::NearestEven).0;
        }
        matches!(
            abs_diff.partial_cmp(&bound).0,
            Some(Ordering::Less | Ordering::Equal)
        )
    }

    #[test]
    fn pi_constants_self_consistent() {
        // π · (2/π) should equal 2 at the table's working precision.
        let pi = super::super::pi_at(1024);
        let two_over_pi = super::super::two_over_pi_at(4096);
        let (prod, _) = pi.mul(&two_over_pi, RoundingMode::NearestEven);
        let two = BigFloat::try_from_i64_exact(2, 4096).unwrap();
        assert!(
            close_at(&prod, &two, 1000),
            "π · (2/π) = {prod}, expected 2 within 1000 bits",
        );
    }

    #[test]
    fn reduce_zero_returns_quadrant_zero() {
        let z = BigFloat::try_new_zero(Sign::Positive, 113).unwrap();
        let r = reduce(&z, 113).unwrap();
        assert_eq!(r.quadrant, 0);
        assert!(r.r.is_zero());
    }

    #[test]
    fn reduce_small_x_passthrough() {
        // x = 0.1, |x| < π/4, so r = x at working precision.
        let x = BigFloat::parse_str("0.1", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let r = reduce(&x, 113).unwrap();
        assert_eq!(r.quadrant, 0);
        assert!(close_at(&r.r, &x, 113 - 12));
    }

    #[test]
    fn reduce_pi_over_2_gives_quadrant_one() {
        let pi_2 = super::super::pi_over_2_at(113);
        let r = reduce(&pi_2, 113).unwrap();
        assert_eq!(r.quadrant, 1);
        // r should be near zero (since x = π/2 maps to q=1, r=0).
        let zero = BigFloat::try_new_zero(Sign::Positive, 113).unwrap();
        // r could be tiny but nonzero from rounding; allow loose tolerance.
        assert!(close_at(&r.r, &zero, 100));
    }

    #[test]
    fn reduce_pi_gives_quadrant_two() {
        let pi = super::super::pi_at(113);
        let r = reduce(&pi, 113).unwrap();
        assert_eq!(r.quadrant, 2);
    }

    #[test]
    fn reduce_neg_pi_over_2_gives_quadrant_three() {
        let pi_2 = super::super::pi_over_2_at(113);
        let neg = pi_2.negated();
        let r = reduce(&neg, 113).unwrap();
        assert_eq!(r.quadrant, 3);
    }

    #[test]
    fn reduce_huge_x_returns_none() {
        // Build x with exponent 5000, exceeding the table budget.
        // Easiest: 2^5000.
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let mut big = one;
        for _ in 0..5000 {
            big = big.mul(&two, RoundingMode::NearestEven).0;
        }
        assert!(reduce(&big, 113).is_none());
    }
}
