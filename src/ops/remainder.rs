//! IEEE 754-2019 §5.3.1 `remainder`.
//!
//! `remainder(x, y) = x - n·y`, where `n` is the integer quotient
//! `x / y` rounded to nearest with ties to even. The result is exact
//! (it never rounds) and satisfies `|remainder| <= |y|/2`. This is the
//! operation `core::ops::Rem` (`%`) resolves to for [`BigFloat`] and
//! `FixedFloat<PREC>`; it differs from C `fmod` (which truncates the
//! quotient), and pfloat documents that choice rather than aliasing
//! `%` to a non-IEEE operation.
//!
//! The quotient `n` can have an exponent-sized integer part (up to the
//! `i64` exponent range), so forming `n` and computing `x - n·y`
//! directly would be a denial-of-service vector. The kernel instead
//! reduces the mantissa integers modulo each other: with
//! `|x| = Mx·2^ax` and `|y| = My·2^ay` (mantissa integers, bottom
//! exponents), the truncated remainder is `(X mod Y)·2^s` where
//! `s = min(ax, ay)` and the larger operand is scaled to `s`. When `x`
//! dominates (`ax > ay`) the scaling factor `2^(ax-ay)` is reduced by
//! modular exponentiation (`O(log)` multiplications), never
//! materialized. When `y` dominates, the early exit `2|x| < |y| → x`
//! bounds the opposite shift, so that branch is materialized safely.
//! The round-to-nearest-even adjustment then compares `2R` to the
//! scaled `|y|`, consulting the quotient parity only on an exact tie.

#[cfg(feature = "alloc")]
extern crate alloc;

use alloc::vec;
use alloc::vec::Vec;

use crate::big::BigFloat;
use crate::class::Class;
use crate::mantissa::limbs_for;
use crate::ops::limbs::{
    cmp_limbs, divmod_limbs, extract_as_integer, multiply_limbs, or_left_shifted_into, top_set_bit,
};
use crate::rounding::{round_finite_to_precision, RoundingMode};
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::propagate_nan2;

impl BigFloat {
    /// IEEE 754-2019 §5.3.1 `remainder(self, other)`.
    ///
    /// Returns `self - n·other` where `n` is `self / other` rounded to
    /// the nearest integer, ties to even. The result is exact, with
    /// magnitude at most `|other|/2`, at a precision of
    /// `max(self.precision, other.precision)`.
    ///
    /// Special cases follow IEEE 754-2019: `remainder(±∞, y)` and
    /// `remainder(x, ±0)` are quiet NaN with `INVALID`;
    /// `remainder(x, ±∞) = x` for finite `x`; `remainder(±0, y) = ±0`
    /// for `y ≠ 0`; NaN operands propagate.
    #[must_use]
    pub fn remainder(&self, other: &Self) -> (Self, Status) {
        remainder_kernel(self, other)
    }
}

fn remainder_kernel(x: &BigFloat, y: &BigFloat) -> (BigFloat, Status) {
    let target = x.precision.max(y.precision);

    if let Some(propagated) = propagate_nan2(x, y, target) {
        return propagated;
    }

    // remainder(±∞, y): invalid. (Inf/Inf and Inf/NaN already handled.)
    if x.is_infinite() {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target, &[]).expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }
    // remainder(x, ±0): invalid (covers 0 % 0 too, x finite here).
    if y.is_zero() {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target, &[]).expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }
    // remainder(±0, y) = ±0 for finite-or-infinite y ≠ 0.
    if x.is_zero() {
        let z = BigFloat::try_new_zero(x.sign(), target).expect("precision >= 1");
        return (z, Status::OK);
    }
    // remainder(x, ±∞) = x for finite nonzero x.
    if y.is_infinite() {
        let (widened, _) = x
            .round_to_precision(target, RoundingMode::NearestEven)
            .expect("target >= 1");
        return (widened, Status::OK);
    }

    remainder_finite(x, y, target)
}

fn remainder_finite(x: &BigFloat, y: &BigFloat, target: u32) -> (BigFloat, Status) {
    let (ex, mx, px) = match &x.class {
        Class::Normal {
            exponent, mantissa, ..
        } => (*exponent, mantissa.as_slice(), x.precision),
        _ => unreachable!("remainder_finite: x is Normal"),
    };
    let (ey, my, py) = match &y.class {
        Class::Normal {
            exponent, mantissa, ..
        } => (*exponent, mantissa.as_slice(), y.precision),
        _ => unreachable!("remainder_finite: y is Normal"),
    };
    let sx = x.sign();

    // Early exit: 2|x| < |y| ⟹ the nearest integer quotient is 0, so
    // the remainder is x itself. The MSB of |x| sits at `ex`, of |y|
    // at `ey`, so 2|x| has its MSB at `ex + 1`. This also bounds the
    // `ay - ax` shift in the y-dominant branch below.
    if i128::from(ex) + 1 < i128::from(ey) {
        let (widened, _) = x
            .round_to_precision(target, RoundingMode::NearestEven)
            .expect("target >= 1");
        return (widened, Status::OK);
    }

    let mx_int = extract_as_integer(mx, px);
    let my_int = extract_as_integer(my, py);
    let ax = i128::from(ex) - i128::from(px) + 1;
    let ay = i128::from(ey) - i128::from(py) + 1;

    // R: the truncated (fmod-style) remainder magnitude as an integer
    // at scale 2^s; y_scaled: |y| at the same scale; quo_parity_odd:
    // parity of the truncated quotient, needed only on a round-to-even
    // tie. `s` is the common bottom scale min(ax, ay).
    let (r, y_scaled, s, quo_parity): (Vec<u64>, Vec<u64>, i128, QuotientParity) = if ax >= ay {
        // x-dominant: s = ay, divide Mx·2^(ax-ay) by My.
        let k = (ax - ay) as u128;
        if k == 0 {
            let (q, r) = divmod_limbs(&mx_int, &my_int);
            (r, my_int.clone(), ay, QuotientParity::Known(low_bit(&q)))
        } else {
            let r = mul_mod(&mod_reduce(&mx_int, &my_int), &modpow2(k, &my_int), &my_int);
            // Defer the parity computation; it costs a second modular
            // exponentiation and is consulted only on an exact tie.
            (r, my_int.clone(), ay, QuotientParity::Deferred { k })
        }
    } else {
        // y-dominant: s = ax, divide Mx by My·2^(ay-ax). The early
        // exit above guarantees ay - ax <= px - py + 1, so the shifted
        // divisor is bounded (at most px + 1 bits) and safe to form.
        let d = (ay - ax) as u32;
        let y_scaled = shl_bits(&my_int, py, d);
        let (q, r) = divmod_limbs(&mx_int, &y_scaled);
        (r, y_scaled, ax, QuotientParity::Known(low_bit(&q)))
    };

    // Round-to-nearest-even adjustment: compare 2R with the scaled |y|.
    let two_r = multiply_limbs(&r, &[2]);
    let subtract = match cmp_limbs(&two_r, &y_scaled) {
        core::cmp::Ordering::Greater => true,
        core::cmp::Ordering::Less => false,
        core::cmp::Ordering::Equal => match quo_parity {
            QuotientParity::Known(odd) => odd,
            QuotientParity::Deferred { k } => quotient_is_odd(&mx_int, &my_int, k, &r),
        },
    };

    let (mag, result_sign) = if subtract {
        // remainder = |x| - (q+1)|y| < 0: magnitude |y|_scaled - R,
        // and the sign of the signed remainder flips relative to x.
        (sub_limbs(&y_scaled, &r), sx.flip())
    } else {
        (r, sx)
    };

    build_scaled(result_sign, &mag, s, target, sx)
}

/// Quotient parity for the tie case: known outright from a direct
/// division, or deferred (recovered with a second modular
/// exponentiation modulo `2·My`) when the dividend was too large to
/// form.
enum QuotientParity {
    Known(bool),
    Deferred { k: u128 },
}

/// `quotient_is_odd` recovers the parity of `q = (Mx·2^k) div My` on a
/// tie. `X mod 2My = (q mod 2)·My + R`, so `(X mod 2My) - R` is `0`
/// (q even) or `My` (q odd).
fn quotient_is_odd(mx_int: &[u64], my_int: &[u64], k: u128, r: &[u64]) -> bool {
    let two_my = multiply_limbs(my_int, &[2]);
    let x_mod = mul_mod(&mod_reduce(mx_int, &two_my), &modpow2(k, &two_my), &two_my);
    let diff = sub_limbs(&x_mod, r);
    !is_zero(&diff)
}

/// Build the exact result `sign · mag · 2^s` at `target` precision.
/// `zero_sign` is the sign of `x`, used for the signed-zero result of
/// an exact multiple (IEEE: `remainder` of a multiple is `±0` with the
/// dividend's sign).
fn build_scaled(
    sign: Sign,
    mag: &[u64],
    s: i128,
    target: u32,
    zero_sign: Sign,
) -> (BigFloat, Status) {
    let Some(top_bit) = top_set_bit(mag) else {
        let z = BigFloat::try_new_zero(zero_sign, target).expect("precision >= 1");
        return (z, Status::OK);
    };

    let intermediate_precision = (top_bit + 1) as u32;
    let intermediate_limbs = limbs_for(intermediate_precision);
    let mut intermediate = vec![0u64; intermediate_limbs];
    let dst_low_zero =
        ((intermediate_limbs as u64) * 64 - u64::from(intermediate_precision)) as u32;
    or_left_shifted_into(&mut intermediate, mag, intermediate_precision, dst_low_zero);

    // Result MSB sits at bit position s + top_bit. The value is exact
    // (it fits in <= max(px, py) <= target bits), but saturate the
    // exponent into i64 for totality, the same contract div uses.
    let result_exp_wide = s + i128::from(top_bit as i64);
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
        sign,
        result_exp,
        &intermediate,
        intermediate_precision,
        false,
        target,
        RoundingMode::NearestEven,
    );
    let status = round_status | exp_saturation;
    auto_raise(status);
    (value, status)
}

// ----- limb helpers local to the remainder kernel -----

fn is_zero(v: &[u64]) -> bool {
    v.iter().all(|&l| l == 0)
}

/// Low bit of a little-endian limb integer (its parity).
fn low_bit(v: &[u64]) -> bool {
    v.first().is_some_and(|&l| l & 1 == 1)
}

/// `a mod m` via the shared limb division.
fn mod_reduce(a: &[u64], m: &[u64]) -> Vec<u64> {
    divmod_limbs(a, m).1
}

/// `(a · b) mod m`.
fn mul_mod(a: &[u64], b: &[u64], m: &[u64]) -> Vec<u64> {
    let product = multiply_limbs(a, b);
    divmod_limbs(&product, m).1
}

/// `2^k mod m`, left-to-right square-and-multiply. `O(log k)`
/// multiplications, so a huge `k` (up to the exponent range) is cheap.
fn modpow2(k: u128, m: &[u64]) -> Vec<u64> {
    let mut result = mod_reduce(&[1u64], m);
    if k == 0 {
        return result;
    }
    let base = mod_reduce(&[2u64], m);
    let nbits = 128 - k.leading_zeros();
    for i in (0..nbits).rev() {
        result = mul_mod(&result, &result, m);
        if (k >> i) & 1 == 1 {
            result = mul_mod(&result, &base, m);
        }
    }
    result
}

/// `v << bits` as a fresh little-endian integer. `v` is `value_bits`
/// wide; the caller bounds `bits` so the result stays small.
fn shl_bits(v: &[u64], value_bits: u32, bits: u32) -> Vec<u64> {
    let out_bits = value_bits + bits;
    let mut out = vec![0u64; limbs_for(out_bits.max(1))];
    or_left_shifted_into(&mut out, v, value_bits, bits);
    out
}

/// `a - b` for `a >= b`, little-endian, with borrow propagation.
fn sub_limbs(a: &[u64], b: &[u64]) -> Vec<u64> {
    let n = a.len().max(b.len());
    let mut out = vec![0u64; n];
    let mut borrow: i128 = 0;
    for (i, slot) in out.iter_mut().enumerate() {
        let av = i128::from(a.get(i).copied().unwrap_or(0));
        let bv = i128::from(b.get(i).copied().unwrap_or(0));
        let mut d = av - bv - borrow;
        if d < 0 {
            d += 1i128 << 64;
            borrow = 1;
        } else {
            borrow = 0;
        }
        *slot = d as u64;
    }
    debug_assert_eq!(borrow, 0, "sub_limbs requires a >= b");
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bf(n: i64) -> BigFloat {
        BigFloat::try_from_i64_exact(n, 53).expect("precision >= 1")
    }

    fn eq(a: &BigFloat, b: &BigFloat) -> bool {
        matches!(a.partial_cmp(b).0, Some(core::cmp::Ordering::Equal))
    }

    #[test]
    fn integer_values() {
        // 5/3 → 2, 5 - 6 = -1.
        assert!(eq(&bf(5).remainder(&bf(3)).0, &bf(-1)));
        // 7/3 → 2, 7 - 6 = 1.
        assert!(eq(&bf(7).remainder(&bf(3)).0, &bf(1)));
        // Exact multiple → 0.
        assert!(eq(&bf(6).remainder(&bf(3)).0, &bf(0)));
    }

    #[test]
    fn round_half_to_even_ties() {
        // 5/2 = 2.5 → 2 (even), 5 - 4 = 1.
        assert!(eq(&bf(5).remainder(&bf(2)).0, &bf(1)));
        // 3/2 = 1.5 → 2 (even), 3 - 4 = -1.
        assert!(eq(&bf(3).remainder(&bf(2)).0, &bf(-1)));
        // 7/2 = 3.5 → 4 (even), 7 - 8 = -1.
        assert!(eq(&bf(7).remainder(&bf(2)).0, &bf(-1)));
        // 1/2 = 0.5 → 0 (even), 1 - 0 = 1.
        assert!(eq(&bf(1).remainder(&bf(2)).0, &bf(1)));
    }

    #[test]
    fn sign_rules() {
        // remainder is odd in the dividend, independent of divisor sign.
        assert!(eq(&bf(-5).remainder(&bf(3)).0, &bf(1)));
        assert!(eq(&bf(5).remainder(&bf(-3)).0, &bf(-1)));
        assert!(eq(&bf(-5).remainder(&bf(-3)).0, &bf(1)));
    }

    #[test]
    fn special_cases() {
        let zero = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let neg_zero = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let inf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();

        // remainder(x, ±0) = qNaN + INVALID.
        let (r, s) = bf(5).remainder(&zero);
        assert!(r.is_nan() && s == Status::INVALID);
        // remainder(±∞, y) = qNaN + INVALID.
        let (r, s) = inf.remainder(&bf(3));
        assert!(r.is_nan() && s == Status::INVALID);
        // remainder(x, ±∞) = x.
        assert!(eq(&bf(5).remainder(&inf).0, &bf(5)));
        // remainder(±0, y) = ±0 (sign of x).
        let (r, _) = zero.remainder(&bf(3));
        assert!(r.is_zero() && !r.is_sign_negative());
        let (r, _) = neg_zero.remainder(&bf(3));
        assert!(r.is_zero() && r.is_sign_negative());
        // NaN propagates from either operand.
        assert!(nan.remainder(&bf(3)).0.is_nan());
        assert!(bf(3).remainder(&nan).0.is_nan());
    }

    #[cfg(feature = "ops")]
    #[test]
    fn operator_matches_method() {
        // The `%` overload is the IEEE remainder.
        let v = bf(7) % bf(3);
        assert!(eq(&v, &bf(7).remainder(&bf(3)).0));
    }
}
