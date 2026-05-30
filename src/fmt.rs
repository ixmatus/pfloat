//! Decimal formatting for [`BigFloat`] and [`FixedFloat<PREC>`].
//!
//! `to_decimal_string(digits, mode)` produces a decimal
//! representation with `digits` significant digits, rounded to
//! nearest under the given mode. [`Display`] uses
//! `round_trip_digit_count(precision)` digits — enough that
//! [`BigFloat::parse_str`] at the same precision recovers the
//! exact value.
//!
//! The algorithm:
//!
//! 1. Estimate `decimal_exp ≈ floor(log10(|v|))` from the binary
//!    exponent using the rational approximation
//!    `log10(2) ≈ 30103/100000`.
//! 2. Compute `scaled = round(v × 10^(N-1-decimal_exp))` as a
//!    multi-precision integer. This is done by expressing the
//!    target as `numerator / denominator` where
//!    `numerator = m × 5^p5_num × 2^p2_num` and
//!    `denominator = 5^p5_den × 2^p2_den`, with the p5/p2 exponents
//!    derived from the signs of `scale = e - p + 1` and `shift =
//!    N - 1 - decimal_exp`.
//! 3. Divmod `numerator / denominator`; round per the user's
//!    rounding mode.
//! 4. Convert `scaled` to decimal digits via repeated `divmod 10`.
//! 5. Verify that `scaled` has exactly `N` digits; if off by one
//!    (rare: the log10 estimate's residual), adjust `decimal_exp`
//!    and retry.
//! 6. Compose either a fixed-point or scientific format depending
//!    on `decimal_exp`'s magnitude.
//!
//! Slice 2b ships fixed-digit-count formatting. Shortest
//! round-trip (Dragon4 / Steele-White) can be a follow-up if and
//! when callers want it; the current Display output rounds-trips
//! at the operand's own precision, which is what most callers want.

use alloc::string::String;
use alloc::vec;
use alloc::vec::Vec;
use core::cmp::Ordering;
use core::fmt::{self, Display, Formatter, Write};

use crate::big::BigFloat;
use crate::class::Class;
use crate::mantissa::limbs_for;
use crate::ops::limbs::{
    cmp_limbs, divmod_limbs, extract_as_integer, limbs_add_assign, multiply_limbs,
    or_left_shifted_into, top_set_bit,
};
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;

/// The decimal-exponent magnitude past which formatting saturates
/// instead of rendering exact digits (ADR-0051).
///
/// Rendering `m · 2^E` exactly needs a `5^|shift|` bridge whose size is
/// linear in `E` (`10^d = 2^d · 5^d`, and the 5-part does not cancel),
/// so a finite value with an astronomically large binary exponent has no
/// bounded-cost decimal rendering — and the bounded scientific shortcut
/// (`E · log10(2)` for the exponent, `exp10(frac)` for the digits) is
/// unavailable, because `log10`/`exp10` are `exp-log`-gated while this
/// module builds under `big` alone. The cap is value-matched to parse's
/// `MAX_DECIMAL_EXPONENT` (`src/parse.rs`) so the exactly-renderable
/// range is exactly the range parse round-trips; past it a value
/// saturates, mirroring parse's own `±inf` / `±0` saturation.
const MAX_FORMAT_DECIMAL_EXPONENT: i64 = 1_000_000;

/// How a finite value past [`MAX_FORMAT_DECIMAL_EXPONENT`] renders: too
/// large to print reads back as `inf`, too small as `0`, matching parse.
#[derive(Clone, Copy)]
enum Saturation {
    Infinite,
    Zero,
}

/// Render a sign-carrying saturation token (`inf` / `-inf` / `0` / `-0`)
/// for a finite value whose magnitude is past the format cap. The token
/// is a resource bound, not an arithmetic claim that the value is
/// infinite or zero; it is the value parse would itself produce for any
/// decimal that large or small (ADR-0051).
fn saturated_string(sign: Sign, kind: Saturation) -> String {
    let mut s = String::new();
    if matches!(sign, Sign::Negative) {
        s.push('-');
    }
    s.push_str(match kind {
        Saturation::Infinite => "inf",
        Saturation::Zero => "0",
    });
    s
}

impl BigFloat {
    /// Returns a decimal string with `digits` significant digits,
    /// rounded under `mode`.
    ///
    /// Special values render as `"nan"`, `"-nan"`, `"inf"`,
    /// `"-inf"`, `"0"`, `"-0"`. Finite normals use fixed-point
    /// formatting when the decimal exponent is in `[-4, 15]` and
    /// scientific notation otherwise.
    #[must_use]
    pub fn to_decimal_string(&self, digits: u32, mode: RoundingMode) -> String {
        debug_assert!(digits >= 1, "digit count must be >= 1");
        match &self.class {
            Class::Nan { sign, .. } => {
                let mut s = String::new();
                if matches!(sign, Sign::Negative) {
                    s.push('-');
                }
                s.push_str("nan");
                s
            }
            Class::Infinity { sign } => {
                let mut s = String::new();
                if matches!(sign, Sign::Negative) {
                    s.push('-');
                }
                s.push_str("inf");
                s
            }
            Class::Zero { sign } => {
                let mut s = String::new();
                if matches!(sign, Sign::Negative) {
                    s.push('-');
                }
                s.push('0');
                s
            }
            Class::Normal {
                sign,
                exponent,
                mantissa,
            } => format_normal(mantissa, self.precision, *exponent, digits, *sign, mode),
        }
    }

    /// Number of decimal digits required for a round-trip-safe
    /// representation at the given precision:
    /// `ceil(p × log10(2)) + 1`. For `p = 53` this is 17 (the f64
    /// round-trip), for `p = 113` it is 35.
    #[inline]
    #[must_use]
    pub fn round_trip_digit_count(precision: u32) -> u32 {
        // p × 30103 / 100000, rounded up, then +1.
        let approx = (u64::from(precision) * 30103).div_ceil(100_000);
        (approx as u32) + 1
    }
}

impl Display for BigFloat {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let digits = Self::round_trip_digit_count(self.precision);
        let s = self.to_decimal_string(digits, RoundingMode::NearestEven);
        f.write_str(&s)
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// Returns a decimal string with `digits` significant digits.
    /// Delegates to [`BigFloat::to_decimal_string`].
    #[must_use]
    pub fn to_decimal_string(&self, digits: u32, mode: RoundingMode) -> String {
        self.to_big().to_decimal_string(digits, mode)
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> Display for FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let digits = BigFloat::round_trip_digit_count(PREC);
        let s = self.to_decimal_string(digits, RoundingMode::NearestEven);
        f.write_str(&s)
    }
}

fn format_normal(
    mantissa: &[u64],
    precision: u32,
    exponent: i64,
    digit_count: u32,
    sign: Sign,
    mode: RoundingMode,
) -> String {
    let m_int = extract_as_integer(mantissa, precision);
    // Saturating exponent arithmetic: `exponent` can sit at i64::MIN/MAX
    // for a value built by repeated squaring (mul saturates the exponent;
    // pfloat has no emax), so the bare `exponent - precision + 1` could
    // overflow. Within the cap below the result is far from the bounds,
    // so saturation never changes a renderable value's scale.
    let scale = exponent
        .saturating_sub(i64::from(precision))
        .saturating_add(1);

    // Magnitude cap (ADR-0051): past MAX_FORMAT_DECIMAL_EXPONENT a finite
    // value has no bounded-cost decimal rendering, so it saturates like
    // parse. For a normalized mantissa `top_set_bit == precision - 1`, so
    // `log2_value == exponent`; the rational log10 estimate is exact to
    // far better than one part in the cap at this magnitude.
    let log2_value = top_set_bit(&m_int)
        .map_or(0, |t| t as i64)
        .saturating_add(scale);
    let decimal_exp_estimate = approximate_log10_floor(log2_value);
    if decimal_exp_estimate > MAX_FORMAT_DECIMAL_EXPONENT {
        return saturated_string(sign, Saturation::Infinite);
    }
    if decimal_exp_estimate < -MAX_FORMAT_DECIMAL_EXPONENT {
        return saturated_string(sign, Saturation::Zero);
    }

    let (digits, decimal_exp) = extract_digits(&m_int, scale, digit_count, mode, sign);
    compose(&digits, decimal_exp, sign)
}

/// Returns `(digits, decimal_exp)` where `digits` is a Vec of bytes
/// in `0..=9` of length exactly `digit_count`, and the represented
/// value is approximately
/// `digits[0].digits[1..] × 10^decimal_exp`.
fn extract_digits(
    m_int: &[u64],
    scale: i64,
    digit_count: u32,
    mode: RoundingMode,
    sign: Sign,
) -> (Vec<u8>, i64) {
    let top_bit_m = top_set_bit(m_int).expect("non-zero m");
    let log2_value = (top_bit_m as i64) + scale;

    // log10(value) ≈ log2_value × log10(2). Use the rational
    // approximation 30103/100000; this is exact to one part in
    // 10^4 over any practical range. Float the residual into the
    // adjustment loop below.
    let mut decimal_exp = approximate_log10_floor(log2_value);

    // Up to ~3 iterations: the estimate is usually correct or
    // off by 1 in either direction. We cap the loop just in case
    // of pathological pathway.
    for _ in 0..10 {
        let shift = i64::from(digit_count) - 1 - decimal_exp;
        let scaled = compute_scaled(m_int, scale, shift, mode, sign);
        let mut digits = int_to_decimal(&scaled);
        let observed = digits.len() as i64;
        let target = i64::from(digit_count);
        if observed == target {
            return (digits, decimal_exp);
        } else if observed == target + 1 {
            // One extra leading digit; bump decimal_exp by the
            // overshoot and retry. (For example: estimate gave
            // dec_exp = 5 but the true value is 6.)
            decimal_exp += 1;
        } else if observed + 1 == target {
            decimal_exp -= 1;
        } else if observed > target {
            // Far off: jump directly.
            decimal_exp += observed - target;
        } else if observed < target {
            // Pad: extend digits with zeros and adjust decimal_exp.
            // This is the "value rounded down past a power of 10"
            // case (e.g., 9.9999 → 10.000 with digit count 5).
            let padding = (target - observed) as usize;
            digits.extend(core::iter::repeat_n(0u8, padding));
            return (digits, decimal_exp);
        }
    }
    // Fallback: shouldn't be reached. Return what we have.
    let shift = i64::from(digit_count) - 1 - decimal_exp;
    let scaled = compute_scaled(m_int, scale, shift, mode, sign);
    let digits = int_to_decimal(&scaled);
    (digits, decimal_exp)
}

/// Approximate `floor(n × log10(2))` using the integer rational
/// approximation `log10(2) ≈ 30103/100000`. Handles negative `n`
/// with proper floor semantics (Rust's `/` truncates toward zero).
fn approximate_log10_floor(n: i64) -> i64 {
    let num = n.saturating_mul(30103);
    if num >= 0 {
        num / 100_000
    } else {
        // For negative numerator, integer division truncates
        // toward zero (i.e., ceil for negative). We want floor.
        let q = num / 100_000;
        if num % 100_000 == 0 {
            q
        } else {
            q - 1
        }
    }
}

/// Compute `round(value × 10^shift)` where `value = m × 2^scale`,
/// under the chosen rounding mode and sign.
///
/// The caller gates this behind [`MAX_FORMAT_DECIMAL_EXPONENT`]
/// (`format_normal`), which bounds `|shift|` and `|scale + shift|` well
/// within `u32` for any normal-precision operand, so the power-of-two
/// bit counts below never approach `u32::MAX`. The `saturating_add`s and
/// the saturating `try_into`s are a defense-in-depth backstop for the
/// astronomically-high-precision corner (a multi-hundred-MB mantissa),
/// where the result is bounded-but-truncated rather than a panic.
fn compute_scaled(
    m_int: &[u64],
    scale: i64,
    shift: i64,
    mode: RoundingMode,
    sign: Sign,
) -> Vec<u64> {
    let combined = scale.saturating_add(shift);
    let num_p2 = combined.max(0).try_into().unwrap_or(u32::MAX);
    let den_p2 = (-combined).max(0).try_into().unwrap_or(u32::MAX);
    let num_p5 = shift.max(0).try_into().unwrap_or(u32::MAX);
    let den_p5 = (-shift).max(0).try_into().unwrap_or(u32::MAX);

    let mut num: Vec<u64> = m_int.to_vec();
    if num_p5 > 0 {
        num = multiply_limbs(&num, &pow5(num_p5));
    }
    if num_p2 > 0 {
        let bits = top_set_bit(&num).map_or(0u32, |t| (t + 1) as u32);
        let total = bits.saturating_add(num_p2);
        let mut shifted = vec![0u64; limbs_for(total)];
        or_left_shifted_into(&mut shifted, &num, bits, num_p2);
        num = shifted;
    }

    let mut den: Vec<u64> = if den_p5 > 0 { pow5(den_p5) } else { vec![1] };
    if den_p2 > 0 {
        let bits = top_set_bit(&den).map_or(0u32, |t| (t + 1) as u32);
        let total = bits.saturating_add(den_p2);
        let mut shifted = vec![0u64; limbs_for(total)];
        or_left_shifted_into(&mut shifted, &den, bits, den_p2);
        den = shifted;
    }

    let (quotient, remainder) = divmod_limbs(&num, &den);

    if should_round_up(&quotient, &remainder, &den, mode, sign) {
        increment_owned(quotient)
    } else {
        quotient
    }
}

fn should_round_up(
    quotient: &[u64],
    remainder: &[u64],
    divisor: &[u64],
    mode: RoundingMode,
    sign: Sign,
) -> bool {
    if remainder.iter().all(|&l| l == 0) {
        return false; // exact, no rounding needed
    }
    // Compare 2 × remainder with divisor.
    let two_r = multiply_limbs(remainder, &[2]);
    let cmp = cmp_limbs(&two_r, divisor);
    match mode {
        RoundingMode::TowardZero => false,
        RoundingMode::TowardPositive => matches!(sign, Sign::Positive),
        RoundingMode::TowardNegative => matches!(sign, Sign::Negative),
        RoundingMode::NearestAway => !matches!(cmp, Ordering::Less),
        RoundingMode::NearestEven => match cmp {
            Ordering::Less => false,
            Ordering::Greater => true,
            Ordering::Equal => quotient.first().copied().unwrap_or(0) & 1 == 1,
        },
    }
}

fn increment_owned(mut v: Vec<u64>) -> Vec<u64> {
    if v.is_empty() {
        v.push(1);
        return v;
    }
    let carried = limbs_add_assign(&mut v, &[1]);
    if carried {
        v.push(1);
    }
    v
}

/// Repeated divmod by 10 to extract decimal digits, most
/// significant first.
fn int_to_decimal(value: &[u64]) -> Vec<u8> {
    if value.iter().all(|&l| l == 0) {
        return vec![0];
    }
    let mut digits: Vec<u8> = Vec::new();
    let mut v: Vec<u64> = value.to_vec();
    while !v.iter().all(|&l| l == 0) {
        let (q, r) = divmod_limbs(&v, &[10u64]);
        digits.push(r.first().copied().unwrap_or(0) as u8);
        v = q;
        // Trim trailing zeros to keep divmod cost down.
        while v.len() > 1 && *v.last().unwrap() == 0 {
            v.pop();
        }
    }
    digits.reverse();
    digits
}

/// Compute `5^exp` as a multi-precision integer. Shared with the
/// parse module's `pow5`. Slice 2b duplicates the helper here to
/// avoid making `pow5` a `pub(crate)` surface in `parse.rs`.
fn pow5(exp: u32) -> Vec<u64> {
    if exp == 0 {
        return vec![1];
    }
    let mut base: Vec<u64> = vec![5];
    let mut result: Vec<u64> = vec![1];
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            result = multiply_limbs(&result, &base);
            trim_high_zeros(&mut result);
        }
        e >>= 1;
        if e > 0 {
            base = multiply_limbs(&base, &base);
            trim_high_zeros(&mut base);
        }
    }
    result
}

fn trim_high_zeros(v: &mut Vec<u64>) {
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
}

/// Choose between fixed-point and scientific notation based on the
/// decimal exponent, then compose the string.
fn compose(digits: &[u8], decimal_exp: i64, sign: Sign) -> String {
    let use_fixed = (-4..=15).contains(&decimal_exp);
    if use_fixed {
        compose_fixed(digits, decimal_exp, sign)
    } else {
        compose_scientific(digits, decimal_exp, sign)
    }
}

fn compose_fixed(digits: &[u8], decimal_exp: i64, sign: Sign) -> String {
    let mut s = String::new();
    if matches!(sign, Sign::Negative) {
        s.push('-');
    }
    if decimal_exp >= 0 {
        let int_count = (decimal_exp + 1) as usize;
        if int_count >= digits.len() {
            // All integer; pad with trailing zeros.
            for &d in digits {
                s.push(char::from(b'0' + d));
            }
            for _ in 0..(int_count - digits.len()) {
                s.push('0');
            }
        } else {
            for &d in &digits[..int_count] {
                s.push(char::from(b'0' + d));
            }
            s.push('.');
            for &d in &digits[int_count..] {
                s.push(char::from(b'0' + d));
            }
            trim_trailing_zeros_after_point(&mut s);
        }
    } else {
        // Value < 1: "0.0...0<digits>"
        s.push_str("0.");
        for _ in 0..((-decimal_exp) - 1) {
            s.push('0');
        }
        for &d in digits {
            s.push(char::from(b'0' + d));
        }
        trim_trailing_zeros_after_point(&mut s);
    }
    s
}

fn compose_scientific(digits: &[u8], decimal_exp: i64, sign: Sign) -> String {
    let mut s = String::new();
    if matches!(sign, Sign::Negative) {
        s.push('-');
    }
    s.push(char::from(b'0' + digits[0]));
    if digits.len() > 1 {
        s.push('.');
        for &d in &digits[1..] {
            s.push(char::from(b'0' + d));
        }
        trim_trailing_zeros_after_point(&mut s);
    }
    let _ = write!(s, "e{decimal_exp}");
    s
}

fn trim_trailing_zeros_after_point(s: &mut String) {
    if !s.contains('.') {
        return;
    }
    while s.ends_with('0') {
        s.pop();
    }
    if s.ends_with('.') {
        s.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::string::ToString;

    fn parse(s: &str) -> BigFloat {
        BigFloat::parse_str(s, 53, RoundingMode::NearestEven)
            .expect("parse should succeed")
            .0
    }

    #[test]
    fn format_zero() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        assert_eq!(pz.to_string(), "0");
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        assert_eq!(nz.to_string(), "-0");
    }

    #[test]
    fn format_inf_nan() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        assert_eq!(pi.to_string(), "inf");
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        assert_eq!(ni.to_string(), "-inf");
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        assert_eq!(q.to_string(), "nan");
        let nn = BigFloat::try_new_quiet_nan(Sign::Negative, 53, &[]).unwrap();
        assert_eq!(nn.to_string(), "-nan");
    }

    #[test]
    fn format_small_integers() {
        for n in [0i64, 1, 2, 3, 5, 7, 10, 42, 100, 12345] {
            let v = BigFloat::try_from_i64_exact(n, 53).unwrap();
            let formatted = v.to_string();
            let parsed = parse(&formatted);
            assert_eq!(
                parsed.partial_cmp(&v).0,
                Some(Ordering::Equal),
                "round-trip failed: {n} formatted as {formatted}"
            );
        }
    }

    #[test]
    fn format_negative_integers() {
        for n in [-1i64, -2, -42, -1000] {
            let v = BigFloat::try_from_i64_exact(n, 53).unwrap();
            let s = v.to_string();
            assert!(s.starts_with('-'), "negative should start with '-': {s}");
        }
    }

    #[test]
    fn format_exact_half() {
        // 0.5 is exact at any precision.
        let half = parse("0.5");
        let s = half.to_string();
        // Parse it back; should be exact.
        let back = parse(&s);
        assert_eq!(back.partial_cmp(&half).0, Some(Ordering::Equal));
    }

    #[test]
    fn format_round_trip_one_third() {
        // 1/3 at precision 53; format and re-parse should round-trip.
        let third = parse("0.333333333333333");
        let s = third.to_string();
        let back = BigFloat::parse_str(&s, 53, RoundingMode::NearestEven)
            .unwrap()
            .0;
        assert_eq!(back.partial_cmp(&third).0, Some(Ordering::Equal));
    }

    #[test]
    fn format_scientific_for_large() {
        // 1e20 should use scientific notation.
        let big = parse("1e20");
        let s = big.to_string();
        assert!(s.contains('e'), "expected scientific notation, got {s}");
    }

    #[test]
    fn format_scientific_for_tiny() {
        // 1e-10 should use scientific notation.
        let tiny = parse("1e-10");
        let s = tiny.to_string();
        assert!(s.contains('e'), "expected scientific notation, got {s}");
    }

    #[test]
    fn format_fixed_for_moderate() {
        // 0.001 (exponent ~ -3): fixed-point.
        let v = parse("0.001");
        let s = v.to_string();
        assert!(!s.contains('e'), "expected fixed-point, got {s}");
        assert!(s.starts_with("0.001"));
    }

    #[test]
    fn round_trip_digit_count_table() {
        // ceil(p × log10(2)) + 1 using log10(2) ≈ 30103/100000.
        assert_eq!(BigFloat::round_trip_digit_count(53), 17);
        assert_eq!(BigFloat::round_trip_digit_count(113), 36);
        assert_eq!(BigFloat::round_trip_digit_count(256), 79);
    }

    #[test]
    fn format_via_to_decimal_string() {
        let v = parse("1.5");
        let s = v.to_decimal_string(5, RoundingMode::NearestEven);
        // Either "1.5000" trimmed → "1.5", or some scientific form.
        let back = parse(&s);
        assert_eq!(back.partial_cmp(&v).0, Some(Ordering::Equal));
    }

    #[test]
    fn max_format_decimal_exponent_matches_parse_cap() {
        // Value-matched to parse's `MAX_DECIMAL_EXPONENT` (private to
        // `parse.rs`) so the renderable range equals the round-trippable
        // range. The literal is asserted here so a change to one side is
        // noticed against the other (ADR-0051).
        assert_eq!(MAX_FORMAT_DECIMAL_EXPONENT, 1_000_000);
    }

    #[test]
    fn format_cap_saturates_finite_huge_to_inf_and_zero() {
        // 2^(2^40): finite (exponent ≈ 1.1e12 ≪ i64::MAX), decimal
        // exponent ≈ 3.3e11 — far past the 1e6 cap. Must saturate to a
        // bounded token, never panic or OOM (regression for finding 11).
        let ne = RoundingMode::NearestEven;
        let mut x = BigFloat::try_from_i64_exact(2, 53).unwrap();
        for _ in 0..40 {
            x = x.mul(&x, ne).0;
        }
        assert!(!x.is_infinite(), "2^(2^40) is finite, not Inf");
        assert_eq!(x.to_decimal_string(17, ne), "inf");
        assert_eq!(x.to_string(), "inf");
        assert_eq!(x.negated().to_string(), "-inf");

        // 2^(-2^40): finite nonzero, decimal exponent ≈ -3.3e11.
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let tiny = one.div(&x, ne).0;
        assert!(!tiny.is_zero(), "2^(-2^40) is finite nonzero");
        assert_eq!(tiny.to_string(), "0");
        assert_eq!(tiny.negated().to_string(), "-0");
    }

    #[test]
    fn format_cap_boundary_under_renders_over_saturates() {
        // Cheap boundary straddle via squaring: 2^(2^k) has decimal
        // exponent ≈ 0.30103·2^k. k = 21 → ≈ 6.3e5 (inside the cap), one
        // more squaring (k = 22) → ≈ 1.26e6 (past it).
        let ne = RoundingMode::NearestEven;
        let mut under = BigFloat::try_from_i64_exact(2, 53).unwrap();
        for _ in 0..21 {
            under = under.mul(&under, ne).0;
        }
        let s = under.to_string();
        assert!(s.contains('e'), "in-cap large value renders scientific: {s}");
        let back = BigFloat::parse_str(&s, 53, ne).unwrap().0;
        assert_eq!(
            back.partial_cmp(&under).0,
            Some(Ordering::Equal),
            "in-cap render round-trips"
        );

        let over = under.mul(&under, ne).0; // 2^(2^22)
        assert_eq!(over.to_string(), "inf", "just past the cap saturates");
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixed_format_basic() {
        let v = FixedFloat::<53>::try_from_i64_exact(42).unwrap();
        let s = v.to_string();
        assert_eq!(s, "42");
    }

    #[test]
    fn parse_format_roundtrip_at_high_precision() {
        // Pick a few values whose binary representation is exact at
        // precision 113 (which can hold them as integers).
        for s in ["1", "-1", "2.5", "0.5", "0.25", "100", "1e10"] {
            let v = BigFloat::parse_str(s, 113, RoundingMode::NearestEven)
                .unwrap()
                .0;
            let formatted = v.to_string();
            let back = BigFloat::parse_str(&formatted, 113, RoundingMode::NearestEven)
                .unwrap()
                .0;
            assert_eq!(
                back.partial_cmp(&v).0,
                Some(Ordering::Equal),
                "round-trip failed: {s} → {formatted}",
            );
        }
    }

    #[test]
    fn display_large_exponent_round_trips() {
        // Regression for the libFuzzer `parse` out-of-memory (slice
        // parse-oom-divmod). A large decimal exponent drove the
        // bit-at-a-time `divmod_limbs` into O(exponent²) work and
        // gigabytes of transient allocation inside the Display digit
        // extraction (`compute_scaled`). Algorithm D makes that linear
        // in the operand size. The exponent here is small enough that
        // the test is fast, yet large enough that the old quadratic
        // routine took seconds per value. Covers the positive and
        // negative exponent paths and a non-power-of-ten mantissa.
        // The exponents stay well inside `MAX_DECIMAL_EXPONENT` so the
        // exact `pow5` itself is cheap and the test runs in well under
        // a second.
        for s in ["1e100000", "-1e100000", "1e-100000", "1.25e50000"] {
            let v = BigFloat::parse_str(s, 113, RoundingMode::NearestEven)
                .unwrap()
                .0;
            let formatted = v.to_string();
            let back = BigFloat::parse_str(&formatted, 113, RoundingMode::NearestEven)
                .unwrap()
                .0;
            assert_eq!(
                back.partial_cmp(&v).0,
                Some(Ordering::Equal),
                "round-trip failed: {s} → {formatted}",
            );
        }
    }

    #[test]
    fn pow5_basic() {
        assert_eq!(pow5(0), vec![1]);
        assert_eq!(pow5(1), vec![5]);
        assert_eq!(pow5(2), vec![25]);
        assert_eq!(pow5(5), vec![3125]);
    }

    #[test]
    fn int_to_decimal_basic() {
        assert_eq!(int_to_decimal(&[0]), vec![0]);
        assert_eq!(int_to_decimal(&[5]), vec![5]);
        assert_eq!(int_to_decimal(&[123]), vec![1, 2, 3]);
        assert_eq!(int_to_decimal(&[1_000_000]), vec![1, 0, 0, 0, 0, 0, 0]);
    }
}
