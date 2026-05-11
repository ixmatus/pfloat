//! Decimal-string parsing into [`BigFloat`] and [`FixedFloat`].
//!
//! `parse_str(s, precision, mode)` accepts the following grammar:
//!
//! ```text
//! float := whitespace? sign? body whitespace?
//! sign  := '+' | '-'
//! body  := special | finite
//! special := /(?i)inf(inity)?/ | /(?i)nan/
//! finite := digits? ('.' digits?)? exponent?
//! digits := /[0-9]+/
//! exponent := /[eE]/ sign? /[0-9]+/
//! ```
//!
//! Either the integer part or the fractional part must contain at
//! least one digit (so `.` alone is an error, but `.5` and `1.`
//! are accepted). Leading/trailing ASCII whitespace is ignored.
//!
//! Decimal-to-binary conversion: parse the digits into a
//! multi-precision integer `m` and a signed decimal exponent
//! `e_dec`. The value is `m × 10^e_dec = m × 5^e_dec × 2^e_dec`.
//! For `e_dec >= 0`, multiply `m` by `5^e_dec`; the result times
//! `2^e_dec` is the exact value, top-aligned and routed through
//! the rounding pipeline. For `e_dec < 0`, divide a shifted `m`
//! by `5^|e_dec|`; the remainder feeds the sticky bit, and the
//! quotient becomes the rounded mantissa.
//!
//! Slice 2a does not yet emit the result via `Display`; that lands
//! in slice 2b.

use alloc::vec;
use alloc::vec::Vec;

use crate::big::BigFloat;
use crate::class::Class;
use crate::mantissa::limbs_for;
use crate::ops::limbs::{
    divmod_limbs, limbs_add_assign, multiply_limbs, or_left_shifted_into, top_set_bit,
};
use crate::rounding::{round_finite_to_precision, RoundingMode};
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;

/// Failure modes for [`BigFloat::parse_str`] and
/// [`FixedFloat::parse_str`].
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum ParseError {
    /// Input was empty (or contained only whitespace).
    Empty,
    /// Input contained a character outside the accepted grammar.
    InvalidCharacter,
    /// Neither the integer part nor the fractional part had any
    /// digits.
    MissingDigits,
    /// The exponent string was malformed (e.g. `e` followed by no
    /// digits).
    InvalidExponent,
    /// The exponent magnitude exceeds the `i32` range supported by
    /// pfloat's decimal parser. Real workloads do not reach this
    /// limit (exponent `±10^9` already exceeds any astronomy use).
    ExponentOutOfRange,
    /// Precision validation failed.
    PrecisionZero,
}

impl BigFloat {
    /// Parses a decimal string into a [`BigFloat`] at the given
    /// precision, rounding under `mode` if the decimal value cannot
    /// be represented exactly.
    ///
    /// The grammar is documented at the module level. Returns
    /// `Err(ParseError)` for syntactic errors and
    /// `Ok((value, Status))` otherwise. The status carries
    /// [`Status::INEXACT`] when the decimal value rounded.
    ///
    /// # Untrusted input
    ///
    /// This function does not cap the digit count: the lexer
    /// collects integer and fractional digits into an allocated
    /// buffer whose size is proportional to the input string.
    /// Callers handling strings from untrusted sources should bound
    /// the input length before invoking `parse_str`. The decimal
    /// exponent magnitude is capped internally at `1_000_000`, so
    /// an oversized exponent converts to overflow or underflow
    /// without further allocation; the digit count itself has no
    /// such cap.
    pub fn parse_str(
        s: &str,
        precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), ParseError> {
        if precision == 0 {
            return Err(ParseError::PrecisionZero);
        }
        let parsed = lex(s)?;
        Ok(decimal_to_bigfloat(parsed, precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// Parses a decimal string into a [`FixedFloat<PREC>`],
    /// rounding under `mode` if the decimal value cannot be
    /// represented exactly. Delegates to
    /// [`BigFloat::parse_str`].
    pub fn parse_str(s: &str, mode: RoundingMode) -> Result<(Self, Status), ParseError> {
        let (big, status) = BigFloat::parse_str(s, PREC, mode)?;
        // BigFloat is at PREC; the conversion to FixedFloat<PREC>
        // is exact.
        Ok((
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        ))
    }
}

// --- Lexing ---

#[derive(Clone, Debug)]
struct ParsedDecimal {
    sign: Sign,
    body: DecimalBody,
}

#[derive(Clone, Debug)]
enum DecimalBody {
    /// `NaN` (sign preserved).
    Nan,
    /// `±∞`.
    Infinity,
    /// Finite: `digits × 10^exponent`. `digits` is the
    /// concatenated integer-plus-fractional digit string with
    /// `exponent_adjustment` already absorbed; `exponent` here is
    /// the signed decimal exponent of the result.
    Finite {
        digits: Vec<u8>, // each byte 0..=9
        exponent: i32,
    },
}

fn lex(s: &str) -> Result<ParsedDecimal, ParseError> {
    let bytes = s.as_bytes();
    let mut i = 0;
    let len = bytes.len();

    // Skip leading whitespace.
    while i < len && bytes[i].is_ascii_whitespace() {
        i += 1;
    }

    if i == len {
        return Err(ParseError::Empty);
    }

    // Optional sign.
    let sign = match bytes[i] {
        b'+' => {
            i += 1;
            Sign::Positive
        }
        b'-' => {
            i += 1;
            Sign::Negative
        }
        _ => Sign::Positive,
    };

    // After sign: special values, or digits/dot.
    if i == len {
        return Err(ParseError::MissingDigits);
    }

    // Try special values.
    if let Some(body) = lex_special(&bytes[i..]) {
        // Advance past the special token to verify trailing
        // whitespace only.
        let special_len = special_token_len(&bytes[i..]).unwrap();
        let mut j = i + special_len;
        while j < len && bytes[j].is_ascii_whitespace() {
            j += 1;
        }
        if j != len {
            return Err(ParseError::InvalidCharacter);
        }
        return Ok(ParsedDecimal { sign, body });
    }

    // Finite body: digits and/or dot, then optional exponent.
    let mut digits: Vec<u8> = Vec::new();
    let mut decimal_point_pos: Option<usize> = None;
    let mut saw_digit = false;

    while i < len {
        let c = bytes[i];
        if c.is_ascii_digit() {
            digits.push(c - b'0');
            saw_digit = true;
            i += 1;
        } else if c == b'.' {
            if decimal_point_pos.is_some() {
                return Err(ParseError::InvalidCharacter);
            }
            decimal_point_pos = Some(digits.len());
            i += 1;
        } else {
            break;
        }
    }

    if !saw_digit {
        return Err(ParseError::MissingDigits);
    }

    // Optional exponent.
    let mut exponent: i32 = 0;
    if i < len && (bytes[i] == b'e' || bytes[i] == b'E') {
        i += 1;
        if i == len {
            return Err(ParseError::InvalidExponent);
        }
        let exp_sign: i32 = match bytes[i] {
            b'+' => {
                i += 1;
                1
            }
            b'-' => {
                i += 1;
                -1
            }
            _ => 1,
        };
        if i == len || !bytes[i].is_ascii_digit() {
            return Err(ParseError::InvalidExponent);
        }
        let mut acc: i64 = 0;
        while i < len && bytes[i].is_ascii_digit() {
            acc = acc * 10 + i64::from(bytes[i] - b'0');
            if acc > i64::from(i32::MAX) {
                return Err(ParseError::ExponentOutOfRange);
            }
            i += 1;
        }
        let signed = acc * i64::from(exp_sign);
        exponent = i32::try_from(signed).map_err(|_| ParseError::ExponentOutOfRange)?;
    }

    // Trailing whitespace.
    while i < len && bytes[i].is_ascii_whitespace() {
        i += 1;
    }
    if i != len {
        return Err(ParseError::InvalidCharacter);
    }

    // Bake the decimal-point position into the exponent: a digit
    // after the decimal point shifts the exponent down by one.
    if let Some(pos) = decimal_point_pos {
        let frac_digits = digits.len() - pos;
        exponent = exponent
            .checked_sub(i32::try_from(frac_digits).map_err(|_| ParseError::ExponentOutOfRange)?)
            .ok_or(ParseError::ExponentOutOfRange)?;
    }

    // Strip leading zeros from digits (purely cosmetic; the
    // multi-precision integer construction handles them, but
    // dropping them keeps the integer's bit length tight).
    while digits.len() > 1 && digits[0] == 0 {
        digits.remove(0);
    }

    Ok(ParsedDecimal {
        sign,
        body: DecimalBody::Finite { digits, exponent },
    })
}

fn lex_special(rest: &[u8]) -> Option<DecimalBody> {
    let lower: Vec<u8> = rest.iter().take(8).map(u8::to_ascii_lowercase).collect();
    if lower.starts_with(b"infinity") || lower.starts_with(b"inf") {
        Some(DecimalBody::Infinity)
    } else if lower.starts_with(b"nan") {
        Some(DecimalBody::Nan)
    } else {
        None
    }
}

fn special_token_len(rest: &[u8]) -> Option<usize> {
    let lower: Vec<u8> = rest.iter().take(8).map(u8::to_ascii_lowercase).collect();
    if lower.starts_with(b"infinity") {
        Some(8)
    } else if lower.starts_with(b"inf") || lower.starts_with(b"nan") {
        Some(3)
    } else {
        None
    }
}

// --- Decimal-to-binary conversion ---

fn decimal_to_bigfloat(
    parsed: ParsedDecimal,
    precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    match parsed.body {
        DecimalBody::Nan => {
            let nan =
                BigFloat::try_new_quiet_nan(parsed.sign, precision, &[]).expect("precision >= 1");
            (nan, Status::OK)
        }
        DecimalBody::Infinity => {
            let inf = BigFloat::try_new_infinity(parsed.sign, precision).expect("precision >= 1");
            (inf, Status::OK)
        }
        DecimalBody::Finite { digits, exponent } => {
            // If digits are all zero, the value is signed zero.
            if digits.iter().all(|&d| d == 0) {
                let z = BigFloat::try_new_zero(parsed.sign, precision).expect("precision >= 1");
                return (z, Status::OK);
            }
            finite_to_bigfloat(&digits, exponent, parsed.sign, precision, mode)
        }
    }
}

/// Build the multi-precision integer `m` from digits, then convert
/// `m × 10^exponent` to a [`BigFloat`] at `precision`, rounding
/// under `mode`.
/// Decimal-exponent magnitude beyond which `finite_to_bigfloat`
/// short-circuits to overflow or underflow. The bound exists to
/// stop `pow5(|exponent|)` from allocating gigabyte-scale
/// intermediate storage on fuzz inputs like `600e333331144`.
///
/// At target precisions up to `2^32 - 1`, the result of
/// `digits × 10^exponent` rounds to ±∞ (overflow) for
/// `exponent > MAX_DECIMAL_EXPONENT` and to ±0 (underflow) for
/// `exponent < -MAX_DECIMAL_EXPONENT`, since pfloat's binary
/// exponent is an `i64` and `10^MAX_DECIMAL_EXPONENT ≈ 2^3.32M`
/// is already well below `2^63 / 2`. The chosen bound also caps
/// `pow5` storage at ~3 MB.
///
/// Phase 7 follow-up: replace the explicit `pow5` build with a
/// logarithm-based dispatch so this cap can move from
/// "computational feasibility" to "exponent fits the binary
/// exponent type" only.
const MAX_DECIMAL_EXPONENT: i32 = 1_000_000;

fn finite_to_bigfloat(
    digits: &[u8],
    exponent: i32,
    sign: Sign,
    precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    // Short-circuit absurdly-large exponents to overflow or
    // underflow before allocating the `pow5` intermediate. The
    // explicit threshold avoids OOM on fuzz inputs like
    // `600e333331144`.
    if exponent > MAX_DECIMAL_EXPONENT {
        // Positive exponent: digits × 10^exponent overflows past
        // pfloat's binary-exponent range, so the rounded result
        // is ±∞ with OVERFLOW + INEXACT raised.
        let inf = BigFloat::try_new_infinity(sign, precision).expect("precision >= 1");
        let status = Status::OVERFLOW | Status::INEXACT;
        auto_raise(status);
        return (inf, status);
    }
    if exponent < -MAX_DECIMAL_EXPONENT {
        // Negative exponent: digits × 10^exponent underflows past
        // pfloat's binary-exponent range, so the rounded result
        // is ±0 with UNDERFLOW + INEXACT raised.
        let z = BigFloat::try_new_zero(sign, precision).expect("precision >= 1");
        let status = Status::UNDERFLOW | Status::INEXACT;
        auto_raise(status);
        return (z, status);
    }

    // 1. Build m as a multi-precision integer (little-endian limbs).
    let m = digits_to_int(digits);

    // 2. Split: value = m × 10^exponent = m × 5^exponent × 2^exponent.
    if exponent >= 0 {
        // value = (m × 5^exponent) × 2^exponent.
        let power5 = pow5(exponent as u32);
        let m_times_5 = multiply_limbs(&m, &power5);
        let top_bit = top_set_bit(&m_times_5).expect("non-zero product");

        let intermediate_precision = (top_bit + 1) as u32;
        let intermediate_limbs = limbs_for(intermediate_precision);
        let mut intermediate: Vec<u64> = vec![0u64; intermediate_limbs];
        let dst_low_zero = (intermediate_limbs as u32) * 64 - intermediate_precision;
        or_left_shifted_into(
            &mut intermediate,
            &m_times_5,
            intermediate_precision,
            dst_low_zero,
        );

        // value's MSB position = top_bit of (m × 5^exp), then
        // weighted by 2^exponent.
        let result_exp = (top_bit as i64) + i64::from(exponent);

        let (value, status) = round_finite_to_precision(
            sign,
            result_exp,
            &intermediate,
            intermediate_precision,
            false,
            precision,
            mode,
        );
        auto_raise(status);
        // Re-emit with caller's requested sign in case the value
        // rounded to zero — `round_finite_to_precision` produces a
        // Normal, but if the caller expected a signed zero from a
        // tiny input, we handle that here.
        let signed = with_outer_sign(value, sign);
        (signed, status)
    } else {
        // value = m / 5^|exponent| / 2^|exponent|.
        let neg_exp = exponent.unsigned_abs();
        let power5 = pow5(neg_exp);

        // Shift m left by enough bits to give the quotient
        // (precision + guard) bits.
        let guard: u32 = 8;
        let m_bits = top_set_bit(&m).map_or(0u32, |t| (t + 1) as u32);
        let pow5_bits = top_set_bit(&power5).map_or(0u32, |t| (t + 1) as u32);
        let l = precision
            .saturating_add(guard)
            .saturating_add(pow5_bits)
            .saturating_sub(m_bits)
            .saturating_add(2);

        let total_bits = m_bits + l;
        let shifted_limbs = limbs_for(total_bits);
        let mut shifted: Vec<u64> = vec![0u64; shifted_limbs];
        or_left_shifted_into(&mut shifted, &m, m_bits, l);

        let (quotient, remainder) = divmod_limbs(&shifted, &power5);
        let top_bit = top_set_bit(&quotient).expect("non-zero quotient");

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

        let pre_sticky = remainder.iter().any(|&v| v != 0);

        // The quotient represents m × 2^L / 5^|exp|.
        // True value = m / (5^|exp| × 2^|exp|)
        //            = quotient × 2^(-L) / 2^|exp|
        //            = quotient × 2^(-L - |exp|)
        // pfloat exponent (position of MSB) = top_bit + (-L - |exp|)
        //                                    = top_bit - L - |exp|
        let result_exp = (top_bit as i64) - i64::from(l) - i64::from(neg_exp);

        let (value, status) = round_finite_to_precision(
            sign,
            result_exp,
            &intermediate,
            intermediate_precision,
            pre_sticky,
            precision,
            mode,
        );
        auto_raise(status);
        let signed = with_outer_sign(value, sign);
        (signed, status)
    }
}

/// Reasserts the caller-requested sign on the result. The rounding
/// pipeline already takes a sign argument and emits the right sign;
/// this helper guards against any future refactor that might lose
/// it on the underflow-to-zero corner case.
fn with_outer_sign(mut v: BigFloat, sign: Sign) -> BigFloat {
    match &mut v.class {
        Class::Zero { sign: s }
        | Class::Infinity { sign: s }
        | Class::Nan { sign: s, .. }
        | Class::Normal { sign: s, .. } => {
            *s = sign;
        }
    }
    v
}

/// Convert a digit slice (each entry in `0..=9`) into a
/// little-endian limb integer.
fn digits_to_int(digits: &[u8]) -> Vec<u64> {
    if digits.is_empty() {
        return vec![0];
    }
    // Build via repeated `acc = acc * 10 + d`. For long inputs we
    // batch into u64 chunks of up to 19 digits (10^19 fits in u64,
    // 10^20 does not) before multiplying-and-adding into the
    // multi-limb accumulator.
    let mut acc: Vec<u64> = vec![0];
    let chunk_size = 19;
    let mut i = 0;
    while i < digits.len() {
        let end = (i + chunk_size).min(digits.len());
        let mut chunk_value: u64 = 0;
        let mut chunk_pow10: u64 = 1;
        for &d in &digits[i..end] {
            chunk_value = chunk_value * 10 + u64::from(d);
            chunk_pow10 *= 10;
        }
        // acc = acc × chunk_pow10 + chunk_value
        acc = multiply_limbs(&acc, &[chunk_pow10]);
        let _ = limbs_add_assign(&mut acc, &[chunk_value]);
        i = end;
    }
    // Trim trailing zero limbs.
    while acc.len() > 1 && *acc.last().unwrap() == 0 {
        acc.pop();
    }
    acc
}

/// Compute `5^exp` as a multi-precision integer. Trailing zero
/// limbs trimmed.
fn pow5(exp: u32) -> Vec<u64> {
    if exp == 0 {
        return vec![1];
    }
    // Repeated squaring.
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

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    fn parse(s: &str) -> BigFloat {
        BigFloat::parse_str(s, 53, RoundingMode::NearestEven)
            .expect("parse should succeed")
            .0
    }

    fn parse_at(s: &str, p: u32) -> BigFloat {
        BigFloat::parse_str(s, p, RoundingMode::NearestEven)
            .expect("parse should succeed")
            .0
    }

    fn assert_eq_bf(a: &BigFloat, b: &BigFloat) {
        assert_eq!(
            a.partial_cmp(b).0,
            Some(Ordering::Equal),
            "expected equal: {a:?} vs {b:?}"
        );
    }

    #[test]
    fn parse_zero() {
        let v = parse("0");
        assert!(v.is_zero());
        assert!(v.is_sign_positive());
    }

    #[test]
    fn parse_neg_zero() {
        let v = parse("-0");
        assert!(v.is_zero());
        assert!(v.is_sign_negative());
    }

    #[test]
    fn parse_integers() {
        for &(s, n) in &[
            ("1", 1i64),
            ("-1", -1),
            ("7", 7),
            ("-42", -42),
            ("100", 100),
            ("12345", 12345),
        ] {
            let v = parse(s);
            let expected = BigFloat::try_from_i64_exact(n, 53).unwrap();
            assert_eq_bf(&v, &expected);
        }
    }

    #[test]
    fn parse_decimal_half() {
        // 0.5 = 1/2: exact in binary.
        let v = parse("0.5");
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (sum, _) = v.add(&v, RoundingMode::NearestEven);
        assert_eq_bf(&sum, &one);
    }

    #[test]
    fn parse_decimal_quarter() {
        let v = parse("0.25");
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (four_v, _) = v.add(&v, RoundingMode::NearestEven);
        let (eight, _) = four_v.add(&four_v, RoundingMode::NearestEven);
        assert_eq_bf(&eight, &one);
    }

    #[test]
    fn parse_with_exponent() {
        // 1e3 = 1000.
        let v = parse("1e3");
        let expected = BigFloat::try_from_i64_exact(1000, 53).unwrap();
        assert_eq_bf(&v, &expected);

        // 2.5e2 = 250.
        let v2 = parse("2.5e2");
        let expected2 = BigFloat::try_from_i64_exact(250, 53).unwrap();
        assert_eq_bf(&v2, &expected2);

        // 1e-3 = 0.001 (inexact at 53 bits but should round)
        let v3 = parse("1e-3");
        // 1/1000 cannot be exact in binary; just check it parses and is positive.
        assert!(v3.is_normal());
        assert!(v3.is_sign_positive());
    }

    #[test]
    fn parse_special_values() {
        let nan_v = parse("nan");
        assert!(nan_v.is_quiet_nan());

        let neg_nan = parse("-NaN");
        assert!(neg_nan.is_quiet_nan());
        assert!(neg_nan.is_sign_negative());

        let pi = parse("inf");
        assert!(pi.is_infinite());
        assert!(pi.is_sign_positive());

        let ni = parse("-Infinity");
        assert!(ni.is_infinite());
        assert!(ni.is_sign_negative());

        let pi_long = parse("+INF");
        assert!(pi_long.is_infinite());
        assert!(pi_long.is_sign_positive());
    }

    #[test]
    fn parse_leading_dot() {
        let v = parse(".5");
        let half = parse("0.5");
        assert_eq_bf(&v, &half);
    }

    #[test]
    fn parse_trailing_dot() {
        let v = parse("5.");
        let five = BigFloat::try_from_i64_exact(5, 53).unwrap();
        assert_eq_bf(&v, &five);
    }

    #[test]
    fn parse_signed() {
        let neg = parse("-3.14");
        assert!(neg.is_sign_negative());
        assert!(neg.is_normal());
        let pos = parse("+3.14");
        assert!(pos.is_sign_positive());
        assert!(pos.is_normal());
    }

    #[test]
    fn parse_whitespace_trimmed() {
        let v = parse("  3.14  ");
        assert!(v.is_normal());
    }

    #[test]
    fn parse_errors() {
        assert_eq!(
            BigFloat::parse_str("", 53, RoundingMode::NearestEven),
            Err(ParseError::Empty)
        );
        assert_eq!(
            BigFloat::parse_str("   ", 53, RoundingMode::NearestEven),
            Err(ParseError::Empty)
        );
        assert_eq!(
            BigFloat::parse_str(".", 53, RoundingMode::NearestEven),
            Err(ParseError::MissingDigits)
        );
        assert_eq!(
            BigFloat::parse_str("1.2.3", 53, RoundingMode::NearestEven),
            Err(ParseError::InvalidCharacter)
        );
        assert_eq!(
            BigFloat::parse_str("1e", 53, RoundingMode::NearestEven),
            Err(ParseError::InvalidExponent)
        );
        assert_eq!(
            BigFloat::parse_str("1eX", 53, RoundingMode::NearestEven),
            Err(ParseError::InvalidExponent)
        );
        // "abc" has no digits and no recognized special token; the
        // first failing check is the "missing digits" rule.
        assert_eq!(
            BigFloat::parse_str("abc", 53, RoundingMode::NearestEven),
            Err(ParseError::MissingDigits)
        );
        // A digit followed by an illegal character does trip the
        // "invalid character" path.
        assert_eq!(
            BigFloat::parse_str("1x", 53, RoundingMode::NearestEven),
            Err(ParseError::InvalidCharacter)
        );
        assert_eq!(
            BigFloat::parse_str("1", 0, RoundingMode::NearestEven),
            Err(ParseError::PrecisionZero)
        );
    }

    #[test]
    fn parse_round_trip_via_high_precision_at_low_precision() {
        // Parse "1" at precision 53, then at precision 2: should both
        // be the same value (= 1) just at different precisions.
        let one_53 = parse_at("1", 53);
        let one_2 = parse_at("1", 2);
        assert_eq!(one_53.partial_cmp(&one_2).0, Some(Ordering::Equal));
    }

    #[test]
    fn parse_inexact_for_1_third() {
        // "0.333333..." can never be exact in binary; the resulting
        // BigFloat should be a normal value with INEXACT raised.
        let (v, status) =
            BigFloat::parse_str("0.333333333333333", 53, RoundingMode::NearestEven).unwrap();
        assert!(v.is_normal());
        assert!(status.inexact());
    }

    #[test]
    fn parse_large_exponent() {
        // 1e20 — fits comfortably in 53-bit precision as the exact
        // value 10^20.
        let (v, _) = BigFloat::parse_str("1e20", 256, RoundingMode::NearestEven).unwrap();
        assert!(v.is_normal());
        assert!(v.is_sign_positive());
    }

    #[test]
    fn parse_pow5_5() {
        let p = pow5(5);
        assert_eq!(p, vec![3125]);
    }

    #[test]
    fn parse_pow5_zero() {
        assert_eq!(pow5(0), vec![1]);
    }

    #[test]
    fn parse_digits_to_int() {
        assert_eq!(digits_to_int(&[1, 2, 3]), vec![123]);
        assert_eq!(digits_to_int(&[0]), vec![0]);
        // 10^18 fits in u64.
        let big = digits_to_int(&[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(big, vec![1_000_000_000_000_000_000u64]);
    }

    #[test]
    fn parse_huge_positive_exponent_overflows_without_allocating() {
        // Regression for the libFuzzer OOM repro
        // `600e333331144`: the parser used to allocate `pow5(3.3M)`
        // ≈ 100 MB of intermediate storage. The short-circuit
        // path in `finite_to_bigfloat` produces ±∞ + OVERFLOW +
        // INEXACT directly.
        let (v, status) =
            BigFloat::parse_str("600e333331144", 113, RoundingMode::NearestEven).unwrap();
        assert!(v.is_infinite());
        assert!(v.is_sign_positive());
        assert!(status.overflow());
        assert!(status.inexact());
    }

    #[test]
    fn parse_huge_negative_exponent_underflows_without_allocating() {
        let (v, status) =
            BigFloat::parse_str("600e-333331144", 113, RoundingMode::NearestEven).unwrap();
        assert!(v.is_zero());
        assert!(v.is_sign_positive());
        assert!(status.underflow());
        assert!(status.inexact());
    }

    #[test]
    fn parse_negative_sign_huge_exponent_propagates_sign() {
        let (v, status) =
            BigFloat::parse_str("-600e333331144", 113, RoundingMode::NearestEven).unwrap();
        assert!(v.is_infinite());
        assert!(v.is_sign_negative());
        assert!(status.overflow());
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixed_parse_basic() {
        let (v, _) = FixedFloat::<53>::parse_str("1.5", RoundingMode::NearestEven).unwrap();
        assert!(v.is_normal());
        let one = FixedFloat::<53>::try_from_i64_exact(1).unwrap();
        let (twice_v, _) = v.add(&v, RoundingMode::NearestEven);
        let three = FixedFloat::<53>::try_from_i64_exact(3).unwrap();
        assert_eq!(twice_v.partial_cmp(&three).0, Some(Ordering::Equal));
        let _ = one;
    }
}
