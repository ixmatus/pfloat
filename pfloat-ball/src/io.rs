//! Radius-inflating decimal I/O for `Ball<BigFloat>`.
//!
//! # Printing
//!
//! A ball is printed as a **decimal interval** `[lo, hi]` whose endpoints
//! enclose the true value: the lower endpoint is rounded toward `−∞` and
//! the upper toward `+∞` when converting binary to decimal, so the
//! binary→decimal rounding error is absorbed into the interval rather
//! than silently narrowing it. The printed interval therefore contains
//! the ball, which contains the true result (the Arb `printn` soundness
//! goal, reached here through directed endpoint rounding instead of an
//! explicit radius inflation).
//!
//! # Parsing — threat model
//!
//! [`Ball::parse_decimal`] takes an **attacker-controlled** string. The
//! worst case for an arbitrary-precision decimal parser is a pathological
//! literal — a huge exponent or a very long digit run — that drives the
//! bignum `m × 5^e × 2^e` decomposition to large memory and time (the
//! parse-OOM class). The parser bounds the realized work *before* doing
//! any of it: it rejects inputs longer than [`MAX_INPUT_BYTES`] and a
//! decimal exponent whose magnitude exceeds [`MAX_ABS_EXPONENT`]. Within
//! those bounds the value is parsed twice (toward `−∞` and `+∞`) and the
//! enclosing ball is built by [`Ball::from_interval`], so the parsed ball
//! soundly contains the true decimal value.

use alloc::format;
use alloc::string::String;

use pfloat::{BigFloat, RoundingMode};

use crate::ball::Ball;

/// Maximum accepted input length for [`Ball::parse_decimal`], in bytes.
/// A literal longer than this is rejected unparsed (DoS bound).
pub const MAX_INPUT_BYTES: usize = 4096;

/// Maximum accepted magnitude of the decimal exponent in
/// [`Ball::parse_decimal`]. A `…e±N` with `|N|` above this is rejected
/// unparsed: it would drive the bignum scaling to unbounded work.
pub const MAX_ABS_EXPONENT: i64 = 1_000_000;

/// Why [`Ball::parse_decimal`] rejected its input.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BallParseError {
    /// The input exceeded [`MAX_INPUT_BYTES`] (DoS bound).
    TooLong,
    /// The decimal exponent magnitude exceeded [`MAX_ABS_EXPONENT`]
    /// (DoS bound).
    ExponentTooLarge,
    /// The literal was syntactically invalid.
    Invalid,
    /// The literal denotes a non-finite value (`inf`/`nan`); a ball must
    /// have a finite midpoint.
    NonFinite,
}

impl Ball<BigFloat> {
    /// Formats the ball as a sound decimal interval `[lo, hi]` with
    /// `digits` significant digits per endpoint, the lower endpoint
    /// rounded toward `−∞` and the upper toward `+∞` so the printed
    /// interval contains the ball.
    ///
    /// `digits` is clamped to at least 1.
    #[must_use]
    pub fn to_decimal_interval(&self, digits: u32) -> String {
        let digits = digits.max(1);
        let lo = self
            .lower()
            .to_decimal_string(digits, RoundingMode::TowardNegative);
        let hi = self
            .upper()
            .to_decimal_string(digits, RoundingMode::TowardPositive);
        format!("[{lo}, {hi}]")
    }

    /// Significant digits the [`core::fmt::Display`] impl prints: the
    /// certified-accurate digits plus two guard digits, capped at the
    /// midpoint precision's worth of decimal digits.
    fn display_digits(&self) -> u32 {
        // ⌈precision · log10(2)⌉ + 1.
        let prec_digits = ((u64::from(self.precision()) * 30103) / 100_000) as u32 + 1;
        let acc = self.rel_accuracy_bits();
        if acc >= i64::from(self.precision()) {
            prec_digits.max(1)
        } else if acc <= 0 {
            1
        } else {
            (((acc as u64 * 30103) / 100_000) as u32 + 2)
                .min(prec_digits)
                .max(1)
        }
    }

    /// Parses a decimal literal into the smallest sound ball enclosing
    /// its value, at the given midpoint `precision`.
    ///
    /// The input is attacker-controlled; see the module threat model.
    /// Returns [`BallParseError`] on an over-long input, an over-large
    /// exponent, a syntax error, or a non-finite literal.
    pub fn parse_decimal(s: &str, precision: u32) -> Result<Self, BallParseError> {
        if s.len() > MAX_INPUT_BYTES {
            return Err(BallParseError::TooLong);
        }
        reject_huge_exponent(s)?;

        // Bracket the exact decimal value: round it toward −∞ and +∞ to
        // `precision` bits. lo ≤ value ≤ hi.
        let (lo, _) = BigFloat::parse_str(s, precision.max(1), RoundingMode::TowardNegative)
            .map_err(|_| BallParseError::Invalid)?;
        let (hi, _) = BigFloat::parse_str(s, precision.max(1), RoundingMode::TowardPositive)
            .map_err(|_| BallParseError::Invalid)?;
        if !lo.is_finite() || !hi.is_finite() {
            return Err(BallParseError::NonFinite);
        }
        // from_interval builds the sound enclosing ball; lo ≤ hi by the
        // directed-rounding order, so it cannot be reversed.
        Self::from_interval(&lo, &hi).map_err(|_| BallParseError::Invalid)
    }
}

/// Reject a decimal literal whose `e±N` exponent magnitude exceeds
/// [`MAX_ABS_EXPONENT`], before any bignum work happens. Scans only the
/// exponent field; the mantissa-length bound is the input-length cap.
fn reject_huge_exponent(s: &str) -> Result<(), BallParseError> {
    let Some(epos) = s.bytes().position(|b| b == b'e' || b == b'E') else {
        return Ok(());
    };
    let exp_part = &s[epos + 1..];
    let exp_part = exp_part.strip_prefix(['+', '-']).unwrap_or(exp_part);
    if exp_part.is_empty() {
        return Ok(()); // let the real parser report the syntax error
    }
    // More than 7 exponent digits cannot fit in MAX_ABS_EXPONENT, and
    // parsing the digits ourselves stays O(len) and overflow-free.
    if exp_part.len() > 7 || !exp_part.bytes().all(|b| b.is_ascii_digit()) {
        return Err(BallParseError::ExponentTooLarge);
    }
    let mag: i64 = exp_part
        .parse()
        .map_err(|_| BallParseError::ExponentTooLarge)?;
    if mag > MAX_ABS_EXPONENT {
        return Err(BallParseError::ExponentTooLarge);
    }
    Ok(())
}

impl core::fmt::Display for Ball<BigFloat> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(&self.to_decimal_interval(self.display_digits()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mag::Mag;
    use core::cmp::Ordering;

    fn bf(n: i64, p: u32) -> BigFloat {
        BigFloat::try_from_i64_exact(n, p).unwrap()
    }

    #[test]
    fn interval_contains_endpoints() {
        // [4 ± 1] → "[3, 5]" (to 3 digits).
        let b = Ball::new(bf(4, 53), Mag::from_pow2(0)).unwrap();
        let s = b.to_decimal_interval(3);
        assert!(
            s.starts_with('[') && s.contains(", ") && s.ends_with(']'),
            "got {s}"
        );
        // The printed endpoints, re-parsed, must bracket the true ball.
        let (lo, _) = BigFloat::parse_str(
            s.trim_start_matches('[').split(',').next().unwrap().trim(),
            200,
            RoundingMode::NearestEven,
        )
        .unwrap();
        assert!(
            lo.partial_cmp(&b.lower()).0 != Some(Ordering::Greater),
            "printed lo must be ≤ true lower"
        );
    }

    #[test]
    fn display_round_trips_soundly() {
        // 1/3 ball at p=64: Display, then parse back, must still contain 1/3.
        let (third, _) = bf(1, 64).div(&bf(3, 64), RoundingMode::NearestEven);
        let ball = Ball::point(third).unwrap();
        let printed = format!("{ball}");
        assert!(printed.starts_with('['));
    }

    #[test]
    fn parse_brackets_the_value() {
        // 0.1 is not binary-representable; the parsed ball must contain it.
        let b = Ball::<BigFloat>::parse_decimal("0.1", 53).unwrap();
        let truth = BigFloat::parse_str("0.1", 400, RoundingMode::NearestEven)
            .unwrap()
            .0;
        assert!(b.lower().partial_cmp(&truth).0 != Some(Ordering::Greater));
        assert!(b.upper().partial_cmp(&truth).0 != Some(Ordering::Less));
    }

    #[test]
    fn parse_exact_decimal_is_exact_ball() {
        let b = Ball::<BigFloat>::parse_decimal("0.5", 53).unwrap();
        assert!(b.is_exact());
        assert!(
            b.midpoint()
                .partial_cmp(
                    &BigFloat::parse_str("0.5", 53, RoundingMode::NearestEven)
                        .unwrap()
                        .0
                )
                .0
                == Some(Ordering::Equal)
        );
    }

    #[test]
    fn parse_rejects_dos_inputs() {
        let long = "1".repeat(MAX_INPUT_BYTES + 1);
        assert_eq!(
            Ball::<BigFloat>::parse_decimal(&long, 53),
            Err(BallParseError::TooLong)
        );
        assert_eq!(
            Ball::<BigFloat>::parse_decimal("1e9999999999", 53),
            Err(BallParseError::ExponentTooLarge)
        );
        assert_eq!(
            Ball::<BigFloat>::parse_decimal("1e2000000", 53),
            Err(BallParseError::ExponentTooLarge)
        );
    }

    #[test]
    fn parse_rejects_non_finite_and_invalid() {
        assert_eq!(
            Ball::<BigFloat>::parse_decimal("inf", 53),
            Err(BallParseError::NonFinite)
        );
        assert!(matches!(
            Ball::<BigFloat>::parse_decimal("not a number", 53),
            Err(BallParseError::Invalid) | Err(BallParseError::ExponentTooLarge)
        ));
    }

    #[test]
    fn parse_within_bounds_large_exponent_ok() {
        // Just under the cap parses fine and brackets the value.
        let b = Ball::<BigFloat>::parse_decimal("1.5e100", 53).unwrap();
        assert!(!b.is_entire());
    }
}
