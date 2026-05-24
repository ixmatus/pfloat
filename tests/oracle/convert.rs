//! Conversion bridges between `f32`, `BigFloat`, and `rug::Float`.
//!
//! The oracle harness's verifier needs three conversions:
//!
//! 1. `f32` bit pattern → `BigFloat` at `p = 24` (the kernel input).
//!    For binary32 normals this matches the bit pattern directly;
//!    for binary32 subnormals the `p = 24` `BigFloat` carries more
//!    precision than `f32` itself does at that exponent range, so
//!    the construction goes through integer mantissa plus
//!    power-of-two scaling rather than a decimal round trip. This
//!    mirrors the `bf53_of_bits` pattern slice p1.2 introduced for
//!    the L-M differential lane.
//!
//! 2. `BigFloat` (any precision) → `f32` bit pattern under `NE`.
//!    The result of a kernel call to `*_round(24, NE)` is at
//!    `p = 24`; converting to `f32` bits goes through the value's
//!    decimal at the f64 round-trip digit count (17) and Rust's
//!    standard `f32` parser (same shape as the L-M lane's
//!    `bf_to_f64_bits` from slice p1.2 at `p = 53`).
//!
//!    Slice p1.4 raised the digit count from the default `Display`
//!    width (`round_trip_digit_count(24) = 9`) to 17 (= f64
//!    round-trip). The original 9-digit width carries enough
//!    decimal precision for every `f32` *normal* result, but for
//!    `f32` *subnormal* results sitting on the subnormal-grid
//!    midpoint, the 9-digit rounding flipped the decimal across
//!    the midpoint and made the `f32` parser pick the wrong
//!    neighbor (slice p1.3 sweep findings on `erf`). 17 digits
//!    captures the exact 24-bit value (`24 < 53`), so the `f32`
//!    parser sees the `BigFloat`'s actual midpoint position and
//!    applies IEEE round-to-nearest tie-to-even correctly.
//!
//! 3. `rug::Float` → `Option<f32>` under any IEEE rounding mode.
//!    The oracle's enclosure endpoints land in `f32` here; the
//!    `Option` is `None` when the value is NaN (NaN has no
//!    `f32`-rounded answer for the certified-rounding check). MPFR
//!    has no roundTiesToAway primitive, so `NearestAway` is
//!    synthesized in the same way the L-M lane's
//!    `round_ties_to_away` (from `tests/differential/mod.rs`)
//!    handles it.

#![cfg(all(unix, feature = "differential-mpfr"))]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Sign};
use rug::float::Round;
use rug::Float;

const F32_PREC: u32 = 24;
const F32_BIAS: i64 = 127;
const F32_MANT_BITS: i64 = 23;
const F32_MIN_NORMAL_EXP: i64 = -126;
const F32_MAX_NORMAL_EXP: i64 = 127;
const F32_SUBNORMAL_MIN_EXP: i64 = F32_MIN_NORMAL_EXP - F32_MANT_BITS;

/// Build a `BigFloat` at `p = 24` representing the binary32 value
/// with the given bit pattern exactly. Mirrors `bf53_of_bits` from
/// the L-M differential driver, sized for binary32.
pub fn bf24_of_bits(bits: u32) -> BigFloat {
    let sign_neg = (bits >> 31) & 1 == 1;
    let exp_field = ((bits >> 23) & 0xFF) as i64;
    let mant_field = (bits & 0x007F_FFFF) as i64;
    let sign = if sign_neg {
        Sign::Negative
    } else {
        Sign::Positive
    };

    if exp_field == 0 && mant_field == 0 {
        return BigFloat::try_new_zero(sign, F32_PREC).expect("p=24 valid");
    }
    if exp_field == 0xFF {
        if mant_field == 0 {
            return BigFloat::try_new_infinity(sign, F32_PREC).expect("p=24 valid");
        }
        return BigFloat::try_new_quiet_nan(sign, F32_PREC, &[]).expect("p=24 valid");
    }

    let (mantissa_int, scale) = if exp_field == 0 {
        // Subnormal: value = mant_field * 2^-149.
        (mant_field, F32_SUBNORMAL_MIN_EXP)
    } else {
        // Normal: value = (2^23 | mant_field) * 2^(exp_field - 127 - 23).
        (
            (1i64 << F32_MANT_BITS) | mant_field,
            exp_field - F32_BIAS - F32_MANT_BITS,
        )
    };
    let signed_mantissa = if sign_neg {
        -mantissa_int
    } else {
        mantissa_int
    };

    let mut bf = BigFloat::try_from_i64_exact(signed_mantissa, F32_PREC)
        .expect("24-bit mantissa fits in i64");

    // Scale by 2^|scale| via chained mul/div by exact powers of 2.
    // Each operation is exact in BigFloat arithmetic (just a binary
    // exponent shift), so the chain preserves the mantissa.
    let chunk_bits: u32 = 24;
    let chunk = BigFloat::try_from_i64_exact(1i64 << chunk_bits, F32_PREC).expect("2^24 fits");
    let abs = scale.unsigned_abs();
    let q = abs / u64::from(chunk_bits);
    let r = (abs % u64::from(chunk_bits)) as u32;
    for _ in 0..q {
        bf = if scale >= 0 {
            bf.mul(&chunk, RoundingMode::NearestEven).0
        } else {
            bf.div(&chunk, RoundingMode::NearestEven).0
        };
    }
    if r > 0 {
        let rem = BigFloat::try_from_i64_exact(1i64 << r, F32_PREC).expect("small power fits");
        bf = if scale >= 0 {
            bf.mul(&rem, RoundingMode::NearestEven).0
        } else {
            bf.div(&rem, RoundingMode::NearestEven).0
        };
    }
    bf
}

/// Convert a `BigFloat` (any precision) to a binary32 bit pattern,
/// rounding under `NearestEven`. Goes through the value's decimal
/// at the f64 round-trip digit count (17 digits for any precision
/// up to 53) and Rust's standard `f32` parser, which rounds the
/// decimal to nearest f32 under IEEE NE (including subnormals).
///
/// The 17-digit width (= `round_trip_digit_count(53)`) captures the
/// exact value of any `BigFloat` with precision ≤ 53. The harness's
/// kernel target is `p = 24`, well below 53, so 17 digits carries
/// the `BigFloat`'s full mantissa information and the `f32` parser
/// sees the true midpoint position on f32-subnormal-grid ties
/// (slice p1.4, closes pf-z0f). For `BigFloat`s above `p = 53` the
/// digit count scales with `bf.precision` so the round-trip stays
/// exact.
pub fn bf_to_f32_bits(bf: &BigFloat) -> u32 {
    let effective_precision = bf.precision().max(53);
    let digits = BigFloat::round_trip_digit_count(effective_precision);
    let s = bf.to_decimal_string(digits, RoundingMode::NearestEven);
    s.parse::<f32>()
        .expect("BigFloat decimal (incl. nan / inf / 0 tokens) parses as f32 under NE")
        .to_bits()
}

/// Round a `rug::Float` to `f32` under the requested rounding mode.
/// Returns `None` when the value is NaN: NaN has no unique f32 it
/// "rounds to" (`NaN != NaN` under IEEE), so the certified-rounding
/// check treats it as inconclusive at the conversion step. Infinity
/// and subnormal values do round to a determined `f32` and return
/// `Some`.
pub fn round_f32(value: &Float, mode: RoundingMode) -> Option<f32> {
    if value.is_nan() {
        return None;
    }
    match mode {
        RoundingMode::NearestEven => Some(value.to_f32_round(Round::Nearest)),
        RoundingMode::TowardZero => Some(value.to_f32_round(Round::Zero)),
        RoundingMode::TowardPositive => Some(value.to_f32_round(Round::Up)),
        RoundingMode::TowardNegative => Some(value.to_f32_round(Round::Down)),
        RoundingMode::NearestAway => Some(round_ties_to_away_f32(value)),
    }
}

/// IEEE 754 roundTiesToAway of a `rug::Float` to `f32`, synthesized
/// because MPFR offers no such mode (`MPFR_RNDA` is directed
/// round-away-from-zero, not ties-to-away). Mirrors the L-M
/// differential lane's `round_ties_to_away` helper, sized for f32.
fn round_ties_to_away_f32(value: &Float) -> f32 {
    let lo = value.to_f32_round(Round::Zero);
    let hi = value.to_f32_round(Round::AwayZero);
    if lo == hi {
        return lo;
    }
    // Distances at the value's working precision, computed without
    // further rounding.
    let g = value.prec();
    let d_lo = Float::with_val(g, value - &Float::with_val(g, lo)).abs();
    let d_hi = Float::with_val(g, &Float::with_val(g, hi) - value).abs();
    match d_lo.partial_cmp(&d_hi) {
        Some(Ordering::Less) => lo,
        Some(Ordering::Greater) => hi,
        // Exact tie: away from zero wins.
        Some(Ordering::Equal) | None => hi,
    }
}
