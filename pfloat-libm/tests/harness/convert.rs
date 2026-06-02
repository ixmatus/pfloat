//! Round a `rug::Float` to a hardware float under any IEEE rounding
//! mode.
//!
//! The libm harness needs only the `rug::Float -> hardware` direction:
//! the shell already produced the hardware float under test, so the
//! oracle's job is to round its enclosure endpoints to the same width
//! and compare. There is no `BigFloat` bridge (the simplifier over
//! pfloat's `tests/oracle/convert.rs`).
//!
//! MPFR has no roundTiesToAway primitive (`MPFR_RNDA` is directed
//! round-away-from-zero, not ties-to-away), so [`RoundingMode::NearestAway`]
//! is synthesized exactly as pfloat's harness does it, here for both
//! widths.

#![cfg(all(unix, feature = "differential-mpfr"))]

use core::cmp::Ordering;

use pfloat_libm::RoundingMode;
use rug::float::Round;
use rug::Float;

/// Round a [`rug::Float`] to `f32` under `mode`. `None` for NaN: NaN has
/// no unique `f32` it "rounds to" (`NaN != NaN`), so the certified
/// rounding treats the conversion step as undecided. Infinities and
/// subnormals round to a determined `f32`.
pub fn round_f32(value: &Float, mode: RoundingMode) -> Option<f32> {
    if value.is_nan() {
        return None;
    }
    Some(match mode {
        RoundingMode::NearestEven => value.to_f32_round(Round::Nearest),
        RoundingMode::TowardZero => value.to_f32_round(Round::Zero),
        RoundingMode::TowardPositive => value.to_f32_round(Round::Up),
        RoundingMode::TowardNegative => value.to_f32_round(Round::Down),
        RoundingMode::NearestAway => round_ties_to_away_f32(value),
    })
}

/// Round a [`rug::Float`] to `f64` under `mode`. The f64 analogue of
/// [`round_f32`].
pub fn round_f64(value: &Float, mode: RoundingMode) -> Option<f64> {
    if value.is_nan() {
        return None;
    }
    Some(match mode {
        RoundingMode::NearestEven => value.to_f64_round(Round::Nearest),
        RoundingMode::TowardZero => value.to_f64_round(Round::Zero),
        RoundingMode::TowardPositive => value.to_f64_round(Round::Up),
        RoundingMode::TowardNegative => value.to_f64_round(Round::Down),
        RoundingMode::NearestAway => round_ties_to_away_f64(value),
    })
}

/// IEEE 754 roundTiesToAway of a `rug::Float` to `f32`, synthesized
/// because MPFR offers no such mode. Overflow above `f32::MAX` rounds to
/// `±inf` under both nearest modes (IEEE 754 §4.3); handle that before
/// the distance comparison, which would otherwise treat `+inf` as
/// infinitely far and wrongly pick `max_finite`.
fn round_ties_to_away_f32(value: &Float) -> f32 {
    if value.is_infinite() {
        return value.to_f32_round(Round::Nearest);
    }
    let lo = value.to_f32_round(Round::Zero);
    let hi = value.to_f32_round(Round::AwayZero);
    if lo == hi {
        return lo;
    }
    if hi.is_infinite() {
        return hi;
    }
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

/// IEEE 754 roundTiesToAway of a `rug::Float` to `f64`. The f64
/// analogue of [`round_ties_to_away_f32`].
fn round_ties_to_away_f64(value: &Float) -> f64 {
    if value.is_infinite() {
        return value.to_f64_round(Round::Nearest);
    }
    let lo = value.to_f64_round(Round::Zero);
    let hi = value.to_f64_round(Round::AwayZero);
    if lo == hi {
        return lo;
    }
    if hi.is_infinite() {
        return hi;
    }
    let g = value.prec();
    let d_lo = Float::with_val(g, value - &Float::with_val(g, lo)).abs();
    let d_hi = Float::with_val(g, &Float::with_val(g, hi) - value).abs();
    match d_lo.partial_cmp(&d_hi) {
        Some(Ordering::Less) => lo,
        Some(Ordering::Greater) => hi,
        Some(Ordering::Equal) | None => hi,
    }
}
