//! Shared helpers for the MPFR differential test lane.
//!
//! Each `tests/differential_<op>.rs` integration-test file does
//! `mod differential;` to pull this module in. Cargo treats the
//! file as a shared submodule of every test crate that imports it
//! (not as a standalone test binary).
//!
//! ADR-0014 records the comparison strategy. Binary radix means
//! one canonical normalized form per finite value, so bit-for-bit
//! equality of `rug::Float` values is the right test.

#![cfg(all(unix, feature = "differential-mpfr"))]

use pfloat::{BigFloat, RoundingMode};
use rug::float::Round;
use rug::Float;

/// Map pfloat's [`RoundingMode`] to MPFR's [`Round`].
pub fn mpfr_round_of(mode: RoundingMode) -> Round {
    match mode {
        RoundingMode::NearestEven => Round::Nearest,
        RoundingMode::NearestAway => Round::AwayZero,
        RoundingMode::TowardZero => Round::Zero,
        RoundingMode::TowardPositive => Round::Up,
        RoundingMode::TowardNegative => Round::Down,
    }
}

/// Convert a [`BigFloat`] to a [`rug::Float`] at the same precision.
///
/// Uses pfloat's [`core::fmt::Display`] output, which renders
/// `round_trip_digit_count(p)` decimal digits — exactly enough that
/// re-parsing at precision `p` recovers the original value. MPFR's
/// `Float::parse` is correctly rounded, so the round-trip preserves
/// the bit pattern.
pub fn bigfloat_to_rug(value: &BigFloat) -> Float {
    let p = value.precision();
    let s = value.to_string();
    let parsed = Float::parse(&s).expect("BigFloat Display must produce valid input");
    Float::with_val(p, parsed)
}

/// Construct a [`BigFloat`] at the given precision from an [`i64`].
///
/// Exact for any `n` whose magnitude fits in `p` bits, which is the
/// only case slice 6a exercises.
pub fn bigfloat_from_i64(n: i64, p: u32) -> BigFloat {
    BigFloat::try_from_i64_exact(n, p).expect("i64 fits in precision")
}

/// Construct a [`rug::Float`] at the given precision from an [`i64`].
pub fn rug_from_i64(n: i64, p: u32) -> Float {
    Float::with_val(p, n)
}

/// Number of random pairs exercised per `(precision, mode)` cell in
/// each cargo-test run. The deep sweep (10⁶) runs locally under
/// `PFLOAT_DEEP=1` per ADR-0014.
pub fn sweep_size() -> u32 {
    if std::env::var("PFLOAT_DEEP").is_ok() {
        1_000_000
    } else {
        10_000
    }
}

/// The four precisions exercised by the CI sweep.
pub const SWEEP_PRECISIONS: &[u32] = &[53, 113, 256, 1024];

/// All five IEEE 754-2019 rounding modes.
pub const ALL_ROUNDING_MODES: &[RoundingMode] = &[
    RoundingMode::NearestEven,
    RoundingMode::NearestAway,
    RoundingMode::TowardZero,
    RoundingMode::TowardPositive,
    RoundingMode::TowardNegative,
];
