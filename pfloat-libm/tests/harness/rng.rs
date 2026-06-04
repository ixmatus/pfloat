//! Deterministic input generation for the differential lanes.
//!
//! splitmix64, the same generator pfloat's `tests/differential/mod.rs`
//! uses, so the sample is reproducible across runs and machines without
//! a `rand` dependency.

#![cfg(all(unix, feature = "differential-mpfr"))]

use pfloat_libm::RoundingMode;

/// splitmix64 step: advance `state` and return the next value.
pub fn next_u64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Draw an `i64` uniformly from `[lo, hi]` (i128 math, overflow-safe).
pub fn next_i64_in(state: &mut u64, lo: i64, hi: i64) -> i64 {
    debug_assert!(lo <= hi);
    let span = (i128::from(hi) - i128::from(lo) + 1) as u64;
    let offset = next_u64(state) % span;
    (i128::from(lo) + i128::from(offset)) as i64
}

/// A finite `f64` spread across magnitude bands: a random significand in
/// `[1, 2)` times `2^e` for `e` drawn uniformly from `[-min_exp,
/// max_exp]`, with a random sign. This concentrates the differential
/// sample on representable finite values across the full exponent range
/// rather than the mostly-infinite uniform-bit-pattern draw.
pub fn next_f64_banded(state: &mut u64, min_exp: i32, max_exp: i32) -> f64 {
    // 52 random mantissa bits give a significand in [1, 2).
    let mant = next_u64(state) & 0x000F_FFFF_FFFF_FFFF;
    let significand = 1.0_f64 + (mant as f64) / (1u64 << 52) as f64;
    let span = (max_exp - min_exp + 1) as u64;
    let e = min_exp + (next_u64(state) % span) as i32;
    let sign = if next_u64(state) & 1 == 0 { 1.0 } else { -1.0 };
    sign * significand * 2.0_f64.powi(e)
}

/// Number of random inputs per `(function, mode)` cell. Deep sweep
/// (`PFLOAT_DEEP=1`) widens it for local stress runs.
pub fn sweep_size() -> u32 {
    if std::env::var("PFLOAT_DEEP").is_ok() {
        1_000_000
    } else {
        10_000
    }
}

/// All five IEEE 754-2019 rounding modes. The shell is correctly
/// rounded under every mode (the directed-pair outer Ziv loop,
/// ADR-0057), so every lane sweeps all five.
pub const ALL_MODES: &[RoundingMode] = &[
    RoundingMode::NearestEven,
    RoundingMode::NearestAway,
    RoundingMode::TowardZero,
    RoundingMode::TowardPositive,
    RoundingMode::TowardNegative,
];
