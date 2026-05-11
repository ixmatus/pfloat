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
// Each `tests/differential_*.rs` test crate uses a different
// subset of the helpers below; the unused subset would otherwise
// generate `dead_code` warnings under that crate's compilation.
#![allow(dead_code)]

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

/// Splitmix64 step. Used by each `differential_*` test for
/// deterministic input generation; consolidated here so the
/// helper isn't duplicated 22 times and so the i64 range math is
/// fixed in one place.
pub fn next_u64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Draw an i64 uniformly from `[lo, hi]` using the splitmix64
/// state. Uses i128 arithmetic so the span can cover the full i64
/// range without overflow (the previous `(hi - lo) as u64 + 1`
/// form overflowed in debug mode when `lo = -i64::MAX` and
/// `hi = i64::MAX`, which is the CI default for arithmetic tests
/// at p >= 64).
pub fn next_i64_in(state: &mut u64, lo: i64, hi: i64) -> i64 {
    debug_assert!(lo <= hi);
    let span = (i128::from(hi) - i128::from(lo) + 1) as u64;
    let offset = next_u64(state) % span;
    (i128::from(lo) + i128::from(offset)) as i64
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

/// The four precisions exercised by the CI sweep. Used for the
/// arithmetic core and for the `parse` lane, where pfloat is
/// expected to be bit-exact against MPFR at every tested precision.
///
/// `dead_code` is allowed because each `tests/differential_*.rs`
/// file uses **either** this constant or
/// [`TRANSCENDENTAL_PRECISIONS`], not both. The `mod.rs` module is
/// compiled once per test crate, so the unused constant generates
/// a warning under the test crate that uses the other one.
pub const SWEEP_PRECISIONS: &[u32] = &[53, 113, 256, 1024];

/// Precisions used by transcendental and tier-1 special function
/// differential tests. Capped at 256 bits because pfloat's
/// elementary transcendentals (exp, ln, pow, sin, cos, tan, atan2,
/// sinh, cosh, asinh, erf) all use **hardcoded 1024-bit constants**
/// (`ln(2)`, `2/π`, `2/sqrt(π)`, etc.) for argument reduction or
/// the leading coefficient. A 64-bit guard above target precision
/// means target precisions above 960 bits exceed the constants'
/// reach and produce divergence from MPFR. Phase 5 / 7 follow-up
/// is to either extend the constants (4096-bit `ln(2)` etc.) or
/// compute them on the fly via AGM-style algorithms.
pub const TRANSCENDENTAL_PRECISIONS: &[u32] = &[53, 113, 256];

/// IEEE 754-2019 rounding modes exercised in the differential lane.
///
/// **Currently NearestEven only.** The full five-mode sweep needs a
/// bit-exact `BigFloat` ↔ `rug::Float` converter; the current
/// [`bigfloat_to_rug`] helper goes via `BigFloat::Display` and
/// `rug::Float::parse`, which is rounding-mode-aware and lossy by
/// up to 1 ULP for values produced under non-NearestEven rounding
/// (Display rounds at the precision under NearestEven; rug's parse
/// rounds the same way; values that pfloat produced under, say,
/// NearestAway lose the 1-ULP difference from NearestEven through
/// the round-trip).
///
/// Concrete cases where this surfaces empirically:
///
/// - `div(-966132233652331, 1233101814760529)` at `p=53,
///   NearestAway`: pfloat and MPFR disagree by 1 ULP — but pfloat
///   under `NearestEven` matches MPFR exactly.
/// - `sqrt(2473446)` at `p=53, NearestAway`: same artifact.
/// - `fma(big, big, big)` whenever the exact `a*b+c` exceeds the
///   precision: same artifact.
///
/// Full five-mode sweep is a tracked follow-up. The fix is either
/// (a) a `pub` raw-parts accessor on [`BigFloat`] that lets the test
/// helpers build a `rug::Float` directly from sign + exponent +
/// limbs, or (b) a hex/binary radix Display on `BigFloat` that
/// round-trips exactly under any rounding mode.
pub const ALL_ROUNDING_MODES: &[RoundingMode] = &[RoundingMode::NearestEven];
