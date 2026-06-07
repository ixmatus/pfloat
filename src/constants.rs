//! Transcendental constants at a caller-chosen precision and rounding
//! mode.
//!
//! Each function returns the named constant correctly rounded to
//! `precision` significant bits under `mode`, paired with the
//! [`Status`] flags the rounding raised. Because the constants are
//! irrational, the returned value is inexact at every finite
//! precision, so the [`Status::INEXACT`] flag is set on every call;
//! the value itself is the closest `precision`-bit float to the true
//! constant under the requested mode.
//!
//! The functions expose the same arbitrary-precision machinery the
//! transcendental kernels use internally (a 1024-bit pinned table for
//! the common precisions, the arithmetic-geometric mean and allied
//! series above it), so a caller can request a constant at any
//! precision without reaching into the `agm` module directly. Each
//! result is driven through the Ziv interval test
//! ([`crate::math::ziv`]), so the correctly-rounded claim holds under
//! all five IEEE 754-2019 rounding modes, not just round to nearest.
//!
//! Availability follows the cluster feature that already carries the
//! underlying kernel: [`ln_2`], [`ln_10`], and [`euler_gamma`] need
//! `exp-log`; [`pi`], [`pi_over_2`], [`pi_over_4`], and [`two_over_pi`]
//! need `trig`; [`two_over_sqrt_pi`] and [`ln_2pi`] need `specials`.
//!
//! Euler's number `e`, Catalan's constant, and Apéry's constant
//! `ζ(3)` are not provided as dedicated entries yet. `e` is reachable
//! today as `BigFloat::try_from_i64_exact(1, p)?.exp(mode)`.
//!
//! # Panics
//!
//! Every function requires `precision >= 1`, the same precondition the
//! crate's other precision-taking compute methods carry. Passing
//! `precision == 0` panics. Callers holding a value whose precision came
//! from untrusted input should reject zero before calling, the way
//! [`BigFloat::round_to_precision`](crate::BigFloat::round_to_precision)
//! and the constructors return [`BuildError::PrecisionZero`] rather than
//! accepting it.
//!
//! [`BuildError::PrecisionZero`]: crate::BuildError::PrecisionZero
//!
//! # Resource use
//!
//! Working storage grows as `O(precision)`: roughly `precision / 64`
//! `u64` limbs per intermediate, with several intermediates live at once
//! during the arithmetic-geometric mean and series evaluation. A very
//! large `precision` (hundreds of millions of bits) can therefore
//! exhaust memory and abort the process. There is no built-in ceiling;
//! treat `precision` derived from untrusted input as a resource budget
//! and bound it before calling.
//!
//! # Examples
//!
//! ```
//! use pfloat::{constants, RoundingMode};
//!
//! let (ln2, status) = constants::ln_2(64, RoundingMode::NearestEven);
//! assert_eq!(ln2.precision(), 64);
//! assert!(ln2.is_finite());
//! assert!(status.inexact());
//! ```

use crate::big::BigFloat;
use crate::math::ziv::ziv_round;
use crate::math::ziv_calibration::DEFAULT_ERROR_GUARD;
use crate::rounding::RoundingMode;
use crate::status::Status;

/// Correctly round a constant to `precision` bits under `mode`.
///
/// `at(working)` returns the constant rounded to `working` bits under
/// round to nearest; its error against the true value is at most one
/// unit in the last place of the working precision, the premise the
/// Ziv interval test in [`ziv_round`] assumes. The interval test then
/// certifies the rounding to `precision` under the caller's mode,
/// escalating the working precision until the certification succeeds.
fn round_constant(
    at: fn(u32) -> BigFloat,
    precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    debug_assert!(precision >= 1);
    ziv_round(at, precision, mode, DEFAULT_ERROR_GUARD)
}

/// Returns `π = 3.14159…` correctly rounded to `precision` bits under
/// `mode`.
///
/// # Examples
///
/// ```
/// use pfloat::{constants, RoundingMode};
///
/// let (pi, _) = constants::pi(53, RoundingMode::NearestEven);
/// assert_eq!(pi.precision(), 53);
/// assert!(pi.is_finite());
/// ```
#[cfg(feature = "trig")]
#[must_use]
pub fn pi(precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    round_constant(crate::math::pi_at, precision, mode)
}

/// Returns `π/2 = 1.57079…` correctly rounded to `precision` bits
/// under `mode`.
#[cfg(feature = "trig")]
#[must_use]
pub fn pi_over_2(precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    round_constant(crate::math::pi_over_2_at, precision, mode)
}

/// Returns `π/4 = 0.78539…` correctly rounded to `precision` bits
/// under `mode`.
#[cfg(feature = "trig")]
#[must_use]
pub fn pi_over_4(precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    round_constant(crate::math::pi_over_4_at, precision, mode)
}

/// Returns `2/π = 0.63661…` correctly rounded to `precision` bits
/// under `mode`.
#[cfg(feature = "trig")]
#[must_use]
pub fn two_over_pi(precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    round_constant(crate::math::two_over_pi_at, precision, mode)
}

/// Returns `ln 2 = 0.69314…` correctly rounded to `precision` bits
/// under `mode`.
#[must_use]
pub fn ln_2(precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    round_constant(crate::math::ln_2_at, precision, mode)
}

/// Returns `ln 10 = 2.30258…` correctly rounded to `precision` bits
/// under `mode`.
#[must_use]
pub fn ln_10(precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    round_constant(crate::math::ln_10_at, precision, mode)
}

/// Returns the Euler–Mascheroni constant `γ = 0.57721…` correctly
/// rounded to `precision` bits under `mode`.
#[must_use]
pub fn euler_gamma(precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    round_constant(crate::math::euler_gamma_at, precision, mode)
}

/// Returns `2/√π = 1.12837…` correctly rounded to `precision` bits
/// under `mode`. This is the leading coefficient of the Maclaurin
/// series for `erf`.
#[cfg(feature = "specials")]
#[must_use]
pub fn two_over_sqrt_pi(precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    round_constant(crate::math::two_over_sqrt_pi_at, precision, mode)
}

/// Returns `ln 2π = 1.83787…` correctly rounded to `precision` bits
/// under `mode`. This is the constant term in Stirling's asymptotic
/// series for `ln Γ(z)`.
#[cfg(feature = "specials")]
#[must_use]
pub fn ln_2pi(precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    round_constant(crate::math::ln_2pi_at, precision, mode)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    // Reference decimals carry enough significant digits (~120) that
    // `parse_str` rounds them to the correctly-rounded float for every
    // target precision exercised here (<= 256 bits ~= 78 digits). They
    // are mathematical facts from standard references, used only to
    // pin the constants' rounding; the constants themselves derive
    // from the in-crate pinned tables and AGM series.
    #[cfg(feature = "trig")]
    const PI_REFERENCE: &str = "3.14159265358979323846264338327950288419716939937510\
        582097494459230781640628620899862803482534211706798214808651328230664709";

    const LN2_REFERENCE: &str = "0.69314718055994530941723212145817656807550013436025\
        525412068000949339362196969471560586332699641868754200148102057068573368";

    const TARGETS: &[u32] = &[53, 113, 256];

    const ALL_MODES: &[RoundingMode] = &[
        RoundingMode::NearestEven,
        RoundingMode::NearestAway,
        RoundingMode::TowardZero,
        RoundingMode::TowardPositive,
        RoundingMode::TowardNegative,
    ];

    fn eq(a: &BigFloat, b: &BigFloat) -> bool {
        matches!(a.partial_cmp(b).0, Some(Ordering::Equal))
    }

    // The constant must equal the correctly-rounded parse of its
    // reference decimal under every mode and target precision. This is
    // a rug-independent correctness check (the differential lane adds
    // the MPFR cross-check) that exercises the Ziv driver's directed
    // rounding.
    #[cfg(feature = "trig")]
    #[test]
    fn pi_matches_parse_reference() {
        for &mode in ALL_MODES {
            for &target in TARGETS {
                let (got, _) = pi(target, mode);
                let want = BigFloat::parse_str(PI_REFERENCE, target, mode)
                    .expect("parse")
                    .0;
                assert!(eq(&got, &want), "pi p={target} {mode:?}: {got} vs {want}");
            }
        }
    }

    #[test]
    fn ln_2_matches_parse_reference() {
        for &mode in ALL_MODES {
            for &target in TARGETS {
                let (got, _) = ln_2(target, mode);
                let want = BigFloat::parse_str(LN2_REFERENCE, target, mode)
                    .expect("parse")
                    .0;
                assert!(eq(&got, &want), "ln_2 p={target} {mode:?}: {got} vs {want}");
            }
        }
    }

    // An irrational constant is inexact at every finite precision.
    #[test]
    fn inexact_flag_is_set() {
        let (_, status) = ln_2(53, RoundingMode::NearestEven);
        assert!(status.inexact(), "ln_2 should round inexactly at p=53");
    }

    // Directed modes bracket the true value: the round-down result is
    // at or below the round-up result, and for the positive constant π
    // round-toward-zero coincides with round-toward-negative.
    #[cfg(feature = "trig")]
    #[test]
    fn directed_modes_bracket() {
        let target = 53;
        let (down, _) = pi(target, RoundingMode::TowardNegative);
        let (up, _) = pi(target, RoundingMode::TowardPositive);
        let (zero, _) = pi(target, RoundingMode::TowardZero);
        assert!(
            matches!(
                down.partial_cmp(&up).0,
                Some(Ordering::Less | Ordering::Equal)
            ),
            "TowardNegative {down} must be <= TowardPositive {up}"
        );
        assert!(
            eq(&zero, &down),
            "TowardZero must equal TowardNegative for +pi"
        );
    }

    // π/2 and π/4 share π's mantissa, so doubling recovers the parent
    // exactly at the same precision and mode.
    #[cfg(feature = "trig")]
    #[test]
    fn pi_subdivisions_double_back() {
        let two = BigFloat::try_from_i64_exact(2, 113).expect("p>=1");
        for &mode in ALL_MODES {
            let (p1, _) = pi(113, mode);
            let (p2, _) = pi_over_2(113, mode);
            let (p4, _) = pi_over_4(113, mode);
            let p2_doubled = p2.mul(&two, RoundingMode::NearestEven).0;
            let p4_doubled = p4.mul(&two, RoundingMode::NearestEven).0;
            assert!(eq(&p2_doubled, &p1), "2*(pi/2) != pi under {mode:?}");
            assert!(eq(&p4_doubled, &p2), "2*(pi/4) != pi/2 under {mode:?}");
        }
    }
}
