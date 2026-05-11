//! Operand bounding helpers for Kani harnesses.
//!
//! Each helper introduces a non-deterministic input through
//! [`kani::any`] and constrains it with [`kani::assume`] so the SAT
//! problem stays tractable. The shape mirrors
//! `ferrodec::verify::addsub::operand`; the constants are adapted to
//! pfloat's binary-float surface.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

/// Number of operands in the canonical selector.
pub(super) const NUM_CANONICAL: u8 = 8;

/// Map a small selector to one of the eight canonical
/// [`BigFloat`] values at the given precision.
///
/// 0: qNaN, 1: sNaN, 2: +∞, 3: −∞, 4: +0, 5: −0, 6: +1, 7: −1.
pub(super) fn canonical_at(idx: u8, precision: u32) -> BigFloat {
    match idx {
        0 => BigFloat::try_new_quiet_nan(Sign::Positive, precision, &[]).expect("precision >= 1"),
        1 => {
            BigFloat::try_new_signaling_nan(Sign::Positive, precision, &[]).expect("precision >= 1")
        }
        2 => BigFloat::try_new_infinity(Sign::Positive, precision).expect("precision >= 1"),
        3 => BigFloat::try_new_infinity(Sign::Negative, precision).expect("precision >= 1"),
        4 => BigFloat::try_new_zero(Sign::Positive, precision).expect("precision >= 1"),
        5 => BigFloat::try_new_zero(Sign::Negative, precision).expect("precision >= 1"),
        6 => BigFloat::try_from_i64_exact(1, precision).expect("1 fits"),
        _ => BigFloat::try_from_i64_exact(-1, precision).expect("-1 fits"),
    }
}

/// Non-deterministic [`BigFloat`] drawn from the eight-constant set.
#[cfg(kani)]
pub(super) fn nondet_constant_at(precision: u32) -> BigFloat {
    let idx: u8 = kani::any();
    kani::assume(idx < NUM_CANONICAL);
    canonical_at(idx, precision)
}

/// Non-deterministic [`RoundingMode`] over all five IEEE modes.
#[cfg(kani)]
pub(super) fn nondet_rounding_mode() -> RoundingMode {
    let idx: u8 = kani::any();
    kani::assume(idx <= 4);
    match idx {
        0 => RoundingMode::NearestEven,
        1 => RoundingMode::NearestAwayFromZero,
        2 => RoundingMode::TowardZero,
        3 => RoundingMode::TowardPositive,
        _ => RoundingMode::TowardNegative,
    }
}
