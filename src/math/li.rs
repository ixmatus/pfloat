//! `li(x)`: the logarithmic integral, `li(x) = Ei(ln x)` for real
//! `x > 0`, `x ≠ 1` (DLMF 6.2.8). pfloat is real-only, so the domain
//! is `x > 0`.
//!
//! The kernel is a composition: it forms `t = ln x` at a boosted
//! working precision and feeds it to the [`super::ei`] kernel, then
//! rounds once to the caller's target. No new series is introduced.
//!
//! Special cases:
//!
//! - `li(0) = 0` (the defining integral over an empty interval; not a
//!   pole, so no flag is raised even though `ln 0 = −∞`).
//! - `li(1) = −∞`, raising `DIV_BY_ZERO` (a pole: `ln 1 = +0`,
//!   `Ei(0) = −∞`).
//! - `li(+∞) = +∞`.
//! - `x < 0` ⇒ `NaN` + `INVALID` (`ln x` is undefined in the reals).
//! - `li(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::ziv::ziv_round;
use super::ziv_calibration::LI_ERROR_GUARD;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `li(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn li(&self, mode: RoundingMode) -> (Self, Status) {
        self.li_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `li(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.30, ADR-0038).
    pub fn li_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(li_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `li(self)` for `FixedFloat`. Delegates to [`BigFloat::li`].
    #[must_use]
    pub fn li(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().li(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn li_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    match &x.class {
        Class::Nan {
            quiet,
            sign,
            payload,
        } => {
            if !*quiet {
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            let nan = BigFloat::try_new_quiet_nan(*sign, target_precision, payload)
                .expect("precision >= 1");
            return (nan, Status::OK);
        }
        Class::Zero { .. } => {
            // li(0) = 0 by definition (not a pole; ln 0 = −∞,
            // Ei(−∞) = 0). No flag.
            let z =
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1");
            return (z, Status::OK);
        }
        Class::Infinity { sign } => {
            if matches!(sign, Sign::Negative) {
                // x = −∞ is outside the real domain.
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            // li(+∞) = Ei(+∞) = +∞.
            let pinf = BigFloat::try_new_infinity(Sign::Positive, target_precision)
                .expect("precision >= 1");
            return (pinf, Status::OK);
        }
        Class::Normal { sign, .. } => {
            if matches!(sign, Sign::Negative) {
                // x < 0: ln x undefined in the reals.
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
        }
    }

    // x = 1 is the li(1) = -∞ + DIV_BY_ZERO pole (ln 1 = +0,
    // Ei(0) = -∞). Handle it before Ziv so the DIV_BY_ZERO flag
    // reaches the caller's status return value, not just the
    // thread-local (the ziv_round signature returns only a rounding
    // status). This special case avoids needing a status-merging
    // multi-arg Ziv driver.
    let one_at_input = BigFloat::try_from_i64_exact(1, x.precision).expect("precision >= 1");
    if matches!(
        x.partial_cmp(&one_at_input).0,
        Some(core::cmp::Ordering::Equal)
    ) {
        let ninf =
            BigFloat::try_new_infinity(Sign::Negative, target_precision).expect("precision >= 1");
        auto_raise(Status::DIV_BY_ZERO);
        return (ninf, Status::DIV_BY_ZERO);
    }

    // x > 0 finite, x ≠ 1. Ziv-driven composition: eval(w) computes
    // t = ln(x) at w then Ei(t) at w. Both ln (Ziv-driven, slice
    // p1.2) and Ei (Ziv-driven, this slice) deliver correctly-rounded
    // NE values at working precision; the outer Ziv envelope drives
    // the rounding mode at target precision.
    let (result, status) = ziv_round(
        |w| {
            // li(x) = Ei(ln x); near the zero (the Ramanujan-Soldner
            // constant) Ei's series cancels to a near-zero value, so
            // boost the precision by the realised cancellation so the
            // Ziv half-width stays sound (review 2026-05-29, root cause
            // 2). Ei's operands are O(1) at its only positive zero, so
            // the operand scale is a small constant.
            super::ziv::cancellation_boosted(w, |ww| {
                let x_w = x
                    .round_to_precision(ww, RoundingMode::NearestEven)
                    .expect("precision >= 1")
                    .0;
                let (t, _) = x_w.ln(RoundingMode::NearestEven);
                (t.ei(RoundingMode::NearestEven).0, 4)
            })
        },
        target_precision,
        mode,
        LI_ERROR_GUARD,
    );
    auto_raise(status);
    (result, status)
}
