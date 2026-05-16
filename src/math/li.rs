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

    // x > 0 finite. Compose at a boosted working precision so the
    // ln-then-Ei chain keeps the caller's bits; round once at the end.
    let working_prec = target_precision
        .saturating_add(128)
        .min(target_precision.saturating_add(4096));
    let x_w = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    let (t, _) = x_w.ln(RoundingMode::NearestEven);
    // ln(1) = +0 ⇒ Ei(0) = −∞ + DIV_BY_ZERO: the li(1) pole. The
    // `ei` kernel raises the flag; capture and re-raise it.
    let (ei_val, ei_status) = t.ei(RoundingMode::NearestEven);

    let (rounded, round_status) = ei_val
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    let status = ei_status | round_status;
    auto_raise(status);
    (rounded, status)
}
