//! `Ci(x) = γ + ln x + ∫₀ˣ (cos t − 1)/t dt`: the cosine integral
//! (DLMF §6.2). pfloat is real-only and `Ci(−x) = Ci(x) − iπ` is
//! complex, so the domain is `x > 0`.
//!
//! Two regimes, dispatched on the binary exponent of `x` like
//! [`super::erf`]:
//!
//! - Small `x`: the convergent alternating series (DLMF 6.6.6)
//!   ```text
//!   Ci(x) = γ + ln x + Σ_{k≥1} (−1)ᵏ x^{2k} / ((2k)·(2k)!)
//!   ```
//!   with the working precision boosted by `≈ x·log₂ e` for the
//!   alternating cancellation.
//!
//! - Large `x`: the auxiliary-function form (DLMF 6.12.4)
//!   ```text
//!   Ci(x) = f(x)·sin x − g(x)·cos x
//!   ```
//!   reusing the shared [`super::si::si_ci_f`]/[`super::si::si_ci_g`]
//!   asymptotic auxiliaries and the shared threshold.
//!
//! Special cases:
//!
//! - `Ci(+0) = −∞`, raising `DIV_BY_ZERO` (a pole: `Ci(x) ∼ γ + ln x`
//!   as `x → 0⁺`).
//! - `Ci(+∞) = +0`.
//! - `x < 0` (and `−∞`, `−0`) ⇒ `NaN` + `INVALID` (complex in the
//!   reals; supersedes any "Ci is even" shorthand).
//! - `Ci(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::euler_gamma_at;
use super::si::{asymptotic_threshold_exponent, si_ci_f, si_ci_g};
use super::ziv::ziv_round;
use super::ziv_calibration::CI_ERROR_GUARD;

impl BigFloat {
    /// `Ci(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn ci(&self, mode: RoundingMode) -> (Self, Status) {
        self.ci_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Ci(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.30, ADR-0038).
    pub fn ci_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(ci_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `Ci(self)` for `FixedFloat`. Delegates to [`BigFloat::ci`].
    #[must_use]
    pub fn ci(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().ci(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn ci_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            // Ci(+0) = −∞ + DIV_BY_ZERO (a pole: γ + ln x → −∞).
            let ninf = BigFloat::try_new_infinity(Sign::Negative, target_precision)
                .expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            return (ninf, Status::DIV_BY_ZERO);
        }
        Class::Infinity { sign } => {
            if matches!(sign, Sign::Negative) {
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            // Ci(+∞) = +0.
            let z =
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1");
            return (z, Status::OK);
        }
        Class::Normal { sign, .. } => {
            if matches!(sign, Sign::Negative) {
                // x < 0: Ci(−x) = Ci(x) − iπ is complex.
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
        }
    }

    // Regime decision pinned from target_precision so it does not
    // flip across Ziv retries (slice p1.4 erf precedent). The
    // working-precision boost inside ci_series (the +extra term for
    // x·log₂ e alternating cancellation) is applied INSIDE the eval
    // closure.
    let e_x = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    let use_asymptotic = e_x >= asymptotic_threshold_exponent(target_precision);

    let (result, status) = ziv_round(
        |w| {
            if use_asymptotic {
                ci_asymptotic(x, w)
            } else {
                // Ci has a real zero (~0.61650549) where γ + ln|x| + Σ
                // is a near-total cancellation of O(1) terms; boost by
                // the realised cancellation so the Ziv half-width stays
                // sound (review 2026-05-29, root cause 2). The operands
                // are O(1) near that zero, so the operand scale is a
                // small constant.
                super::ziv::cancellation_boosted(w, |ww| (ci_series(x, ww), 4))
            }
        },
        target_precision,
        mode,
        CI_ERROR_GUARD,
    );
    // Ci(x) for finite-normal x > 0 is transcendental (the cosine
    // integral takes transcendental values at nonzero algebraic
    // arguments), so the rounded result is INEXACT even where the
    // working-precision evaluation lands on a grid value (ADR-0064).
    // Ci's real zero at x ≈ 0.6165 is itself transcendental, so no
    // dyadic input yields an exact zero. Ci(+0) = −∞ (pole), the
    // exact limit Ci(+∞) = +0, and x < 0 (INVALID) are dispatched
    // above and are not Class::Normal, so they keep their status.
    let status = if matches!(result.class, Class::Normal { .. }) {
        status | Status::INEXACT
    } else {
        status
    };
    auto_raise(status);
    (result, status)
}

/// Convergent alternating series for `Ci(x)`, `x > 0` (DLMF 6.6.6).
/// Working precision boosted by `≈ x·log₂ e` for the alternating
/// cancellation (the [`super::erf`] guard idiom).
fn ci_series(x: &BigFloat, target_precision: u32) -> BigFloat {
    let e_x = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    let extra = if e_x <= 0 {
        64
    } else {
        let shift = (e_x + 1).min(20) as u32;
        let mag: u64 = 1u64 << shift;
        (mag.saturating_mul(23) / 16).min(4096) as u32
    };
    let working_prec = target_precision
        .saturating_add(64)
        .saturating_add(extra)
        .min(target_precision.saturating_add(4096));

    let x_w = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let (x_sq, _) = x_w.mul(&x_w, RoundingMode::NearestEven);
    let (ln_x, _) = x_w.ln(RoundingMode::NearestEven);
    let gamma = euler_gamma_at(working_prec);
    let (mut acc, _) = gamma.add(&ln_x, RoundingMode::NearestEven);

    // U_k = (−1)ᵏ x^{2k}/((2k)·(2k)!), k ≥ 1. U_1 = −x²/4.
    // (2k)! = (2k)(2k−1)(2k−2)! gives the ratio
    // U_k/U_{k-1} = −x²·(2k−2) / [(2k)²·(2k−1)] (verified at
    // k=2: −x²/24).
    let four = BigFloat::try_from_i64_exact(4, working_prec).expect("precision >= 1");
    let (mut term, _) = x_sq.div(&four, RoundingMode::NearestEven);
    term = term.negated();
    let (acc1, _) = acc.add(&term, RoundingMode::NearestEven);
    acc = acc1;

    let floor = -i64::from(working_prec) - 8;
    let max_iter: i64 = 1 << 22;
    for k in 2..=max_iter {
        let num = BigFloat::try_from_i64_exact(2 * k - 2, working_prec).expect("precision >= 1");
        let d1 = BigFloat::try_from_i64_exact(2 * k, working_prec).expect("precision >= 1");
        let d2 = BigFloat::try_from_i64_exact(2 * k, working_prec).expect("precision >= 1");
        let d3 = BigFloat::try_from_i64_exact(2 * k - 1, working_prec).expect("precision >= 1");
        let (t1, _) = term.mul(&x_sq, RoundingMode::NearestEven);
        let (t2, _) = t1.mul(&num, RoundingMode::NearestEven);
        let (t3, _) = t2.div(&d1, RoundingMode::NearestEven);
        let (t4, _) = t3.div(&d2, RoundingMode::NearestEven);
        let (t5, _) = t4.div(&d3, RoundingMode::NearestEven);
        term = t5.negated();
        let (acc_next, _) = acc.add(&term, RoundingMode::NearestEven);
        acc = acc_next;

        if k > (1i64 << e_x.max(0)).saturating_add(2) {
            let small = match &term.class {
                Class::Zero { .. } => true,
                Class::Normal { exponent, .. } => *exponent < acc_exponent(&acc) + floor,
                _ => false,
            };
            if small {
                break;
            }
        }
    }
    acc
}

/// `Ci(x) = f(x)·sin x − g(x)·cos x` (DLMF 6.12.4) for large `x > 0`,
/// reusing the shared `Si`/`Ci` auxiliaries.
fn ci_asymptotic(x: &BigFloat, target_precision: u32) -> BigFloat {
    let working_prec = target_precision
        .saturating_add(64)
        .min(target_precision.saturating_add(512));
    let x_w = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    let f = si_ci_f(&x_w, working_prec);
    let g = si_ci_g(&x_w, working_prec);
    let (cos_x, _) = x_w.cos(RoundingMode::NearestEven);
    let (sin_x, _) = x_w.sin(RoundingMode::NearestEven);

    let (f_sin, _) = f.mul(&sin_x, RoundingMode::NearestEven);
    let (g_cos, _) = g.mul(&cos_x, RoundingMode::NearestEven);
    let (result, _) = f_sin.sub(&g_cos, RoundingMode::NearestEven);
    result
}

fn acc_exponent(v: &BigFloat) -> i64 {
    match &v.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    }
}
