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
//! - `Ci(±0) = −∞` + `DIV_BY_ZERO` (a pole for BOTH zero signs:
//!   `γ + ln x → −∞`, and `log(±0) = −∞` groups `−0` with the pole per
//!   IEEE 754-2019 §9.2 / C11 F.10.3.7; pf-k8ax, ADR-0123).
//! - `x < 0` (and `−∞`) ⇒ `NaN` + `INVALID` (complex in the
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
    // The asymptotic path is cheap but has an irreducible truncation
    // floor (`ci_asymptotic`); near a Ci zero the result falls below it
    // and the asymptotic certifies a wrong value. Use the asymptotic
    // only where it can correctly round, otherwise the convergent series
    // (which cancellation_boosted drives to any depth). ci_series carries
    // Ci's real zero at x ≈ 0.6165 for the small-x path and every
    // large-x near-zero the asymptotic hands off (pf-1vzg, ADR-0125).
    let use_asymptotic = e_x >= asymptotic_threshold_exponent(target_precision)
        && ci_asymptotic_reliable(x, target_precision);

    let (result, status) = ziv_round(
        |w| {
            if use_asymptotic {
                // Reliable (above the floor), but still boost: at large x
                // the near-zero cancellation can exceed the Ziv guard cap,
                // and cancellation_boosted resolves it (the asymptotic is
                // above its floor here, so it converges to the true value,
                // not the floor).
                super::ziv::cancellation_boosted(w, |ww| {
                    let (v, op_scale, _floor) = ci_asymptotic(x, ww);
                    (v, op_scale)
                })
            } else {
                super::ziv::cancellation_boosted(w, |ww| ci_series(x, ww))
            }
        },
        target_precision,
        mode,
        CI_ERROR_GUARD,
    );
    // Ci(x) at nonzero algebraic arguments is γ-entangled
    // (Ci(x) = γ + ln x + an E-function value), so its irrationality
    // there is NOT a theorem — the INEXACT force is conditionally
    // sound, the ADR-0066 γ/ζ(5) posture (ADR-0064's recorded scope
    // split; ADR-0105). The same conditionality covers Ci's real zero
    // at x ≈ 0.6165: not proven irrational, so "no dyadic input
    // yields an exact zero" is conditional on that constant ∉ ℚ₂. Ci(+0) = −∞ (pole), the
    // exact limit Ci(+∞) = +0, and x < 0 (INVALID) are dispatched
    // above and are not Class::Normal, so they keep their status.
    let status = super::force_transcendental_inexact(&result, status);
    auto_raise(status);
    (result, status)
}

/// Convergent alternating series for `Ci(x)`, `x > 0` (DLMF 6.6.6).
/// Working precision boosted by `≈ x·log₂ e` for the alternating
/// cancellation (the [`super::erf`] guard idiom).
///
/// Returns `(value, op_scale)` for [`super::ziv::cancellation_boosted`],
/// where `op_scale` is the largest partial-term exponent. Near the
/// small-`x` root `x ≈ 0.6165` the terms are `O(1)`, but as the
/// convergent (non-asymptotic) representation this series is also the
/// large-`x` near-zero fallback, where the terms peak at `≈ 2^{x·log₂ e}`
/// and the true operand scale is what charges the deep cancellation
/// (pf-1vzg, ADR-0125). The prior hardcoded `op_scale = 4` undercharged
/// the large-`x` regime.
fn ci_series(x: &BigFloat, target_precision: u32) -> (BigFloat, i64) {
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

    // Largest partial term seen: the operand scale for the alternating
    // cancellation (charged by cancellation_boosted at the call site).
    let mut max_term_exp =
        super::ziv::value_exponent(&gamma).max(super::ziv::value_exponent(&ln_x));

    // U_k = (−1)ᵏ x^{2k}/((2k)·(2k)!), k ≥ 1. U_1 = −x²/4.
    // (2k)! = (2k)(2k−1)(2k−2)! gives the ratio
    // U_k/U_{k-1} = −x²·(2k−2) / [(2k)²·(2k−1)] (verified at
    // k=2: −x²/24).
    let four = BigFloat::try_from_i64_exact(4, working_prec).expect("precision >= 1");
    let (mut term, _) = x_sq.div(&four, RoundingMode::NearestEven);
    term = term.negated();
    max_term_exp = max_term_exp.max(super::ziv::value_exponent(&term));
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
        max_term_exp = max_term_exp.max(super::ziv::value_exponent(&term));
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
    (acc, max_term_exp)
}

/// Whether the asymptotic path can correctly round `Ci(x)` at
/// `target_precision`, or whether `x` sits near a Ci zero below the
/// asymptotic truncation floor and must hand off to [`ci_series`]. The
/// shared [`super::ziv::asymptotic_reliable`] driver grows the working
/// precision until the result resolves, then applies the soundness test
/// (pf-1vzg, ADR-0125).
fn ci_asymptotic_reliable(x: &BigFloat, target_precision: u32) -> bool {
    super::ziv::asymptotic_reliable(target_precision, CI_ERROR_GUARD, |w| ci_asymptotic(x, w))
}

/// `Ci(x) = f(x)·sin x − g(x)·cos x` (DLMF 6.12.4) for large `x > 0`,
/// reusing the shared `Si`/`Ci` auxiliaries.
///
/// Returns `(value, op_scale, floor_exp)`. `f·sin` and `g·cos` cancel
/// near each large-`x` Ci zero, so `op_scale` (their larger exponent) is
/// the realised cancellation scale. But `f`/`g` are DIVERGENT asymptotic
/// series with an irreducible truncation floor `≈ e^{−x}`: below it the
/// value cannot be computed at any working precision (pf-1vzg,
/// ADR-0125). `floor_exp` carries that floor to the result scale so the
/// caller can detect a near-zero result that has fallen below it and
/// route to the convergent [`ci_series`] instead. `cancellation_boosted`
/// on the asymptotic would only converge to the wrong floor value.
fn ci_asymptotic(x: &BigFloat, target_precision: u32) -> (BigFloat, i64, i64) {
    let working_prec = target_precision
        .saturating_add(64)
        .min(target_precision.saturating_add(512));
    let x_w = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    let (f, floor_f) = si_ci_f(&x_w, working_prec);
    let (g, floor_g) = si_ci_g(&x_w, working_prec);
    let (cos_x, _) = x_w.cos(RoundingMode::NearestEven);
    let (sin_x, _) = x_w.sin(RoundingMode::NearestEven);

    let (f_sin, _) = f.mul(&sin_x, RoundingMode::NearestEven);
    let (g_cos, _) = g.mul(&cos_x, RoundingMode::NearestEven);
    let (result, _) = f_sin.sub(&g_cos, RoundingMode::NearestEven);
    let op_scale = super::ziv::value_exponent(&f_sin).max(super::ziv::value_exponent(&g_cos));
    // Floor of f·sin − g·cos: each series' truncation floor scaled by its
    // trig factor (exponents ≤ 0). The larger (shallower) dominates.
    let floor_exp = floor_f
        .saturating_add(super::ziv::value_exponent(&sin_x))
        .max(floor_g.saturating_add(super::ziv::value_exponent(&cos_x)));
    (result, op_scale, floor_exp)
}

fn acc_exponent(v: &BigFloat) -> i64 {
    match &v.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    }
}
