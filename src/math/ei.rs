//! `Ei(x)`: the exponential integral, the Cauchy principal value of
//! `∫_{-∞}^{x} eᵗ/t dt` for real `x ≠ 0` (DLMF §6.2).
//!
//! Two regimes, dispatched on the binary exponent of `|x|` exactly
//! like [`super::erf`]:
//!
//! - Small `|x|`: the convergent series (DLMF 6.6.2)
//!   ```text
//!   Ei(x) = γ + ln|x| + Σ_{k≥1} xᵏ / (k · k!)
//!   ```
//!   The working precision is boosted to absorb the alternating
//!   cancellation that dominates for `x < 0` (the terms alternate in
//!   sign while `Ei(x) → 0⁻`).
//!
//! - Large `|x|`: the divergent asymptotic series (DLMF 6.12.2)
//!   ```text
//!   Ei(x) ∼ (eˣ / x) · Σ_{k≥0} k! / xᵏ
//!   ```
//!   summed to its smallest term (optimal truncation near `k ≈ |x|`).
//!   The same closed form covers `x → +∞` and `x → −∞`
//!   (`Ei(x) = −E₁(−x)` reduces to it).
//!
//! The threshold is the smallest `|x|` at which the asymptotic's
//! smallest term already delivers `target + 32` bits, mirroring
//! [`super::erf::asymptotic_threshold_exponent`].
//!
//! Special cases:
//!
//! - `Ei(±0) = −∞`, raising `DIV_BY_ZERO` (a pole, the `ln(±0)`
//!   convention; DLMF 6.2.5 has `Ei(x) → −∞` as `x → 0`).
//! - `Ei(+∞) = +∞`; `Ei(−∞) = −0` (`Ei(x) → 0⁻`).
//! - `Ei(NaN) = NaN`; `sNaN` raises `INVALID`.

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
use super::ziv::ziv_round;

impl BigFloat {
    /// `Ei(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn ei(&self, mode: RoundingMode) -> (Self, Status) {
        self.ei_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Ei(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.30, ADR-0038).
    pub fn ei_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(ei_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `Ei(self)` for `FixedFloat`. Delegates to [`BigFloat::ei`].
    #[must_use]
    pub fn ei(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().ei(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn ei_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            // Ei(±0) = −∞ + DIV_BY_ZERO (a pole, the ln(±0) rule).
            let ninf = BigFloat::try_new_infinity(Sign::Negative, target_precision)
                .expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            return (ninf, Status::DIV_BY_ZERO);
        }
        Class::Infinity { sign } => {
            if matches!(sign, Sign::Negative) {
                // Ei(−∞) = 0⁻.
                let z = BigFloat::try_new_zero(Sign::Negative, target_precision)
                    .expect("precision >= 1");
                return (z, Status::OK);
            }
            let pinf = BigFloat::try_new_infinity(Sign::Positive, target_precision)
                .expect("precision >= 1");
            return (pinf, Status::OK);
        }
        Class::Normal { .. } => {}
    }

    // Regime decision is pinned from target_precision so it does not
    // flip across Ziv retries (slice p1.4 erf precedent). The
    // working-precision boost inside ei_series (the +extra term for
    // |x|·log₂ e cancellation) is applied INSIDE the eval closure.
    let e_x = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    let use_asymptotic = e_x >= asymptotic_threshold_exponent(target_precision);

    let (result, status) = ziv_round(
        |w| {
            if use_asymptotic {
                ei_asymptotic(x, w)
            } else {
                ei_series(x, w)
            }
        },
        target_precision,
        mode,
    );
    auto_raise(status);
    (result, status)
}

/// Smallest binary exponent of `|x|` at which the asymptotic series
/// already gives `target_precision + 32` bits.
///
/// Derivation: the asymptotic's smallest term is `≈ e^{−|x|}`, so a
/// relative error below `2^{−(p+32)}` needs `|x| ≥ (p+32)·ln 2`.
/// With `|x| ≥ 2^{e_x}` the requirement is `2^{e_x} ≥ (p+32)·ln 2`.
/// `ln 2 ≈ 0.6932` is taken slightly high (rational `6932/10000`) so
/// the threshold is never placed below where the asymptotic is
/// genuinely accurate; the integer search keeps the helper
/// no_std-clean.
pub(super) fn asymptotic_threshold_exponent(target_precision: u32) -> i64 {
    let bits_needed: u64 = u64::from(target_precision) + 32;
    let need: u64 = (bits_needed * 6932).div_ceil(10000);
    let mut e: i64 = 1;
    let mut pow_2: u64 = 2;
    while pow_2 < need && e < 60 {
        e += 1;
        pow_2 = pow_2.saturating_mul(2);
    }
    e
}

/// Convergent series `Ei(x) = γ + ln|x| + Σ_{k≥1} xᵏ/(k·k!)`
/// (DLMF 6.6.2). The peak term sits near `k ≈ |x|`; for `x < 0` the
/// terms alternate while the result decays, so the working precision
/// is boosted by roughly `2·|x|·log₂ e` to absorb the cancellation
/// (the [`super::erf::erf_maclaurin`] guard idiom, doubled for the
/// alternating case).
fn ei_series(x: &BigFloat, target_precision: u32) -> BigFloat {
    let e_x = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    let extra = if e_x <= 0 {
        64
    } else {
        let shift = (e_x + 1).min(20) as u32;
        let mag: u64 = 1u64 << shift;
        (mag.saturating_mul(23) / 16).saturating_mul(2).min(4096) as u32
    };
    let working_prec = target_precision
        .saturating_add(64)
        .saturating_add(extra)
        .min(target_precision.saturating_add(4096));

    let x_w = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let abs_x = x_w.abs();
    let (ln_abs, _) = abs_x.ln(RoundingMode::NearestEven);
    let gamma = euler_gamma_at(working_prec);
    let (g_plus_ln, _) = gamma.add(&ln_abs, RoundingMode::NearestEven);

    // T_k = xᵏ/(k·k!). T_1 = x; T_k = T_{k-1}·x·(k−1)/k² (since
    // k! = k·(k−1)!, so the per-term 1/k is already carried). The
    // sum is Σ_{k≥1} T_k with no extra factor.
    let mut term = x_w.clone();
    let (mut acc, _) = g_plus_ln.add(&term, RoundingMode::NearestEven);

    let floor = -i64::from(working_prec) - 8;
    let max_iter: i64 = 1 << 22;
    for k in 2..=max_iter {
        let k_big = BigFloat::try_from_i64_exact(k, working_prec).expect("precision >= 1");
        let km1 = BigFloat::try_from_i64_exact(k - 1, working_prec).expect("precision >= 1");
        let (k_sq, _) = k_big.mul(&k_big, RoundingMode::NearestEven);
        let (tx, _) = term.mul(&x_w, RoundingMode::NearestEven);
        let (txm, _) = tx.mul(&km1, RoundingMode::NearestEven);
        let (t_next, _) = txm.div(&k_sq, RoundingMode::NearestEven);
        term = t_next;
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

/// Divergent asymptotic `Ei(x) ∼ (eˣ/x)·Σ_{k≥0} k!/xᵏ` (DLMF 6.12.2),
/// summed to its smallest term. The same form is valid for both
/// signs of `x`.
fn ei_asymptotic(x: &BigFloat, target_precision: u32) -> BigFloat {
    let working_prec = target_precision
        .saturating_add(64)
        .min(target_precision.saturating_add(512));
    let x_w = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    // Σ_{k≥0} k!/xᵏ : a_0 = 1, a_k = a_{k-1} · k / x. Stop at the
    // smallest term (|a_k| ≥ |a_{k-1}|) or once a term is negligible.
    let mut term = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let mut sum = term.clone();
    let mut prev_mag = term_magnitude(&term);
    let floor = -i64::from(working_prec) - 8;
    let max_iter: i64 = 1 << 22;
    for k in 1..=max_iter {
        let k_big = BigFloat::try_from_i64_exact(k, working_prec).expect("precision >= 1");
        let (tk, _) = term.mul(&k_big, RoundingMode::NearestEven);
        let (t_next, _) = tk.div(&x_w, RoundingMode::NearestEven);
        term = t_next;

        let mag = term_magnitude(&term);
        if mag > prev_mag {
            // Asymptotic series turned divergent: stop before adding.
            break;
        }
        prev_mag = mag;

        let (sum_next, _) = sum.add(&term, RoundingMode::NearestEven);
        sum = sum_next;

        let negligible = match &term.class {
            Class::Zero { .. } => true,
            Class::Normal { exponent, .. } => *exponent < acc_exponent(&sum) + floor,
            _ => false,
        };
        if negligible {
            break;
        }
    }

    let (exp_x, _) = x_w.exp(RoundingMode::NearestEven);
    let (exp_over_x, _) = exp_x.div(&x_w, RoundingMode::NearestEven);
    let (result, _) = exp_over_x.mul(&sum, RoundingMode::NearestEven);
    result
}

/// Binary exponent of a normal accumulator, `i64::MIN`-safe for the
/// zero/non-normal cases (treated as "no magnitude").
fn acc_exponent(v: &BigFloat) -> i64 {
    match &v.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    }
}

/// Coarse magnitude proxy (the binary exponent) for the
/// smallest-term test; `None`-like values map below everything.
fn term_magnitude(v: &BigFloat) -> i64 {
    match &v.class {
        Class::Normal { exponent, .. } => *exponent,
        Class::Zero { .. } => i64::MIN,
        _ => i64::MAX,
    }
}
