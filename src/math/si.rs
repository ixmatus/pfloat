//! `Si(x) = ∫₀ˣ (sin t)/t dt`: the sine integral (DLMF §6.2). Entire
//! and odd; defined for all real `x`.
//!
//! Two regimes, dispatched on the binary exponent of `|x|` like
//! [`super::erf`], computed on `|x|` with the sign reapplied
//! afterwards (oddness):
//!
//! - Small `|x|`: the convergent alternating series (DLMF 6.6.5)
//!   ```text
//!   Si(x) = Σ_{k≥0} (−1)ᵏ x^{2k+1} / ((2k+1) · (2k+1)!)
//!   ```
//!   with the working precision boosted by `≈ |x|·log₂ e` to absorb
//!   the alternating cancellation (the [`super::erf`] guard idiom).
//!
//! - Large `|x|`: the auxiliary-function form (DLMF 6.12.3)
//!   ```text
//!   Si(x) = π/2 − f(x)·cos x − g(x)·sin x
//!   ```
//!   with `f`, `g` the shared asymptotic auxiliaries
//!   [`si_ci_f`]/[`si_ci_g`] (DLMF 6.12.1–2), summed to their
//!   smallest term.
//!
//! Special cases:
//!
//! - `Si(±0) = ±0` (odd).
//! - `Si(+∞) = π/2`, `Si(−∞) = −π/2`.
//! - `Si(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::pi_over_2_at;
use super::ziv::ziv_round;
use super::ziv_calibration::SI_ERROR_GUARD;
use super::{pi_over_2_at_round, signed_constant_at_round};

impl BigFloat {
    /// `Si(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn si(&self, mode: RoundingMode) -> (Self, Status) {
        self.si_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Si(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.30, ADR-0038).
    pub fn si_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(si_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `Si(self)` for `FixedFloat`. Delegates to [`BigFloat::si`].
    #[must_use]
    pub fn si(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().si(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn si_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
        Class::Zero { sign } => {
            // Si(±0) = ±0 (odd).
            let z = BigFloat::try_new_zero(*sign, target_precision).expect("precision >= 1");
            return (z, Status::OK);
        }
        Class::Infinity { sign } => {
            // Si(±∞) = ±π/2 (odd). Mode-aware: the negative case rounds
            // π/2 under the mirrored mode before negating (Phase 4
            // directed-mode constant audit; Si(−∞, TowardNegative) used to
            // land above −π/2). See signed_constant_at_round.
            let (result, status) =
                signed_constant_at_round(pi_over_2_at_round, *sign, target_precision, mode);
            auto_raise(status);
            return (result, status);
        }
        Class::Normal { .. } => {}
    }

    // Tiny x: Si(x) = x − x³/18 + … shrinks toward x in magnitude (the
    // −x³/18 correction opposes x's sign), so round x with that
    // magnitude-shrinking infinitesimal directly — the ADR-0059/0104
    // tiny-x pattern (Si is the asinh analogue, asinh(x) = x − x³/6).
    // Past the Ziv guard cap (target + 1024) the series correction
    // (position ~2·|e_x| below x) is unreachable, the driver collapses
    // to x and rounds it as if exact, and directed modes returned x
    // itself — 1 ulp wrong toward zero where Si(x) < x (pf-31ql,
    // ADR-0113). The two-part depth clears both the target ulp and the
    // INPUT's grid (pf-fbjn, ADR-0104); arm-failing inputs go to the
    // driver's deep rung.
    let e = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => unreachable!("special classes dispatched above"),
    };
    if e <= -(i64::from(target_precision) + 2)
        && e.saturating_mul(-2) >= i64::from(x.precision).saturating_add(6)
    {
        return crate::rounding::round_with_infinitesimal(
            x,
            x.sign(),
            true, // magnitude shrinks: the −x³/18 correction opposes x's sign
            target_precision,
            mode,
        );
    }

    // Regime decision pinned from target_precision so it does not
    // flip across Ziv retries. Si is odd; compute on |x| and reapply
    // sign inside the eval closure so the Ziv interval test sees the
    // signed value.
    let sign = x.sign();
    let abs_x = x.abs();
    let e_x = match &abs_x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    };
    let use_asymptotic = e_x >= asymptotic_threshold_exponent(target_precision);

    let (result, status) = ziv_round(
        |w| {
            let result_abs = if use_asymptotic {
                si_asymptotic(&abs_x, w)
            } else {
                si_series(&abs_x, w)
            };
            if matches!(sign, Sign::Negative) {
                result_abs.negated()
            } else {
                result_abs
            }
        },
        target_precision,
        mode,
        SI_ERROR_GUARD,
    );
    // Si(x) for finite-normal x ≠ 0 is transcendental (the sine
    // integral takes transcendental values at nonzero algebraic
    // arguments), so the rounded result is INEXACT even where the
    // working-precision evaluation collapses onto a grid value — e.g.
    // Si(2⁻ᵏ) → 2⁻ᵏ when the −x³/18 term falls below the working
    // precision, yet the true value differs (pf-njs5 under-report,
    // ADR-0064). Si(±0) = ±0 and the exact limits Si(±∞) = ±π/2 are
    // dispatched above; only finite-normal x ≠ 0 reaches here.
    let status = super::force_transcendental_inexact(&result, status);
    auto_raise(status);
    (result, status)
}

/// Smallest binary exponent of `|x|` at which the asymptotic
/// auxiliaries already give `target_precision + 32` bits. The `f`/`g`
/// series' smallest term is `≈ e^{−|x|}`, so the requirement is
/// `2^{e_x} ≥ (p+32)·ln 2`, the same form (and the same
/// conservative rational `ln 2 ≈ 6932/10000`) as
/// [`super::ei::asymptotic_threshold_exponent`].
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

/// Convergent alternating series for `Si(|x|)` (DLMF 6.6.5). Working
/// precision boosted by `≈ |x|·log₂ e` for the alternating
/// cancellation (the [`super::erf::erf_maclaurin`] idiom).
fn si_series(abs_x: &BigFloat, target_precision: u32) -> BigFloat {
    let e_x = match &abs_x.class {
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

    let x_w = abs_x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let (x_sq, _) = x_w.mul(&x_w, RoundingMode::NearestEven);

    // T_k = (−1)ᵏ x^{2k+1}/((2k+1)·(2k+1)!), T_0 = x. Using
    // (2k+1)! = (2k+1)(2k)(2k−1)!, the ratio is
    // T_k/T_{k-1} = −x²·(2k−1) / [(2k+1)²·(2k)] (verified at
    // k=1: −x²/18). Carried below as ·(−x²)·(2k−1) /(2k) /(2k+1)
    // /(2k+1).
    let mut term = x_w.clone();
    let mut acc = x_w.clone();
    let floor = -i64::from(working_prec) - 8;
    let max_iter: i64 = 1 << 22;
    for k in 1..=max_iter {
        // Move T_{k-1} (= (-1)^{k-1} x^{2k-1}/((2k-1)(2k-1)!)) to
        // T_k by multiplying by −x² and the factorial/denominator
        // adjustment (2k-1)/( (2k)(2k+1)(2k+1) ).
        let two_k = 2 * k;
        let d1 = BigFloat::try_from_i64_exact(2 * k - 1, working_prec).expect("precision >= 1");
        let d2 = BigFloat::try_from_i64_exact(two_k, working_prec).expect("precision >= 1");
        let d3 = BigFloat::try_from_i64_exact(2 * k + 1, working_prec).expect("precision >= 1");
        let d4 = BigFloat::try_from_i64_exact(2 * k + 1, working_prec).expect("precision >= 1");
        let (t1, _) = term.mul(&x_sq, RoundingMode::NearestEven);
        let (t2, _) = t1.mul(&d1, RoundingMode::NearestEven);
        let (t3, _) = t2.div(&d2, RoundingMode::NearestEven);
        let (t4, _) = t3.div(&d3, RoundingMode::NearestEven);
        let (t5, _) = t4.div(&d4, RoundingMode::NearestEven);
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

/// `Si(|x|) = π/2 − f·cos|x| − g·sin|x|` (DLMF 6.12.3) for large
/// `|x|`.
fn si_asymptotic(abs_x: &BigFloat, target_precision: u32) -> BigFloat {
    let working_prec = target_precision
        .saturating_add(64)
        .min(target_precision.saturating_add(512));
    let x_w = abs_x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    // Si(x) → π/2 with no real zero for x > 0, so its asymptotic never
    // sits below the truncation floor; the floor return is for Ci only.
    let (f, _) = si_ci_f(&x_w, working_prec);
    let (g, _) = si_ci_g(&x_w, working_prec);
    let (cos_x, _) = x_w.cos(RoundingMode::NearestEven);
    let (sin_x, _) = x_w.sin(RoundingMode::NearestEven);
    let half_pi = pi_over_2_at(working_prec);

    let (f_cos, _) = f.mul(&cos_x, RoundingMode::NearestEven);
    let (g_sin, _) = g.mul(&sin_x, RoundingMode::NearestEven);
    let (a, _) = half_pi.sub(&f_cos, RoundingMode::NearestEven);
    let (result, _) = a.sub(&g_sin, RoundingMode::NearestEven);
    result
}

/// Shared asymptotic auxiliary `f(x) ∼ (1/x)·Σ_{k≥0} (−1)ᵏ (2k)!/x^{2k}`
/// (DLMF 6.12.1), summed to its smallest term. Used by both `Si` and
/// `Ci`. `x > 0` is assumed (callers pass `|x|`).
pub(super) fn si_ci_f(x: &BigFloat, working_prec: u32) -> (BigFloat, i64) {
    let (x_sq, _) = x.mul(x, RoundingMode::NearestEven);
    let mut term = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let mut sum = term.clone();
    let mut prev_mag = magnitude(&term);
    let max_iter: i64 = 1 << 22;
    for k in 1..=max_iter {
        // c_k/c_{k-1} = −(2k)(2k−1)/x².
        let a = BigFloat::try_from_i64_exact(2 * k, working_prec).expect("precision >= 1");
        let b = BigFloat::try_from_i64_exact(2 * k - 1, working_prec).expect("precision >= 1");
        let (t1, _) = term.mul(&a, RoundingMode::NearestEven);
        let (t2, _) = t1.mul(&b, RoundingMode::NearestEven);
        let (t3, _) = t2.div(&x_sq, RoundingMode::NearestEven);
        let cand = t3.negated();
        let mag = magnitude(&cand);
        if mag > prev_mag {
            break; // smallest term passed: optimal truncation.
        }
        prev_mag = mag;
        term = cand;
        let (sum_next, _) = sum.add(&term, RoundingMode::NearestEven);
        sum = sum_next;
        if negligible(&term, &sum, working_prec) {
            break;
        }
    }
    let (f, _) = sum.div(x, RoundingMode::NearestEven);
    // The asymptotic series is divergent; its irreducible truncation
    // floor is the smallest retained term `prev_mag` (at sum scale),
    // carried to the value scale by the final `/x`. Ci uses this to
    // detect that a near-zero result has fallen below what the
    // asymptotic can compute (pf-1vzg, ADR-0125). Si ignores it.
    let floor_exp = prev_mag.saturating_sub(magnitude(x));
    (f, floor_exp)
}

/// Shared asymptotic auxiliary
/// `g(x) ∼ (1/x²)·Σ_{k≥0} (−1)ᵏ (2k+1)!/x^{2k}` (DLMF 6.12.2),
/// summed to its smallest term. Used by both `Si` and `Ci`.
pub(super) fn si_ci_g(x: &BigFloat, working_prec: u32) -> (BigFloat, i64) {
    let (x_sq, _) = x.mul(x, RoundingMode::NearestEven);
    let mut term = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let mut sum = term.clone();
    let mut prev_mag = magnitude(&term);
    let max_iter: i64 = 1 << 22;
    for k in 1..=max_iter {
        // c_k/c_{k-1} = −(2k+1)(2k)/x².
        let a = BigFloat::try_from_i64_exact(2 * k + 1, working_prec).expect("precision >= 1");
        let b = BigFloat::try_from_i64_exact(2 * k, working_prec).expect("precision >= 1");
        let (t1, _) = term.mul(&a, RoundingMode::NearestEven);
        let (t2, _) = t1.mul(&b, RoundingMode::NearestEven);
        let (t3, _) = t2.div(&x_sq, RoundingMode::NearestEven);
        let cand = t3.negated();
        let mag = magnitude(&cand);
        if mag > prev_mag {
            break;
        }
        prev_mag = mag;
        term = cand;
        let (sum_next, _) = sum.add(&term, RoundingMode::NearestEven);
        sum = sum_next;
        if negligible(&term, &sum, working_prec) {
            break;
        }
    }
    let (g1, _) = sum.div(&x_sq, RoundingMode::NearestEven);
    // Truncation floor of the divergent series at the value scale: the
    // smallest retained term `prev_mag` carried through the final `/x²`.
    let floor_exp = prev_mag.saturating_sub(2i64.saturating_mul(magnitude(x)));
    (g1, floor_exp)
}

fn acc_exponent(v: &BigFloat) -> i64 {
    match &v.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => 0,
    }
}

fn magnitude(v: &BigFloat) -> i64 {
    match &v.class {
        Class::Normal { exponent, .. } => *exponent,
        Class::Zero { .. } => i64::MIN,
        _ => i64::MAX,
    }
}

fn negligible(term: &BigFloat, sum: &BigFloat, working_prec: u32) -> bool {
    match &term.class {
        Class::Zero { .. } => true,
        Class::Normal { exponent, .. } => {
            *exponent < acc_exponent(sum) - i64::from(working_prec) - 8
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use crate::big::BigFloat;
    use crate::rounding::RoundingMode;
    use crate::sign::Sign;
    use core::cmp::Ordering;

    #[test]
    fn si_inf_directed_rounding_is_sound() {
        // Regression (Phase 4 directed-mode constant audit): Si(−∞) = −π/2
        // (Si is odd) used to round on the wrong side of −π/2.
        let hp = crate::math::pi_over_2_at(600);
        let nhp = hp.negated();
        for &p in &[24u32, 53, 113, 200] {
            let pi = BigFloat::try_new_infinity(Sign::Positive, p).unwrap();
            let ni = BigFloat::try_new_infinity(Sign::Negative, p).unwrap();
            assert_ne!(
                pi.si(RoundingMode::TowardNegative).0.partial_cmp(&hp).0,
                Some(Ordering::Greater),
                "Si(+inf, TN) ≤ π/2 at p={p}"
            );
            assert_ne!(
                pi.si(RoundingMode::TowardPositive).0.partial_cmp(&hp).0,
                Some(Ordering::Less),
                "Si(+inf, TP) ≥ π/2 at p={p}"
            );
            assert_ne!(
                ni.si(RoundingMode::TowardNegative).0.partial_cmp(&nhp).0,
                Some(Ordering::Greater),
                "Si(-inf, TN) ≤ −π/2 at p={p}"
            );
            assert_ne!(
                ni.si(RoundingMode::TowardPositive).0.partial_cmp(&nhp).0,
                Some(Ordering::Less),
                "Si(-inf, TP) ≥ −π/2 at p={p}"
            );
        }
    }
}
