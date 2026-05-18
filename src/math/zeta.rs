//! Riemann zeta function `ζ(s)` for real argument `s` (DLMF
//! Chapter 25). Real-valued, single argument (no order parameter,
//! unlike Bessel). The structural template is
//! [`super::gamma_stirling`] — a Bernoulli-table asymptotic
//! correction — not the Bessel recurrence kernels.
//!
//! Two evaluation regimes (filled by slices 6r.2 / 6r.3):
//!
//! - `s ≥ 1/2`: the Euler–Maclaurin summation (DLMF 25.11), the
//!   Dirichlet series accelerated by the Bernoulli `B_{2k}`
//!   correction reused from [`super::gamma_stirling`].
//! - `s < 1/2`: the functional equation DLMF 25.4.2
//!   `ζ(s) = 2·(2π)^{s−1}·sin(πs/2)·Γ(1−s)·ζ(1−s)`, reflecting into
//!   the well-conditioned `1−s > 1/2` Euler–Maclaurin region. The
//!   reflection point `1/2` is the symmetry point of the completed
//!   functional equation `ξ(s) = ξ(1−s)` (DLMF 25.4.3/25.4.4).
//!   Routes through the in-crate `π`, `pow`, `sin`, `Γ`.
//!
//! Special values are handled directly per this domain table
//! (DLMF 25.2 / 25.6, derived not recalled; ADR-0026):
//!
//! - `ζ(NaN) = NaN`; `sNaN` raises `INVALID`.
//! - `ζ(1) = +∞`, raising `DIV_BY_ZERO`: the only singularity is a
//!   simple pole at `s = 1` with residue `+1` (DLMF 25.2). The
//!   `s → 1⁺` side is `+∞`; `+∞` is the documented pole convention,
//!   the [`super::ci`] / [`super::li`] / [`super::bessel_k`]
//!   precedent.
//! - `ζ(0) = −1/2` exact (DLMF 25.6.1), for `±0`.
//! - `ζ(−2n) = +0` exact for `n ≥ 1` (DLMF 25.6.4): the trivial
//!   zeros at the negative even integers. Special-cased here so the
//!   functional-equation path's `sin(πs/2) = 0` cancellation does
//!   not have to produce an exact zero.
//! - `ζ(+∞) = 1`, `Status::OK`: a genuine limit (the Dirichlet
//!   series DLMF 25.2.1 collapses to its first term).
//! - `ζ(−∞) = NaN`, raising `INVALID`. Via the functional equation
//!   `Γ(1−s) → +∞` super-exponentially while `sin(πs/2)` oscillates
//!   in `[−1, 1]`, so `|ζ(s)|` grows without bound *and* does not
//!   converge. This is an **unbounded non-converging oscillation**,
//!   explicitly **not** the `J`/`Y`/Airy decaying-envelope
//!   convention (a *bounded* non-converging oscillation, where `+0`
//!   is a total-keeping choice). With no limit and no finite
//!   total-keeping value, the honest convention is `NaN` +
//!   `INVALID`. ADR-0026 records the distinction (the K-vs-Y `+∞`
//!   precedent, ADR-0025).

use super::lgamma::is_integer_test;
use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use core::cmp::Ordering;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `ζ(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn zeta(&self, mode: RoundingMode) -> (Self, Status) {
        self.zeta_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `ζ(self)` with explicit result precision.
    pub fn zeta_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(zeta_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `ζ(self)` for `FixedFloat`. Delegates to [`BigFloat::zeta`].
    #[must_use]
    pub fn zeta(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().zeta(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

/// Integer `v` as a `BigFloat` at precision `p` (exact for the small
/// integers this kernel forms).
fn ci(v: i64, p: u32) -> BigFloat {
    BigFloat::try_from_i64_exact(v, p).expect("precision >= 1")
}

/// `true` if `x` is a negative even integer (`−2, −4, −6, …`), the
/// trivial zeros `ζ(−2n) = 0`, `n ≥ 1`. `x` is integral and `x/2`
/// is also integral (division by two is exact, a pure exponent
/// shift, so the test is total and exact). `x = 0` is excluded by
/// the caller (the `Class::Zero` arm); odd negative integers fail
/// the `x/2` integrality test and route to the functional equation.
fn is_negative_even_integer(x: &BigFloat) -> bool {
    if !matches!(x.sign(), Sign::Negative) || !is_integer_test(x) {
        return false;
    }
    let two = ci(2, x.precision());
    let (half, _) = x.div(&two, RoundingMode::NearestEven);
    is_integer_test(&half)
}

/// `ζ(s)` for real `s`.
///
/// Special values are handled directly per the module-level domain
/// table; a finite non-special argument routes to [`zeta_finite`]
/// (the Euler–Maclaurin / functional-equation regimes, slices
/// 6r.2 / 6r.3), then the single final round is applied.
fn zeta_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            (nan, Status::OK)
        }
        Class::Zero { .. } => {
            // ζ(0) = −1/2 exact (DLMF 25.6.1), for ±0. −1 ÷ 2 is an
            // exact exponent shift, so −1/2 is exact at any precision.
            let two = ci(2, target_precision);
            let (half, _) = ci(-1, target_precision).div(&two, RoundingMode::NearestEven);
            (half, Status::OK)
        }
        Class::Infinity { sign } => {
            if matches!(sign, Sign::Negative) {
                // ζ(−∞): unbounded non-converging oscillation via the
                // functional equation (Γ(1−s) → ∞, sin oscillates).
                // Not the decaying-envelope convention; no limit.
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            // ζ(+∞) = 1, the genuine Dirichlet-series limit (DLMF
            // 25.2.1: every term past the first vanishes).
            (ci(1, target_precision), Status::OK)
        }
        Class::Normal { .. } => {
            // Pole: ζ(1) = +∞ + DIV_BY_ZERO (simple pole, residue +1,
            // DLMF 25.2; +∞ is the s → 1⁺ side, the Ci/li/K pole
            // convention).
            let one = ci(1, x.precision());
            if matches!(x.partial_cmp(&one).0, Some(Ordering::Equal)) {
                let pinf = BigFloat::try_new_infinity(Sign::Positive, target_precision)
                    .expect("precision >= 1");
                auto_raise(Status::DIV_BY_ZERO);
                return (pinf, Status::DIV_BY_ZERO);
            }

            // Trivial zeros: ζ(−2n) = +0 exact, n ≥ 1 (DLMF 25.6.4).
            if is_negative_even_integer(x) {
                let zero = BigFloat::try_new_zero(Sign::Positive, target_precision)
                    .expect("precision >= 1");
                return (zero, Status::OK);
            }

            let value = zeta_finite(x, target_precision);
            let (rounded, status) = value
                .round_to_precision(target_precision, mode)
                .expect("precision >= 1");
            auto_raise(status);
            (rounded, status)
        }
    }
}

/// `ζ(s)` for a finite, non-special real `s` (not `0`, not `1`, not
/// a negative even integer). The Euler–Maclaurin regime (`s ≥ 1/2`)
/// lands in slice 6r.2 and the functional-equation regime
/// (`s < 1/2`) in slice 6r.3; until then this is a typed
/// placeholder that returns a quiet NaN so the crate builds and the
/// special-value dispatch above is exercised in isolation.
fn zeta_finite(x: &BigFloat, target_precision: u32) -> BigFloat {
    let _ = x;
    BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[]).expect("precision >= 1")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zeta_quiet_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, s) = q.zeta(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(!s.invalid());
    }

    #[test]
    fn zeta_signaling_nan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, st) = sn.zeta(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(st.invalid());
    }

    #[test]
    fn zeta_pole_at_one() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, s) = one.zeta(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive(), "ζ(1) = +∞");
        assert!(s.div_by_zero());
    }

    #[test]
    fn zeta_at_zero_is_minus_half() {
        // ζ(0) = −1/2 exact, for both +0 and −0.
        for sign in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(sign, 113).unwrap();
            let (r, s) = z.zeta(RoundingMode::NearestEven);
            assert!(!s.invalid() && !s.div_by_zero());
            let expected = {
                let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
                BigFloat::try_from_i64_exact(-1, 113)
                    .unwrap()
                    .div(&two, RoundingMode::NearestEven)
                    .0
            };
            assert!(
                matches!(r.partial_cmp(&expected).0, Some(Ordering::Equal)),
                "ζ(0) = −1/2"
            );
        }
    }

    #[test]
    fn zeta_trivial_zeros_at_negative_even_integers() {
        for k in [-2i64, -4, -6, -10, -42] {
            let s = BigFloat::try_from_i64_exact(k, 113).unwrap();
            let (r, st) = s.zeta(RoundingMode::NearestEven);
            assert!(
                r.is_zero() && r.is_sign_positive(),
                "ζ({k}) = +0 (trivial zero)"
            );
            assert!(!st.invalid());
        }
    }

    #[test]
    fn zeta_negative_odd_integers_are_not_trivial_zeros() {
        // ζ(−1), ζ(−3), … are nonzero rationals; they must NOT be
        // special-cased to zero (they route to the functional
        // equation in 6r.3, here the stubbed finite path → NaN).
        for k in [-1i64, -3, -5] {
            let s = BigFloat::try_from_i64_exact(k, 53).unwrap();
            let (r, _) = s.zeta(RoundingMode::NearestEven);
            assert!(!r.is_zero(), "ζ({k}) is not the zero special-case");
        }
    }

    #[test]
    fn zeta_plus_infinity_is_one() {
        let inf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, st) = inf.zeta(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert!(
            matches!(r.partial_cmp(&one).0, Some(Ordering::Equal)),
            "ζ(+∞) = 1"
        );
        assert!(!st.invalid());
    }

    #[test]
    fn zeta_negative_infinity_is_invalid() {
        let ninf = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, s) = ninf.zeta(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan(), "ζ(−∞) = NaN (no limit)");
        assert!(s.invalid());
    }

    #[test]
    fn precision_zero_is_rejected() {
        let x = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert!(x.zeta_round(0, RoundingMode::NearestEven).is_err());
    }
}
