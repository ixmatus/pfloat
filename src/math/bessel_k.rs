//! Modified Bessel functions of the second kind `K0`, `K1`, `Kn`
//! (DLMF Chapter 10): integer order, real argument. Like
//! [`super::bessel_y`], `K` is real-valued only for `x > 0`: it has a
//! logarithmic branch point at the origin and is complex for
//! `x < 0`. The domain convention follows the [`super::ci`] /
//! [`super::li`] / [`super::bessel_y`] precedent (real-only, complex
//! off the positive axis):
//!
//! - `Kₙ(+0) = +∞`, raising `DIV_BY_ZERO` (a pole: DLMF 10.30.2
//!   `Kν(z) ∼ ½Γ(ν)(½z)⁻ν` for `ν > 0` and DLMF 10.30.3
//!   `K₀(z) ∼ −ln z` both diverge to **+∞** as `x → 0⁺`). Note the
//!   sign is the opposite of `Yₙ(+0) = −∞`.
//! - `x < 0` (and `−0`, `−∞`) ⇒ `NaN` + `INVALID` (`K` is complex in
//!   the reals there).
//! - `Kₙ(+∞) = +0` for every order, `Status::OK`. This is a
//!   **genuine exponential-decay limit** (DLMF 10.40.2
//!   `√(π/2x)·e⁻ˣ → 0`), **not** the decaying-envelope *convention*
//!   used by `J`/`Y`/Airy. There the function oscillates with a
//!   shrinking but non-converging envelope and `+0` is a
//!   conservative choice that keeps the function total; here `K`
//!   actually converges to `0`, so `+0` is the true mathematical
//!   limit. ADR-0025 records the distinction.
//! - `Kₙ(NaN) = NaN`; `sNaN` raises `INVALID`.
//!
//! The order parity is **even with no sign**: `K₋ₙ(x) = Kₙ(x)`
//! (DLMF 10.27.3), in contrast to `Y₋ₙ = (−1)ⁿ Yₙ`. There is no
//! argument-parity reduction (the domain is `x > 0`; a negative
//! argument is `INVALID`, not folded to `|x|`). So the kernel
//! evaluates `K_m(x)` for `m = |n| ≥ 0` with no sign applied at all.
//!
//! `K` is the *dominant* solution of the modified-Bessel three-term
//! recurrence in order (DLMF 10.30.2 `Kν → ∞` as `ν → ∞` at fixed
//! `z`), so the kernel computes `K₀` and `K₁` directly and climbs to
//! `Kₙ` by **upward** recurrence (DLMF 10.29.1; slice 6q.6), the
//! [`super::bessel_y`] template. The recurrence differs from the
//! ordinary-Bessel one by a sign: `𝒵_{ν−1} − 𝒵_{ν+1} = (2ν/z)𝒵_ν`
//! rather than `𝒞_{ν−1} + 𝒞_{ν+1} = (2ν/z)𝒞_ν`. The base pair
//! `K₀`/`K₁` uses two regimes on the binary exponent of `x`: the
//! DLMF 10.31.1 logarithmic series below the cut (slice 6q.5), the
//! DLMF 10.40.2 asymptotic at/above it (slice 6q.7). ADR-0025
//! records the design and the DLMF provenance.

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
    /// `K₀(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn k0(&self, mode: RoundingMode) -> (Self, Status) {
        self.k0_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `K₀(self)` with explicit result precision.
    pub fn k0_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_k_kernel(0, self, target_precision, mode))
    }

    /// `K₁(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn k1(&self, mode: RoundingMode) -> (Self, Status) {
        self.k1_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `K₁(self)` with explicit result precision.
    pub fn k1_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_k_kernel(1, self, target_precision, mode))
    }

    /// `Kₙ(self)` for integer order `n` (any `i32`, including
    /// negative: `K₋ₙ = Kₙ`), rounded under `mode` to
    /// `self.precision`.
    #[must_use]
    pub fn kn(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        self.kn_round(n, self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Kₙ(self)` with explicit result precision.
    pub fn kn_round(
        &self,
        n: i32,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_k_kernel(n, self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `K₀(self)` for `FixedFloat`. Delegates to [`BigFloat::k0`].
    #[must_use]
    pub fn k0(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().k0(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `K₁(self)` for `FixedFloat`. Delegates to [`BigFloat::k1`].
    #[must_use]
    pub fn k1(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().k1(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Kₙ(self)` for `FixedFloat`. Delegates to [`BigFloat::kn`].
    #[must_use]
    pub fn kn(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().kn(n, mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

/// `Kₙ(x)` for integer order `n` and real `x`.
///
/// Special values are handled directly per the module-level domain
/// table; for a normal positive argument the order is reduced to
/// `K_m(x)`, `m = |n|`, with no sign applied (`K₋ₙ = Kₙ`; the domain
/// is `x > 0`, so there is no argument parity), then the regime
/// evaluator runs.
fn bessel_k_kernel(
    n: i32,
    x: &BigFloat,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    let m = n.unsigned_abs();

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
        Class::Zero { sign } => {
            if matches!(sign, Sign::Negative) {
                // −0: K is complex off the positive axis (the Ci/li
                // convention; −0 groups with x < 0).
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            // Kₙ(+0) = +∞ + DIV_BY_ZERO (a pole, DLMF 10.30.2/10.30.3;
            // +∞, the opposite sign of Yₙ(+0) = −∞).
            let pinf = BigFloat::try_new_infinity(Sign::Positive, target_precision)
                .expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            (pinf, Status::DIV_BY_ZERO)
        }
        Class::Infinity { sign } => {
            if matches!(sign, Sign::Negative) {
                // −∞: complex (groups with x < 0).
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            // Kₙ(+∞) = +0, the genuine exponential-decay limit
            // (DLMF 10.40.2 √(π/2x)·e⁻ˣ → 0), Status::OK. This is a
            // true limit, not the decaying-envelope convention.
            let zero =
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1");
            (zero, Status::OK)
        }
        Class::Normal { sign, .. } => {
            if matches!(sign, Sign::Negative) {
                // x < 0: Kₙ(−x) is complex in the reals.
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }

            // K₋ₙ(x) = Kₙ(x) (DLMF 10.27.3): even in order, no sign.
            // No argument parity (x > 0 only). Evaluate K_m(x).
            let value = bessel_k_eval_normal(m, x, target_precision);

            let (rounded, status) = value
                .round_to_precision(target_precision, mode)
                .expect("precision >= 1");
            auto_raise(status);
            (rounded, status)
        }
    }
}

/// `K_m(x)` for `m ≥ 0`, normal `x > 0`: the regime evaluator.
/// Returns the unrounded working-precision value; [`bessel_k_kernel`]
/// does the single final round.
///
/// Slice 6q.1 placeholder: returns a quiet NaN so the skeleton
/// compiles and the special-value path is exercised in isolation.
/// Slices 6q.5 (log series), 6q.6 (upward recurrence), and 6q.7
/// (asymptotic + regime dispatch) replace this with the real
/// evaluator; no normal-argument tests run until then.
fn bessel_k_eval_normal(_m: u32, x: &BigFloat, target_precision: u32) -> BigFloat {
    let _ = x;
    BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[]).expect("precision >= 1")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn k_positive_zero_is_pole() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, s) = z.k0(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive(), "K0(+0) = +∞");
        assert!(s.div_by_zero());
        let (r1, s1) = z.k1(RoundingMode::NearestEven);
        assert!(r1.is_infinite() && r1.is_sign_positive(), "K1(+0) = +∞");
        assert!(s1.div_by_zero());
        let (rn, sn) = z.kn(3, RoundingMode::NearestEven);
        assert!(rn.is_infinite() && rn.is_sign_positive(), "K3(+0) = +∞");
        assert!(sn.div_by_zero());
    }

    #[test]
    fn k_negative_zero_is_invalid() {
        let z = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, s) = z.k0(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan(), "K0(−0) = NaN");
        assert!(s.invalid());
    }

    #[test]
    fn k_negative_argument_is_invalid() {
        let x = BigFloat::try_from_i64_exact(-3, 53).unwrap();
        for n in [0i32, 1, 2, -2] {
            let (r, s) = x.kn(n, RoundingMode::NearestEven);
            assert!(r.is_quiet_nan(), "K_{n}(−3) = NaN (complex)");
            assert!(s.invalid());
        }
    }

    #[test]
    fn k_positive_infinity_is_zero() {
        let inf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        for n in [0i32, 1, 5, -3] {
            let (r, st) = inf.kn(n, RoundingMode::NearestEven);
            assert!(r.is_zero() && r.is_sign_positive(), "K_{n}(+∞) = +0");
            assert!(!st.invalid());
        }
    }

    #[test]
    fn k_negative_infinity_is_invalid() {
        let ninf = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, s) = ninf.k1(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan(), "K1(−∞) = NaN (complex)");
        assert!(s.invalid());
    }

    #[test]
    fn k_quiet_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.kn(2, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn k_signaling_nan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, st) = sn.k0(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(st.invalid());
    }

    #[test]
    fn precision_zero_is_rejected() {
        let x = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert!(x.k0_round(0, RoundingMode::NearestEven).is_err());
        assert!(x.kn_round(3, 0, RoundingMode::NearestEven).is_err());
    }
}
