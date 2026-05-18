//! Bessel functions of the second kind `Y0`, `Y1`, `Yn` (DLMF
//! Chapter 10): ordinary Bessel of integer order, real argument.
//!
//! Unlike [`super::bessel_j`], `Y` is real-valued only for `x > 0`:
//! `Y_n` has a logarithmic branch point at the origin and is complex
//! for `x < 0`. The domain convention follows the [`super::ci`] /
//! [`super::li`] precedent (cosine / logarithmic integral, same
//! "real-only, complex off the positive axis" shape):
//!
//! - `Y_n(+0) = −∞`, raising `DIV_BY_ZERO` (a pole: the DLMF 10.8.1
//!   `−(½x)^{−n}/π` head, and `(2/π) ln(½x) J_0` for `n = 0`, both
//!   diverge to `−∞` as `x → 0⁺`).
//! - `x < 0` (and `−0`, `−∞`) ⇒ `NaN` + `INVALID` (`Y` is complex in
//!   the reals there).
//! - `Y_n(+∞) = +0` for every order, by the decaying-envelope
//!   convention (ADR-0021/0023, the [`super::airy`] / J precedent):
//!   the true behaviour at `+∞` is a bounded decaying oscillation
//!   with no limit; the conservative total result is `+0`,
//!   `Status::OK`.
//! - `Y_n(NaN) = NaN`; `sNaN` raises `INVALID`.
//!
//! Negative order reduces before evaluation: `Y₋ₙ(x) = (−1)ⁿ Yₙ(x)`
//! (DLMF 10.4.1, the same parity as `J`), so the kernel evaluates
//! `Y_m(x)` for `m = |n| ≥ 0` and applies one parity sign. There is
//! no argument-parity reduction (the domain is `x > 0`; a negative
//! argument is `INVALID`, not folded to `|x|`).
//!
//! `Y` is the *dominant* solution of the Bessel three-term
//! recurrence, so the kernel computes `Y₀` and `Y₁` directly and
//! climbs to `Yₙ` by **upward** recurrence
//! `Y_{k+1}(x) = (2k/x)·Y_k(x) − Y_{k−1}(x)` (DLMF 10.6.1), which is
//! stable for the dominant solution. This is the opposite of `J`'s
//! Miller backward descent; [`super::bessel_j::bessel_j_miller`] is
//! not reused (there is no recessive solution to renormalise). The
//! base pair `Y₀`/`Y₁` uses two regimes on the binary exponent of
//! `x`, sharing [`super::bessel_j::bessel_j_threshold`] with `J`:
//!
//! - Below threshold: the DLMF 10.8.1 logarithmic series, with
//!   working precision boosted `≈ x·log₂e` for the alternating
//!   cancellation (the [`super::ci`] guard idiom). `Y` has no
//!   recessive-normalisation trick, so unlike `J` there is no cheap
//!   middle "moderate" regime; the log series carries everything
//!   below the asymptotic cut.
//! - At/above threshold: the DLMF 10.17.4 Hankel asymptotic, reusing
//!   the J `a_k(ν)` coefficients (DLMF 10.17.1, pinned in ADR-0023)
//!   with `Y`'s trig combination.
//!
//! ADR-0024 records the design and the coefficient provenance.

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
    /// `Y₀(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn y0(&self, mode: RoundingMode) -> (Self, Status) {
        self.y0_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Y₀(self)` with explicit result precision.
    pub fn y0_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_y_kernel(0, self, target_precision, mode))
    }

    /// `Y₁(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn y1(&self, mode: RoundingMode) -> (Self, Status) {
        self.y1_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Y₁(self)` with explicit result precision.
    pub fn y1_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_y_kernel(1, self, target_precision, mode))
    }

    /// `Yₙ(self)` for integer order `n` (any `i32`, including
    /// negative: `Y₋ₙ = (−1)ⁿ Yₙ`), rounded under `mode` to
    /// `self.precision`.
    #[must_use]
    pub fn yn(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        self.yn_round(n, self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Yₙ(self)` with explicit result precision.
    pub fn yn_round(
        &self,
        n: i32,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_y_kernel(n, self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `Y₀(self)` for `FixedFloat`. Delegates to [`BigFloat::y0`].
    #[must_use]
    pub fn y0(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().y0(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Y₁(self)` for `FixedFloat`. Delegates to [`BigFloat::y1`].
    #[must_use]
    pub fn y1(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().y1(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Yₙ(self)` for `FixedFloat`. Delegates to [`BigFloat::yn`].
    #[must_use]
    pub fn yn(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().yn(n, mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

/// `Yₙ(x)` for integer order `n` and real `x`.
///
/// Special values are handled directly per the module-level domain
/// table; for a normal positive argument the order is reduced to
/// `Y_m(x)`, `m = |n|`, with one parity sign, then the regime
/// evaluator runs.
fn bessel_y_kernel(
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
                // −0: Y is complex off the positive axis (the Ci/li
                // convention; −0 groups with x < 0).
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            // Y_n(+0) = −∞ + DIV_BY_ZERO (a pole, DLMF 10.8.1).
            let ninf = BigFloat::try_new_infinity(Sign::Negative, target_precision)
                .expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            (ninf, Status::DIV_BY_ZERO)
        }
        Class::Infinity { sign } => {
            if matches!(sign, Sign::Negative) {
                // −∞: complex (groups with x < 0).
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }
            // Y_n(+∞) = +0 (decaying-envelope convention,
            // ADR-0021/0023).
            let zero =
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1");
            (zero, Status::OK)
        }
        Class::Normal { sign, .. } => {
            if matches!(sign, Sign::Negative) {
                // x < 0: Y_n(−x) is complex in the reals.
                let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                    .expect("precision >= 1");
                auto_raise(Status::INVALID);
                return (nan, Status::INVALID);
            }

            // Y₋ₙ(x) = (−1)ⁿ Yₙ(x) (DLMF 10.4.1): order parity only,
            // negate when m is odd and n < 0.
            let negate = (m % 2 == 1) && (n < 0);
            let value = bessel_y_eval_normal(m, x, target_precision);
            let value = if negate { value.negated() } else { value };

            let (rounded, status) = value
                .round_to_precision(target_precision, mode)
                .expect("precision >= 1");
            auto_raise(status);
            (rounded, status)
        }
    }
}

/// `Y_m(x)` for `m ≥ 0`, normal `x > 0`: the base pair plus upward
/// recurrence. Returns the unrounded working-precision value;
/// [`bessel_y_kernel`] does the single final round.
///
/// Slice 6p.1 placeholder: returns a quiet NaN so the skeleton
/// compiles and the special-value path is exercised in isolation.
/// Slices 6p.2 (log series), 6p.3 (asymptotic), and 6p.4 (upward
/// recurrence + dispatch) replace this with the real evaluator; no
/// normal-argument tests run until then.
fn bessel_y_eval_normal(_m: u32, x: &BigFloat, target_precision: u32) -> BigFloat {
    let _ = x;
    BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[]).expect("precision >= 1")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn y_positive_zero_is_pole() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, s) = z.y0(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative(), "Y0(+0) = −∞");
        assert!(s.div_by_zero());
        let (r1, s1) = z.y1(RoundingMode::NearestEven);
        assert!(r1.is_infinite() && r1.is_sign_negative(), "Y1(+0) = −∞");
        assert!(s1.div_by_zero());
        let (rn, sn) = z.yn(3, RoundingMode::NearestEven);
        assert!(rn.is_infinite() && rn.is_sign_negative(), "Y3(+0) = −∞");
        assert!(sn.div_by_zero());
    }

    #[test]
    fn y_negative_zero_is_invalid() {
        let z = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, s) = z.y0(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan(), "Y0(−0) = NaN");
        assert!(s.invalid());
    }

    #[test]
    fn y_negative_argument_is_invalid() {
        let x = BigFloat::try_from_i64_exact(-3, 53).unwrap();
        for n in [0i32, 1, 2, -2] {
            let (r, s) = x.yn(n, RoundingMode::NearestEven);
            assert!(r.is_quiet_nan(), "Y_{n}(−3) = NaN (complex)");
            assert!(s.invalid());
        }
    }

    #[test]
    fn y_positive_infinity_is_zero() {
        let inf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        for n in [0i32, 1, 5, -3] {
            let (r, st) = inf.yn(n, RoundingMode::NearestEven);
            assert!(r.is_zero() && r.is_sign_positive(), "Y_{n}(+∞) = +0");
            assert!(!st.invalid());
        }
    }

    #[test]
    fn y_negative_infinity_is_invalid() {
        let ninf = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, s) = ninf.y1(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan(), "Y1(−∞) = NaN (complex)");
        assert!(s.invalid());
    }

    #[test]
    fn y_quiet_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.yn(2, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn y_signaling_nan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, st) = sn.y0(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(st.invalid());
    }

    #[test]
    fn precision_zero_is_rejected() {
        let x = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert!(x.y0_round(0, RoundingMode::NearestEven).is_err());
        assert!(x.yn_round(3, 0, RoundingMode::NearestEven).is_err());
    }
}
