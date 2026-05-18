//! Modified Bessel functions of the first kind `I0`, `I1`, `In`
//! (DLMF Chapter 10): integer order, real argument. Like
//! [`super::bessel_j`], `I` is entire on the real line, with no poles
//! and no domain restriction; unlike the oscillatory ordinary Bessel
//! functions it grows monotonically and diverges at infinity.
//!
//! Order and sign reduce before evaluation. The order parity is
//! **even with no sign**: `I₋ₙ(x) = Iₙ(x)` (DLMF 10.27.1), in
//! contrast to `J`/`Y` where `𝒞₋ₙ = (−1)ⁿ 𝒞ₙ`. The argument parity
//! matches `J`: `Iₙ(−x) = (−1)ⁿ Iₙ(x)` (from the `(x/2)ⁿ` prefactor
//! of the DLMF 10.25.2 series, whose remaining sum is even in `x`).
//! So the kernel evaluates `I_m(|x|)` for `m = |n| ≥ 0` and negates
//! exactly when `m` is odd and `x < 0`; the order sign never
//! contributes.
//!
//! Special cases:
//!
//! - `I₀(±0) = 1`, `Iₙ(±0) = 0` for `n ≠ 0` (exact, DLMF 10.30.1
//!   `Iν(z) ∼ (½z)ν/Γ(ν+1)`; entire, both zero signs alike).
//! - `Iₙ(+∞) = +∞`; `Iₙ(−∞) = (−1)ⁿ·∞` (so `−∞` for odd `n`). This
//!   is a **genuine infinite limit** (`I` grows like
//!   `eˣ/√(2πx) → ∞`, DLMF 10.30.4), `Status::OK`, the
//!   `exp(+∞) = +∞` precedent — explicitly **not** the
//!   decaying-envelope convention of `J`/`Y`/Airy (which covers a
//!   bounded non-converging oscillation, a different situation).
//! - `Iₙ(NaN) = NaN`; `sNaN` raises `INVALID`.
//!
//! Three regimes will dispatch on the binary exponent of `|x|` (the
//! [`super::bessel_j`] template, since `I` is the *recessive*
//! solution in order just as `J` is): tiny all-positive Maclaurin
//! (DLMF 10.25.2, slice 6q.2), Miller backward recurrence normalised
//! by the DLMF 10.35.5 sum rule `eˣ = I₀ + 2Σ_{k≥1} Iₖ` (slice 6q.3),
//! and the DLMF 10.40.1 asymptotic reusing the ADR-0023 `aₖ(ν)`
//! coefficients (slice 6q.4). ADR-0025 records the design and the
//! DLMF provenance.

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
    /// `I₀(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn i0(&self, mode: RoundingMode) -> (Self, Status) {
        self.i0_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `I₀(self)` with explicit result precision.
    pub fn i0_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_i_kernel(0, self, target_precision, mode))
    }

    /// `I₁(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn i1(&self, mode: RoundingMode) -> (Self, Status) {
        self.i1_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `I₁(self)` with explicit result precision.
    pub fn i1_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_i_kernel(1, self, target_precision, mode))
    }

    /// `Iₙ(self)` for integer order `n` (any `i32`, including
    /// negative: `I₋ₙ = Iₙ`), rounded under `mode` to
    /// `self.precision`.
    #[must_use]
    pub fn in_(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        self.in_round(n, self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Iₙ(self)` with explicit result precision.
    pub fn in_round(
        &self,
        n: i32,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_i_kernel(n, self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `I₀(self)` for `FixedFloat`. Delegates to [`BigFloat::i0`].
    #[must_use]
    pub fn i0(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().i0(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `I₁(self)` for `FixedFloat`. Delegates to [`BigFloat::i1`].
    #[must_use]
    pub fn i1(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().i1(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Iₙ(self)` for `FixedFloat`. Delegates to [`BigFloat::in_`].
    #[must_use]
    pub fn in_(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().in_(n, mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

/// `Iₙ(x)` for integer order `n` and real `x`.
///
/// Special values are handled directly per the module-level domain
/// table; for a normal argument the order is reduced to `I_m(|x|)`,
/// `m = |n|`, with one argument-parity sign (`Iₙ(−x) = (−1)ⁿ Iₙ(x)`;
/// `I₋ₙ = Iₙ` adds no sign), then the regime evaluator runs.
fn bessel_i_kernel(
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
        Class::Zero { .. } => {
            // I₀(±0) = 1, Iₙ(±0) = 0 for n ≠ 0 (DLMF 10.30.1); exact,
            // both zero signs alike (I is entire).
            let value = if m == 0 {
                BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1")
            } else {
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1")
            };
            (value, Status::OK)
        }
        Class::Infinity { sign } => {
            // Iₙ(+∞) = +∞; Iₙ(−∞) = (−1)ⁿ·∞ (the argument parity, so
            // −∞ for odd m). A genuine infinite limit (DLMF 10.30.4
            // eˣ/√(2πx) → ∞), Status::OK, the exp(+∞) precedent; not
            // the decaying-envelope convention.
            let neg = matches!(sign, Sign::Negative) && (m % 2 == 1);
            let result_sign = if neg { Sign::Negative } else { Sign::Positive };
            let inf =
                BigFloat::try_new_infinity(result_sign, target_precision).expect("precision >= 1");
            (inf, Status::OK)
        }
        Class::Normal { .. } => {
            // I is entire: no domain restriction. Reduce to I_m(|x|)
            // with one argument-parity sign. Iₙ(−x) = (−1)ⁿ Iₙ(x);
            // I₋ₙ(x) = Iₙ(x) (no order sign). Negate exactly when m is
            // odd and x < 0.
            let negate = (m % 2 == 1) && x.is_sign_negative();
            let ax = x.abs();

            let value = bessel_i_eval_normal(m, &ax, target_precision);
            let value = if negate { value.negated() } else { value };

            let (rounded, status) = value
                .round_to_precision(target_precision, mode)
                .expect("precision >= 1");
            auto_raise(status);
            (rounded, status)
        }
    }
}

/// `I_m(ax)` for `m ≥ 0`, normal `ax > 0`: the regime evaluator.
/// Returns the unrounded working-precision value; [`bessel_i_kernel`]
/// does the single final round.
///
/// Slice 6q.1 placeholder: returns a quiet NaN so the skeleton
/// compiles and the special-value path is exercised in isolation.
/// Slices 6q.2 (Maclaurin), 6q.3 (Miller + sum rule), and 6q.4
/// (asymptotic + regime dispatch) replace this with the real
/// evaluator; no normal-argument tests run until then.
fn bessel_i_eval_normal(_m: u32, ax: &BigFloat, target_precision: u32) -> BigFloat {
    let _ = ax;
    BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[]).expect("precision >= 1")
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    #[test]
    fn i0_zero_is_one() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, st) = z.i0(RoundingMode::NearestEven);
            let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
            assert_eq!(r.partial_cmp(&one).0, Some(Ordering::Equal), "I0(±0) = 1");
            assert!(!st.invalid());
        }
    }

    #[test]
    fn in_zero_is_zero_for_nonzero_order() {
        for n in [1i32, 2, 3, -1, -4] {
            for s in [Sign::Positive, Sign::Negative] {
                let z = BigFloat::try_new_zero(s, 53).unwrap();
                let (r, _) = z.in_(n, RoundingMode::NearestEven);
                assert!(r.is_zero(), "I_{n}(±0) = 0");
            }
        }
    }

    #[test]
    fn i_positive_infinity_is_positive_infinity() {
        let inf = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        for n in [0i32, 1, 2, 5, -3] {
            let (r, st) = inf.in_(n, RoundingMode::NearestEven);
            assert!(r.is_infinite() && r.is_sign_positive(), "I_{n}(+∞) = +∞");
            assert!(!st.invalid());
        }
    }

    #[test]
    fn i_negative_infinity_is_signed_by_parity() {
        let ninf = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        // Even order: I_n(−∞) = +∞. Odd order: −∞.
        for n in [0i32, 2, 4, -2] {
            let (r, st) = ninf.in_(n, RoundingMode::NearestEven);
            assert!(r.is_infinite() && r.is_sign_positive(), "I_{n}(−∞) = +∞");
            assert!(!st.invalid());
        }
        for n in [1i32, 3, -1, -3] {
            let (r, st) = ninf.in_(n, RoundingMode::NearestEven);
            assert!(r.is_infinite() && r.is_sign_negative(), "I_{n}(−∞) = −∞");
            assert!(!st.invalid());
        }
    }

    #[test]
    fn i_quiet_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.in_(2, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn i_signaling_nan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, st) = sn.i0(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(st.invalid());
    }

    #[test]
    fn precision_zero_is_rejected() {
        let x = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert!(x.i0_round(0, RoundingMode::NearestEven).is_err());
        assert!(x.in_round(3, 0, RoundingMode::NearestEven).is_err());
    }
}
