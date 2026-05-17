//! Bessel functions of the first kind `J0`, `J1`, `Jn` (DLMF
//! Chapter 10): ordinary Bessel, integer order, real argument.
//! Entire on the real line: no poles, no domain restriction.
//!
//! Order and sign are reduced before evaluation. `Jₙ(−x) =
//! (−1)ⁿ Jₙ(x)` (DLMF 10.11.1) and `J₋ₙ(x) = (−1)ⁿ Jₙ(x)` (DLMF
//! 10.4.1), so the kernel evaluates `J_m(|x|)` for `m = |n| ≥ 0`
//! and applies one parity sign.
//!
//! Three regimes, dispatched on the binary exponent of `|x|` (the
//! [`super::airy`] / [`super::si`] integer-exponent selector idiom);
//! the regime evaluators land in slices 6o.2 and 6o.3:
//!
//! - Tiny `|x|`: the leading Maclaurin terms (DLMF 10.2.2). Keeps
//!   the `2k/x` recurrence away from `x → 0`.
//! - Moderate `|x|`: Miller backward recurrence
//!   `f_{k−1} = (2k/x)·f_k − f_{k+1}` normalised by the sum rule
//!   `J₀ + 2·Σ J_{2k} = 1` (DLMF 10.6.1, 10.12.4). ADR-0023.
//! - Large `|x|`: the Hankel-form asymptotic (DLMF 10.17.3) summed
//!   to its smallest term.
//!
//! Special cases:
//!
//! - `J₀(±0) = 1`, `Jₙ(±0) = 0` for `n ≠ 0` (exact, DLMF 10.2.2).
//! - `Jₙ(±∞) = +0` for every order, by the decaying-envelope
//!   convention (ADR-0021, the [`super::airy`] precedent): the true
//!   behaviour at `±∞` is a bounded decaying oscillation with no
//!   limit; the conservative total result is `+0`, `Status::OK`.
//! - `Jₙ(NaN) = NaN`; `sNaN` raises `INVALID`.

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
    /// `J₀(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn j0(&self, mode: RoundingMode) -> (Self, Status) {
        self.j0_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `J₀(self)` with explicit result precision.
    pub fn j0_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_j_kernel(0, self, target_precision, mode))
    }

    /// `J₁(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn j1(&self, mode: RoundingMode) -> (Self, Status) {
        self.j1_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `J₁(self)` with explicit result precision.
    pub fn j1_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_j_kernel(1, self, target_precision, mode))
    }

    /// `Jₙ(self)` for integer order `n` (any `i32`, including
    /// negative: `J₋ₙ = (−1)ⁿ Jₙ`), rounded under `mode` to
    /// `self.precision`.
    #[must_use]
    pub fn jn(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        self.jn_round(n, self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Jₙ(self)` with explicit result precision.
    pub fn jn_round(
        &self,
        n: i32,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(bessel_j_kernel(n, self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `J₀(self)` for `FixedFloat`. Delegates to [`BigFloat::j0`].
    #[must_use]
    pub fn j0(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().j0(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `J₁(self)` for `FixedFloat`. Delegates to [`BigFloat::j1`].
    #[must_use]
    pub fn j1(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().j1(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Jₙ(self)` for `FixedFloat`. Delegates to [`BigFloat::jn`].
    #[must_use]
    pub fn jn(&self, n: i32, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().jn(n, mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

/// `Jₙ(x)` for integer order `n` and real `x`.
///
/// Special values are handled directly; for a normal argument the
/// order and sign are reduced to `J_m(|x|)`, `m = |n|`, with one
/// parity sign, then the regime evaluator runs.
fn bessel_j_kernel(
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
            return (nan, Status::OK);
        }
        Class::Zero { .. } => {
            // J₀(±0) = 1, Jₙ(±0) = 0 for n ≠ 0 (DLMF 10.2.2); exact.
            let value = if m == 0 {
                BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1")
            } else {
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1")
            };
            return (value, Status::OK);
        }
        Class::Infinity { .. } => {
            // Decaying-envelope convention (ADR-0021): bounded
            // oscillation with no limit → +0, Status::OK.
            let zero =
                BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1");
            return (zero, Status::OK);
        }
        Class::Normal { .. } => {}
    }

    // Normal argument: reduce to J_m(|x|) with one parity sign.
    // Jₙ(−x) = (−1)ⁿ Jₙ(x); J₋ₙ(x) = (−1)ⁿ Jₙ(x). Each negative
    // contributes (−1)^m, so the result is negated exactly when m is
    // odd and exactly one of {n<0, x<0} holds.
    let negate = (m % 2 == 1) && ((n < 0) ^ x.is_sign_negative());
    let ax = x.abs();

    let value = bessel_j_eval_normal(m, &ax, target_precision);
    let value = if negate { value.negated() } else { value };

    let (rounded, status) = value
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    auto_raise(status);
    (rounded, status)
}

/// `J_m(ax)` for `m ≥ 0`, `ax ≥ 0`, normal. The three-regime
/// evaluator (tiny / Miller recurrence / Hankel asymptotic).
///
/// Slice 6o.1 placeholder: returns a quiet NaN so the skeleton
/// compiles and the special-value path is exercised in isolation.
/// Slices 6o.2 (tiny + Miller) and 6o.3 (asymptotic) replace this
/// with the real dispatch; no normal-argument tests run until then.
fn bessel_j_eval_normal(_m: u32, ax: &BigFloat, target_precision: u32) -> BigFloat {
    let _ = ax;
    BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[]).expect("precision >= 1")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn j0_zero_is_one() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, s) = z.j0(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert_eq!(r.partial_cmp(&one).0, Some(core::cmp::Ordering::Equal));
        assert!(!s.invalid());
    }

    #[test]
    fn jn_zero_is_zero_for_nonzero_order() {
        for n in [1i32, 2, 3, -1, -4] {
            let z = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
            let (r, _) = z.jn(n, RoundingMode::NearestEven);
            assert!(r.is_zero(), "J_{n}(0) should be 0");
        }
        // J₀(−0) is still 1.
        let z = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, _) = z.j1(RoundingMode::NearestEven);
        assert!(r.is_zero());
    }

    #[test]
    fn jn_infinity_is_zero() {
        for n in [0i32, 1, 5, -3] {
            for s in [Sign::Positive, Sign::Negative] {
                let inf = BigFloat::try_new_infinity(s, 53).unwrap();
                let (r, st) = inf.jn(n, RoundingMode::NearestEven);
                assert!(r.is_zero() && r.is_sign_positive(), "J_{n}(±∞) = +0");
                assert!(!st.invalid());
            }
        }
    }

    #[test]
    fn jn_quiet_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.jn(2, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn jn_signaling_nan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, st) = sn.j0(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(st.invalid());
    }

    #[test]
    fn precision_zero_is_rejected() {
        let x = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert!(x.j0_round(0, RoundingMode::NearestEven).is_err());
        assert!(x.jn_round(3, 0, RoundingMode::NearestEven).is_err());
    }
}
