//! Airy functions `Ai`, `Bi` and their derivatives `Ai′`, `Bi′`
//! (DLMF Chapter 9). All four are entire on the real line: no
//! poles, no domain restriction.
//!
//! Three regimes, dispatched on the binary exponent of `|x|` and
//! its sign (the [`super::erf`] / [`super::si`] dispatch idiom
//! extended with a sign-aware asymptotic split):
//!
//! - Small `|x|`: the Maclaurin series (DLMF 9.4.1–9.4.6) in the
//!   two entire solutions `f`, `g`, combined with the boundary
//!   constants `Ai(0)`, `Ai′(0)` (DLMF 9.2.3–9.2.6).
//! - Large positive `x`: the exponential asymptotic (DLMF
//!   9.7.5–9.7.8) in `ζ = (2/3)·x^{3/2}`, summed to its smallest
//!   term.
//! - Large negative `x`: the oscillatory asymptotic (DLMF
//!   9.7.9–9.7.12) in `ζ = (2/3)·|x|^{3/2}` with the phase
//!   `ξ = ζ + π/4`.
//!
//! Special cases:
//!
//! - `Ai(±0)`, `Bi(±0)`, `Ai′(±0)`, `Bi′(±0)`: the exact boundary
//!   constants (finite, normal).
//! - `Ai(+∞) = +0`, `Ai′(+∞) = −0`, `Bi(+∞) = +∞`, `Bi′(+∞) = +∞`
//!   (the exact limits at an infinite argument, `Status::OK`, the
//!   `exp(+∞)`/`gamma(+∞)` convention).
//! - `Ai(−∞) = Bi(−∞) = Ai′(−∞) = Bi′(−∞) = +0` by the
//!   decaying-envelope convention: the true behaviour at `−∞` is a
//!   bounded oscillation with no limit; the conservative total
//!   result is `+0` with `Status::OK`. ADR-0021 records this.
//! - `f(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

/// Which Airy function a kernel invocation evaluates. The four share
/// the boundary constants, the `f`/`g` Maclaurin series, the
/// `u_k`/`v_k` asymptotic coefficient recurrence, and `ζ`/`x^{1/4}`,
/// so the kernel is parameterised rather than duplicated.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum AiryFn {
    Ai,
    Bi,
    AiPrime,
    BiPrime,
}

impl BigFloat {
    /// `Ai(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn ai(&self, mode: RoundingMode) -> (Self, Status) {
        self.ai_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Ai(self)` with explicit result precision.
    pub fn ai_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(airy_kernel(AiryFn::Ai, self, target_precision, mode))
    }

    /// `Bi(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn bi(&self, mode: RoundingMode) -> (Self, Status) {
        self.bi_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Bi(self)` with explicit result precision.
    pub fn bi_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(airy_kernel(AiryFn::Bi, self, target_precision, mode))
    }

    /// `Ai′(self)` (derivative of `Ai`) rounded under `mode` to
    /// `self.precision`.
    #[must_use]
    pub fn ai_prime(&self, mode: RoundingMode) -> (Self, Status) {
        self.ai_prime_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Ai′(self)` with explicit result precision.
    pub fn ai_prime_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(airy_kernel(AiryFn::AiPrime, self, target_precision, mode))
    }

    /// `Bi′(self)` (derivative of `Bi`) rounded under `mode` to
    /// `self.precision`.
    #[must_use]
    pub fn bi_prime(&self, mode: RoundingMode) -> (Self, Status) {
        self.bi_prime_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `Bi′(self)` with explicit result precision.
    pub fn bi_prime_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(airy_kernel(AiryFn::BiPrime, self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `Ai(self)` for `FixedFloat`. Delegates to [`BigFloat::ai`].
    #[must_use]
    pub fn ai(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().ai(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Bi(self)` for `FixedFloat`. Delegates to [`BigFloat::bi`].
    #[must_use]
    pub fn bi(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().bi(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Ai′(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::ai_prime`].
    #[must_use]
    pub fn ai_prime(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().ai_prime(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }

    /// `Bi′(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::bi_prime`].
    #[must_use]
    pub fn bi_prime(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().bi_prime(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn airy_kernel(
    which: AiryFn,
    x: &BigFloat,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
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
            // f(±0) = the exact boundary constant (finite, normal).
            let value = airy_zero_value(which, target_precision);
            let (rounded, status) = value
                .round_to_precision(target_precision, mode)
                .expect("precision >= 1");
            auto_raise(status);
            return (rounded, status);
        }
        Class::Infinity { sign } => {
            let result = airy_at_infinity(which, *sign, target_precision);
            return (result, Status::OK);
        }
        Class::Normal { .. } => {}
    }

    let value = airy_eval_normal(which, x, target_precision);
    let (rounded, status) = value
        .round_to_precision(target_precision, mode)
        .expect("precision >= 1");
    auto_raise(status);
    (rounded, status)
}

/// The exact limit of an Airy function at `±∞` (DLMF 9.7). At `+∞`:
/// `Ai → +0`, `Ai′ → −0`, `Bi → +∞`, `Bi′ → +∞` (true limits at an
/// infinite argument, mirroring `exp(+∞)`). At `−∞` all four → `+0`
/// by the decaying-envelope convention (ADR-0021): the true
/// behaviour is a bounded oscillation with no limit; the
/// conservative total result is `+0`.
fn airy_at_infinity(which: AiryFn, sign: Sign, target_precision: u32) -> BigFloat {
    if matches!(sign, Sign::Negative) {
        return BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1");
    }
    match which {
        AiryFn::Ai => {
            BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1")
        }
        AiryFn::AiPrime => {
            BigFloat::try_new_zero(Sign::Negative, target_precision).expect("precision >= 1")
        }
        AiryFn::Bi | AiryFn::BiPrime => {
            BigFloat::try_new_infinity(Sign::Positive, target_precision).expect("precision >= 1")
        }
    }
}

/// The boundary constant `f(0)` for `f ∈ {Ai, Bi, Ai′, Bi′}` (DLMF
/// 9.2.3–9.2.6). Implemented in beads issue pf-w8l.
fn airy_zero_value(_which: AiryFn, _working_prec: u32) -> BigFloat {
    unimplemented!("slice-6n: boundary constants — beads pf-w8l")
}

/// Evaluate an Airy function at a finite non-zero argument via the
/// three-regime sign-aware dispatch (Maclaurin / exponential
/// asymptotic / oscillatory asymptotic). Implemented in beads issues
/// pf-85k (series) and pf-nfb (asymptotic).
fn airy_eval_normal(_which: AiryFn, _x: &BigFloat, _target_precision: u32) -> BigFloat {
    unimplemented!("slice-6n: three-regime dispatch — beads pf-85k / pf-nfb")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn airy_nan_propagates() {
        for which in [AiryFn::Ai, AiryFn::Bi, AiryFn::AiPrime, AiryFn::BiPrime] {
            let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
            let (r, status) = airy_kernel(which, &q, 53, RoundingMode::NearestEven);
            assert!(r.is_quiet_nan(), "{which:?} qNaN");
            assert!(status.is_ok(), "{which:?} qNaN status");
        }
    }

    #[test]
    fn airy_snan_raises_invalid() {
        for which in [AiryFn::Ai, AiryFn::Bi, AiryFn::AiPrime, AiryFn::BiPrime] {
            let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
            let (r, status) = airy_kernel(which, &sn, 53, RoundingMode::NearestEven);
            assert!(r.is_quiet_nan(), "{which:?} sNaN");
            assert!(status.invalid(), "{which:?} sNaN INVALID");
        }
    }

    #[test]
    fn airy_pos_inf_limits() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (ai, _) = pi.ai(RoundingMode::NearestEven);
        assert!(ai.is_zero() && !ai.is_sign_negative(), "Ai(+∞) = +0");
        let (aip, _) = pi.ai_prime(RoundingMode::NearestEven);
        assert!(aip.is_zero() && aip.is_sign_negative(), "Ai′(+∞) = −0");
        let (bi, _) = pi.bi(RoundingMode::NearestEven);
        assert!(bi.is_infinite() && !bi.is_sign_negative(), "Bi(+∞) = +∞");
        let (bip, _) = pi.bi_prime(RoundingMode::NearestEven);
        assert!(bip.is_infinite() && !bip.is_sign_negative(), "Bi′(+∞) = +∞");
    }

    #[test]
    fn airy_neg_inf_is_pos_zero_by_convention() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        for (got, name) in [
            (ni.ai(RoundingMode::NearestEven).0, "Ai"),
            (ni.bi(RoundingMode::NearestEven).0, "Bi"),
            (ni.ai_prime(RoundingMode::NearestEven).0, "Ai′"),
            (ni.bi_prime(RoundingMode::NearestEven).0, "Bi′"),
        ] {
            assert!(
                got.is_zero() && !got.is_sign_negative(),
                "{name}(−∞) = +0 by the decaying-envelope convention"
            );
        }
    }

    #[test]
    fn airy_fn_enum_is_copy() {
        // Guards the parameterised-kernel design: AiryFn must stay a
        // trivial Copy tag so the four entry points share one kernel.
        let a = AiryFn::Ai;
        let b = a;
        assert_eq!(a, b);
        assert_ne!(AiryFn::Ai, AiryFn::Bi);
    }
}
