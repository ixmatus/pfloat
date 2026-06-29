//! `gamma(x) = Γ(x)`: the gamma function. Defined for all real `x`
//! except the non-positive integers (simple poles).
//!
//! Implementation: compose `lgamma` with `exp`, applying the
//! correct sign. For positive `x`, `Γ(x) > 0`, so the sign is
//! always positive. For negative non-integer `x`, the sign
//! alternates with the floor of `|x|`: `Γ(x) > 0` for
//! `x ∈ (−2, −1) ∪ (−4, −3) ∪ …` and `Γ(x) < 0` for
//! `x ∈ (−1, 0) ∪ (−3, −2) ∪ …`. Concretely, the sign equals
//! `−sign(sin(πx))` on the negative reals.
//!
//! For exact positive integer inputs `n ≥ 1`, the kernel returns
//! `(n−1)!` exactly when the result fits in target precision.
//! Otherwise it goes through `exp(lgamma(x))`, which inherits the
//! `lgamma` precision.
//!
//! Special cases per IEEE 754-2019 §9.4:
//!
//! - `gamma(+0) = +∞ + DIV_BY_ZERO`.
//! - `gamma(−0) = −∞ + DIV_BY_ZERO`.
//! - `gamma(negative integer) = qNaN + INVALID` (pole).
//! - `gamma(+∞) = +∞`.
//! - `gamma(−∞) = qNaN + INVALID`.
//! - `gamma(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::lgamma::is_integer_test;
use super::pi_at;
use super::ziv::ziv_round;
use super::ziv_calibration::GAMMA_ERROR_GUARD;

impl BigFloat {
    /// `gamma(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn gamma(&self, mode: RoundingMode) -> (Self, Status) {
        self.gamma_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `gamma(self)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.29, ADR-0038).
    pub fn gamma_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(gamma_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `gamma(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::gamma`].
    #[must_use]
    pub fn gamma(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().gamma(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn gamma_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            // gamma(±0) = ±∞ + DIV_BY_ZERO.
            let inf = BigFloat::try_new_infinity(*sign, target_precision).expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            return (inf, Status::DIV_BY_ZERO);
        }
        Class::Infinity {
            sign: Sign::Positive,
        } => {
            return (
                BigFloat::try_new_infinity(Sign::Positive, target_precision)
                    .expect("precision >= 1"),
                Status::OK,
            );
        }
        Class::Infinity {
            sign: Sign::Negative,
        } => {
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        Class::Normal { .. } => {}
    }

    // Negative integer: pole, qNaN + INVALID.
    if matches!(x.sign(), Sign::Negative) && is_integer_test(x) {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }

    // Positive-integer exact fast path. For positive integer `n`
    // where `(n−1)!` is exactly representable in target_precision
    // bits, gamma(n) = (n−1)! exactly under any rounding mode. The
    // exp(lgamma) chain returns a value epsilon-away from the exact
    // factorial; under NE that rounds correctly, but under directed
    // modes the epsilon's sign tips the rounding to the adjacent
    // ULP. The exact factorial dispatch is mode-independent and
    // sidesteps the Ziv interval test's inability to certify
    // exactly-representable true values.
    if let Some(exact) = try_gamma_pos_integer_exact(x, target_precision) {
        return (exact, Status::OK);
    }

    // Ziv-driven correct rounding under every IEEE mode. The eval
    // closure composes lgamma (Ziv-driven internally per the
    // already-Ziv cohort + Spouge dispatch from pf-l6s5) and exp
    // (Ziv-driven, slice p1.2) at working precision `w`, then
    // applies the sign from gamma_sign_of. gamma_sign_of returns a
    // binary Sign that is robust away from negative-integer poles
    // (those are handled by the special-case dispatch above), so
    // the sign does not participate in the Ziv interval test.
    let (result, status) = ziv_round(
        |w| {
            let (ln_abs_gamma, _) = x
                .lgamma_round(w, RoundingMode::NearestEven)
                .expect("precision >= 1");
            let (abs_gamma, _) = ln_abs_gamma.exp(RoundingMode::NearestEven);
            let result_sign = gamma_sign_of(x, w);
            if matches!(result_sign, Sign::Negative) {
                abs_gamma.negated()
            } else {
                abs_gamma
            }
        },
        target_precision,
        mode,
        GAMMA_ERROR_GUARD,
    );
    // Defensive INEXACT guard (pf-umlm, ADR-0066): a finite-normal
    // fall-through (a non-integer x; the integer points are dispatched
    // above) is irrational. The ADR-0065 sweep showed this path already
    // flags INEXACT everywhere, so the force is a no-op hardening against
    // regression; its worst-case soundness rests on the irrationality of
    // Γ at dyadic non-integers, which is not proven for every dyadic.
    let status = super::force_transcendental_inexact(&result, status);
    auto_raise(status);
    (result, status)
}

/// If `x` is a positive integer `n` and `(n−1)!` is exactly
/// representable in `target_precision` bits, returns `(n−1)!` at
/// target precision; otherwise returns `None`. The fast path
/// sidesteps the `exp(lgamma)` chain's noise floor that tips
/// directed-mode rounding off the exact factorial.
///
/// The iteration accumulates `1 · 2 · 3 · … · (n−1)` at
/// `target_precision`; the first multiplication that rounds (the
/// `Status::INEXACT` flag) signals that the running factorial no
/// longer fits exactly, so the function returns `None` and the
/// caller falls through to the Ziv envelope. The cap on iterations
/// is a safety bound — even at `target_precision = 4096`, only
/// `n ≲ 500` factorials fit, so 512 is generous.
fn try_gamma_pos_integer_exact(x: &BigFloat, target_precision: u32) -> Option<BigFloat> {
    if !matches!(x.sign(), Sign::Positive) || x.is_zero() || !is_integer_test(x) {
        return None;
    }
    let one = BigFloat::try_from_i64_exact(1, target_precision).ok()?;
    let mut acc = one.clone();
    let mut k = BigFloat::try_from_i64_exact(2, target_precision).ok()?;
    for _ in 0..512 {
        // Stop when k ≥ x (so the loop multiplies by 2, 3, …, n−1).
        match k.partial_cmp(x).0 {
            Some(core::cmp::Ordering::Less) => {}
            _ => return Some(acc),
        }
        let (next_acc, status) = acc.mul(&k, RoundingMode::NearestEven);
        if status.inexact() {
            return None;
        }
        acc = next_acc;
        let (next_k, _) = k.add(&one, RoundingMode::NearestEven);
        k = next_k;
    }
    None
}

/// Sign of `Γ(x)` for finite non-zero `x` that is not a negative
/// integer pole. Positive `x` always yields positive sign; negative
/// non-integer `x` alternates per the reflection
/// `Γ(x)·Γ(1−x) = π/sin(πx)` (`Γ(1−x) > 0` for `x < 1`, so the
/// sign of `Γ(x)` matches the sign of `sin(πx)` for negative `x`).
pub(super) fn gamma_sign_of(x: &BigFloat, working_prec: u32) -> Sign {
    if matches!(x.sign(), Sign::Positive) {
        return Sign::Positive;
    }
    let pi = pi_at(working_prec);
    let (pi_x, _) = pi.mul(x, RoundingMode::NearestEven);
    let (sin_val, _) = pi_x.sin(RoundingMode::NearestEven);
    if matches!(sin_val.sign(), Sign::Negative) {
        Sign::Negative
    } else {
        Sign::Positive
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    fn close_at(v: &BigFloat, expected: &BigFloat, bits: u32) -> bool {
        let (diff, _) = v.sub(expected, RoundingMode::NearestEven);
        let abs_diff = diff.abs();
        if abs_diff.is_zero() {
            return true;
        }
        let p = v.precision().max(expected.precision());
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let abs_b = expected.abs();
        let mut bound = if abs_b.is_zero() { one } else { abs_b };
        for _ in 0..bits {
            bound = bound.div(&two, RoundingMode::NearestEven).0;
        }
        matches!(
            abs_diff.partial_cmp(&bound).0,
            Some(Ordering::Less | Ordering::Equal)
        )
    }

    #[test]
    fn gamma_one_is_one() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.gamma(RoundingMode::NearestEven);
        let expected = BigFloat::try_from_i64_exact(1, 113).unwrap();
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn gamma_two_is_one() {
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (r, _) = two.gamma(RoundingMode::NearestEven);
        let expected = BigFloat::try_from_i64_exact(1, 113).unwrap();
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn gamma_five_is_twentyfour() {
        // Γ(5) = 4! = 24.
        let five = BigFloat::try_from_i64_exact(5, 113).unwrap();
        let (r, _) = five.gamma(RoundingMode::NearestEven);
        let expected = BigFloat::try_from_i64_exact(24, 113).unwrap();
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn gamma_half_is_sqrt_pi() {
        // Γ(1/2) = √π ≈ 1.7724538509055160272981674833411.
        let half = BigFloat::parse_str("0.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = half.gamma(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "1.7724538509055160272981674833411451827975494561224",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn gamma_pos_zero_is_pos_inf_div() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, status) = z.gamma(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
        assert!(status.div_by_zero());
    }

    #[test]
    fn gamma_neg_zero_is_neg_inf_div() {
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, status) = nz.gamma(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative());
        assert!(status.div_by_zero());
    }

    #[test]
    fn gamma_negative_integer_is_nan() {
        let neg_three = BigFloat::try_from_i64_exact(-3, 53).unwrap();
        let (r, status) = neg_three.gamma(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn gamma_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.gamma(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn gamma_neg_inf_is_invalid() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, status) = ni.gamma(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn gamma_neg_half_is_neg_two_sqrt_pi() {
        // Γ(-1/2) = -2√π ≈ -3.5449077018110320546.
        let neg_half = BigFloat::parse_str("-0.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = neg_half.gamma(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "-3.5449077018110320545963349666822903655950989122448",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 60));
    }

    #[test]
    fn gamma_neg_1_5_is_positive() {
        // Γ(-1.5) = 4√π/3 ≈ 2.36327180120735... (positive).
        let neg = BigFloat::parse_str("-1.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = neg.gamma(RoundingMode::NearestEven);
        assert!(r.is_sign_positive());
        let expected = BigFloat::parse_str(
            "2.3632718012073547030642233111215269103967326081632",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 60));
    }

    #[test]
    fn gamma_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.gamma(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn gamma_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.gamma(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
