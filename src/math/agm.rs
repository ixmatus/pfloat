//! `agm(a, b)`: the arithmetic-geometric mean of `a` and `b`.
//!
//! The Gauss iteration:
//!
//! ```text
//! a_{n+1} = (a_n + b_n) / 2     (arithmetic mean)
//! b_{n+1} = sqrt(a_n · b_n)      (geometric mean)
//! ```
//!
//! starting from `a_0 = a` and `b_0 = b`, converges quadratically to
//! a common limit: the AGM. `O(log p)` iterations suffice at working
//! precision `p`; the loop terminates once `|a_n − b_n|` falls below
//! `2^(−p_work − 4)`.
//!
//! The kernel computes at a working precision of
//! `target_precision + 64` bits, then rounds back. The 64-bit guard
//! absorbs the per-iteration rounding error of three operations
//! (one add, one mul, one sqrt) compounded over `O(log p_work)`
//! iterations; for any precision pfloat supports (up to
//! `u32::MAX − 64`), this is well over twice the worst-case
//! accumulated error bound for the iteration.
//!
//! Domain: AGM is defined for non-negative real `a` and `b`. The
//! geometric mean of a negative operand is not real, so negative
//! finite operands raise `INVALID` and return qNaN.
//!
//! Special cases:
//!
//! - `agm(NaN, _) = agm(_, NaN) = NaN`; `sNaN` raises `INVALID`.
//! - `agm(negative_finite, _) = agm(_, negative_finite) = qNaN +
//!   INVALID`.
//! - `agm(+0, x) = agm(x, +0) = +0` for `x >= 0` (the geometric
//!   mean kills the sequence after one step). For `x = +∞` the
//!   iteration does not converge; return `qNaN + INVALID`.
//! - `agm(+∞, +∞) = +∞`.
//! - `agm(+∞, finite_positive) = +∞` (`a_n` stays `+∞`; `b_n`
//!   grows without bound at rate `sqrt(+∞ · finite_positive)`).
//! - `agm(x, x) = x` (fixed point).
//!
//! ADR-0015 records the choice of Gauss's iteration over
//! Brent-Salamin's variant (the latter is a specialization for `π`
//! computation that layers extra bookkeeping on top of AGM; the
//! standalone AGM kernel benefits from neither).

use core::cmp::Ordering;

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::ziv::ziv_round;
use super::ziv_calibration::AGM_ERROR_GUARD;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `agm(self, other)` rounded under `mode` to
    /// `max(self.precision, other.precision)`.
    #[must_use]
    pub fn agm(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let target = self.precision.max(other.precision);
        self.agm_round(other, target, mode)
            .expect("max of two valid precisions is valid")
    }

    /// `agm(self, other)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.36, ADR-0038).
    pub fn agm_round(
        &self,
        other: &Self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(agm_kernel(self, other, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `agm(self, other)` for `FixedFloat`. Delegates to
    /// [`BigFloat::agm`].
    #[must_use]
    pub fn agm(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().agm(&other.to_big(), mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn agm_kernel(
    a: &BigFloat,
    b: &BigFloat,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    // Signaling-NaN check first: any sNaN operand raises INVALID.
    if a.is_signaling_nan() || b.is_signaling_nan() {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }
    // Quiet-NaN propagation. The sign of the produced NaN follows the
    // first NaN operand to match the rest of the surface.
    if a.is_nan() {
        let nan =
            BigFloat::try_new_quiet_nan(a.sign(), target_precision, &[]).expect("precision >= 1");
        return (nan, Status::OK);
    }
    if b.is_nan() {
        let nan =
            BigFloat::try_new_quiet_nan(b.sign(), target_precision, &[]).expect("precision >= 1");
        return (nan, Status::OK);
    }

    // Negative finite or negative infinity operand: AGM is undefined.
    if is_negative(a) || is_negative(b) {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }

    // `agm(+∞, +∞) = +∞`; `agm(+∞, finite_positive) = +∞` (the AM
    // stays infinite and the GM grows without bound). The mixed
    // `agm(+∞, +0)` case does not converge; flag it.
    if a.is_infinite() || b.is_infinite() {
        if (a.is_infinite() && b.is_zero()) || (b.is_infinite() && a.is_zero()) {
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        let inf =
            BigFloat::try_new_infinity(Sign::Positive, target_precision).expect("precision >= 1");
        return (inf, Status::OK);
    }

    // Either zero short-circuits to +0: the geometric mean of zero
    // and anything finite is zero, after which the arithmetic mean
    // halves on every step and the b sequence stays at zero.
    if a.is_zero() || b.is_zero() {
        let z = BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1");
        return (z, Status::OK);
    }

    // Equal-argument fixed-point dispatch (pf-kk16, ADR-0039).
    // agm(x, x) = x exactly: the iteration is at its fixed point,
    // so every subsequent (a_n, b_n) pair equals (x, x). Without
    // the dispatch, the Ziv-wrapped iteration would return x +
    // epsilon (from the round-to-w step) and tip rounding under
    // directed modes off the exact value
    // (`feedback_exact_value_defeats_ziv`).
    if matches!(a.partial_cmp(b).0, Some(Ordering::Equal)) {
        let (rounded, status) = a
            .round_to_precision(target_precision, mode)
            .expect("precision >= 1");
        auto_raise(status);
        return (rounded, status);
    }

    // Both operands are now finite, strictly positive, and unequal.
    // Ziv-driven correct rounding under every IEEE mode: the eval
    // closure captures (a, b) and runs the Gauss AGM iteration at
    // working precision w. Quadratic convergence doubles the bit
    // agreement each step, so Ziv adds at most one extra iteration
    // per retry (O(log w) → O(log 2w) ≈ +1).
    let (result, status) = ziv_round(
        |w| {
            let mut a_n = a
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let mut b_n = b
                .round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;

            // Canonicalize so a_n >= b_n. The iteration is symmetric
            // in (a, b), so this only tidies the convergence trace.
            if matches!(a_n.partial_cmp(&b_n).0, Some(Ordering::Less)) {
                core::mem::swap(&mut a_n, &mut b_n);
            }

            let max_iter = 64u32;
            let convergence_exponent_floor = -i64::from(w) - 4;
            let two = BigFloat::try_from_i64_exact(2, w).expect("precision >= 1");

            for _ in 0..max_iter {
                let (diff, _) = a_n.sub(&b_n, RoundingMode::NearestEven);
                let abs_diff = diff.abs();
                let converged = match &abs_diff.class {
                    Class::Zero { .. } => true,
                    Class::Normal { exponent, .. } => *exponent < convergence_exponent_floor,
                    _ => false,
                };
                if converged {
                    break;
                }

                let (sum, _) = a_n.add(&b_n, RoundingMode::NearestEven);
                let (am, _) = sum.div(&two, RoundingMode::NearestEven);
                let (prod, _) = a_n.mul(&b_n, RoundingMode::NearestEven);
                let (gm, _) = prod.sqrt(RoundingMode::NearestEven);
                a_n = am;
                b_n = gm;
            }

            // After convergence the AM and GM agree to working
            // precision; averaging them absorbs any final 1-ULP
            // separation.
            let (sum, _) = a_n.add(&b_n, RoundingMode::NearestEven);
            sum.div(&two, RoundingMode::NearestEven).0
        },
        target_precision,
        mode,
        AGM_ERROR_GUARD,
    );
    auto_raise(status);
    (result, status)
}

fn is_negative(x: &BigFloat) -> bool {
    match &x.class {
        Class::Normal { sign, .. } | Class::Infinity { sign } => matches!(sign, Sign::Negative),
        _ => false,
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
    fn agm_x_x_is_x() {
        let x = BigFloat::try_from_i64_exact(7, 113).unwrap();
        let (r, _) = x.agm(&x, RoundingMode::NearestEven);
        assert_eq!(r.partial_cmp(&x).0, Some(Ordering::Equal));
    }

    #[test]
    fn agm_x_x_is_x_under_every_directed_mode() {
        // pf-kk16 pinning test: the equal-argument fixed-point
        // dispatch returns x exactly under every mode. Without the
        // dispatch, the Ziv iteration's round-to-w step would return
        // x + epsilon and tip directed-mode rounding off the exact
        // value (feedback_exact_value_defeats_ziv). Exercise an
        // integer (exactly representable, so the agm value is the
        // integer) and a non-integer-but-equal pair.
        for &mode in &[
            RoundingMode::NearestEven,
            RoundingMode::TowardPositive,
            RoundingMode::TowardNegative,
            RoundingMode::TowardZero,
            RoundingMode::NearestAway,
        ] {
            for &prec in &[24u32, 53, 113] {
                for &v in &[1i64, 7, 100, -0] {
                    let x = BigFloat::try_from_i64_exact(v.max(0), prec).unwrap();
                    let (r, status) = x.agm(&x, mode);
                    assert!(status.is_ok(), "agm({v},{v}) status under {mode:?}@p{prec}");
                    assert_eq!(
                        r.partial_cmp(&x).0,
                        Some(Ordering::Equal),
                        "agm({v},{v}) = {v} expected under {mode:?}@p{prec}, got {r:?}"
                    );
                    assert_eq!(r.precision(), prec);
                }
            }
        }
    }

    #[test]
    fn agm_zero_zero_is_zero() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, _) = z.agm(&z, RoundingMode::NearestEven);
        assert!(r.is_zero());
        assert!(r.is_sign_positive());
    }

    #[test]
    fn agm_zero_x_is_zero() {
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let x = BigFloat::try_from_i64_exact(5, 53).unwrap();
        let (r, _) = z.agm(&x, RoundingMode::NearestEven);
        assert!(r.is_zero());
        let (r2, _) = x.agm(&z, RoundingMode::NearestEven);
        assert!(r2.is_zero());
    }

    #[test]
    fn agm_one_two() {
        // Reference: agm(1, 2) ≈ 1.4567910310469068691864323832650819749738248292...
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (r, _) = one.agm(&two, RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "1.4567910310469068691864323832650819749738248292",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 100));
    }

    #[test]
    fn agm_is_symmetric() {
        let a = BigFloat::parse_str("3.7", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let b = BigFloat::parse_str("11.25", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (ab, _) = a.agm(&b, RoundingMode::NearestEven);
        let (ba, _) = b.agm(&a, RoundingMode::NearestEven);
        assert!(close_at(&ab, &ba, 100));
    }

    #[test]
    fn agm_step_invariance() {
        // agm(a, b) = agm((a + b) / 2, sqrt(a · b)).
        let a = BigFloat::parse_str("2.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let b = BigFloat::parse_str("9.0", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (sum, _) = a.add(&b, RoundingMode::NearestEven);
        let (am, _) = sum.div(&two, RoundingMode::NearestEven);
        let (prod, _) = a.mul(&b, RoundingMode::NearestEven);
        let (gm, _) = prod.sqrt(RoundingMode::NearestEven);
        let (direct, _) = a.agm(&b, RoundingMode::NearestEven);
        let (one_step, _) = am.agm(&gm, RoundingMode::NearestEven);
        assert!(close_at(&direct, &one_step, 100));
    }

    #[test]
    fn agm_negative_is_invalid() {
        let neg = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, status) = neg.agm(&one, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
        let (r2, status2) = one.agm(&neg, RoundingMode::NearestEven);
        assert!(r2.is_quiet_nan());
        assert!(status2.invalid());
    }

    #[test]
    fn agm_pos_inf_finite_is_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, _) = pi.agm(&one, RoundingMode::NearestEven);
        assert!(r.is_infinite());
        assert!(r.is_sign_positive());
    }

    #[test]
    fn agm_pos_inf_pos_inf_is_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.agm(&pi, RoundingMode::NearestEven);
        assert!(r.is_infinite());
        assert!(r.is_sign_positive());
    }

    #[test]
    fn agm_pos_inf_pos_zero_is_invalid() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let z = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, status) = pi.agm(&z, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn agm_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, _) = q.agm(&one, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        let (r2, _) = one.agm(&q, RoundingMode::NearestEven);
        assert!(r2.is_quiet_nan());
    }

    #[test]
    fn agm_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, status) = sn.agm(&one, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn agm_sandwich_property() {
        // For a > b > 0, b < agm(a, b) < a.
        let a = BigFloat::try_from_i64_exact(10, 113).unwrap();
        let b = BigFloat::try_from_i64_exact(3, 113).unwrap();
        let (r, _) = a.agm(&b, RoundingMode::NearestEven);
        assert_eq!(r.partial_cmp(&b).0, Some(Ordering::Greater));
        assert_eq!(r.partial_cmp(&a).0, Some(Ordering::Less));
    }

    #[test]
    fn agm_round_rejects_zero_precision() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert_eq!(
            one.agm_round(&two, 0, RoundingMode::NearestEven),
            Err(BuildError::PrecisionZero)
        );
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixed_agm() {
        let one = FixedFloat::<113>::try_from_i64_exact(1).unwrap();
        let two = FixedFloat::<113>::try_from_i64_exact(2).unwrap();
        let (r, _) = one.agm(&two, RoundingMode::NearestEven);
        let one_again = FixedFloat::<113>::try_from_i64_exact(1).unwrap();
        let two_again = FixedFloat::<113>::try_from_i64_exact(2).unwrap();
        // AGM strictly between min and max.
        assert_eq!(r.partial_cmp(&one_again).0, Some(Ordering::Greater));
        assert_eq!(r.partial_cmp(&two_again).0, Some(Ordering::Less));
    }
}
