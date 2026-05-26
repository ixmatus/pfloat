//! `atan(x)`: arc tangent, returning a value in `[−π/2, π/2]`.
//!
//! Algorithm: range-reduce by the identity
//! `atan(x) = π/2 − atan(1/x)` for `|x| > 1`, then apply the half-
//! angle identity `atan(y) = 2 · atan(y / (1 + sqrt(1 + y²)))` a
//! handful of times to bring `|y|` below `1/16`. The Taylor series
//! `atan(y) = y − y³/3 + y⁵/5 − …` then converges at roughly four
//! bits per term.
//!
//! Correctly rounded under every IEEE 754-2019 rounding mode via
//! the shared [`crate::math::ziv::ziv_round`] driver (slice p1.25,
//! ADR-0038). The `atan(±∞) = ±π/2` special case rounds via
//! [`super::pi_over_2_at_round`].
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `atan(±0) = ±0`.
//! - `atan(±∞) = ±π/2`.
//! - `atan(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::ziv::ziv_round;
use super::{pi_over_2_at, pi_over_2_at_round};

impl BigFloat {
    /// `atan(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn atan(&self, mode: RoundingMode) -> (Self, Status) {
        self.atan_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `atan(self)` with explicit result precision.
    pub fn atan_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(atan_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `atan(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::atan`].
    #[must_use]
    pub fn atan(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().atan(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn atan_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            let z = BigFloat::try_new_zero(*sign, target_precision).expect("precision >= 1");
            return (z, Status::OK);
        }
        Class::Infinity { sign } => {
            // atan(±∞) = ±π/2. Mode-aware (slice p1.25).
            let (pi_2, status) = pi_over_2_at_round(target_precision, mode);
            let signed = if matches!(sign, Sign::Negative) {
                pi_2.negated()
            } else {
                pi_2
            };
            crate::status::auto_raise(status);
            return (signed, status);
        }
        Class::Normal { .. } => {}
    }

    // Ziv-driven correct rounding under every IEEE mode. The eval
    // closure runs the existing half-angle reduction + Taylor
    // composition (atan_finite_unsigned on |x|) at working precision
    // `w`; the sign-flip for negative x happens inside eval so the
    // returned value's class matches the kernel's domain.
    let is_negative = matches!(x.sign(), Sign::Negative);
    let abs_x = x.abs();
    let (result, status) = ziv_round(
        |w| {
            let result = atan_finite_unsigned(&abs_x, w);
            if is_negative {
                result.negated()
            } else {
                result
            }
        },
        target_precision,
        mode,
    );
    auto_raise(status);
    (result, status)
}

/// `atan(|x|)` for finite normal positive `x`. Returns a value in
/// `[0, π/2]` at `working_prec`.
pub(super) fn atan_finite_unsigned(abs_x: &BigFloat, working_prec: u32) -> BigFloat {
    let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let two = BigFloat::try_from_i64_exact(2, working_prec).expect("precision >= 1");

    // For |x| > 1, atan(x) = π/2 − atan(1/x).
    let (y_initial, subtract_from_pi_half) = match abs_x.partial_cmp(&one).0 {
        Some(core::cmp::Ordering::Greater) => {
            let (recip, _) = one.div(abs_x, RoundingMode::NearestEven);
            (recip, true)
        }
        _ => (abs_x.clone(), false),
    };

    // Half-angle reduction: y ← y / (1 + sqrt(1 + y²)). Each step
    // shrinks |y| by roughly a factor of two. Stop when |y| < 1/16
    // so the Taylor series converges at ~4 bits per term.
    let mut y = y_initial
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;
    let mut k: u32 = 0;
    while should_halve(&y) {
        let (y_sq, _) = y.mul(&y, RoundingMode::NearestEven);
        let (one_plus_sq, _) = one.add(&y_sq, RoundingMode::NearestEven);
        let (s, _) = one_plus_sq.sqrt(RoundingMode::NearestEven);
        let (denom, _) = one.add(&s, RoundingMode::NearestEven);
        let (next, _) = y.div(&denom, RoundingMode::NearestEven);
        y = next;
        k += 1;
        if k >= 64 {
            break;
        }
    }

    let mut sum = atan_taylor(&y, working_prec);

    // Reverse the half-angles: multiply by 2^k.
    for _ in 0..k {
        let (doubled, _) = sum.mul(&two, RoundingMode::NearestEven);
        sum = doubled;
    }

    if subtract_from_pi_half {
        let pi_2 = pi_over_2_at(working_prec);
        let (diff, _) = pi_2.sub(&sum, RoundingMode::NearestEven);
        diff
    } else {
        sum
    }
}

/// Returns `true` if `y` is large enough that another half-angle
/// step would meaningfully accelerate Taylor convergence. Threshold
/// `|y| ≥ 1/16` (exponent ≥ −4).
fn should_halve(y: &BigFloat) -> bool {
    match &y.class {
        Class::Normal { exponent, .. } => *exponent >= -4,
        _ => false,
    }
}

/// `atan(y) = y − y³/3 + y⁵/5 − …` for `|y| < 1`. Best convergence
/// at small `|y|`.
fn atan_taylor(y: &BigFloat, working_prec: u32) -> BigFloat {
    if y.is_zero() {
        return BigFloat::try_new_zero(y.sign(), working_prec).expect("precision >= 1");
    }

    let (y_sq, _) = y.mul(y, RoundingMode::NearestEven);
    let mut x_power = y.clone();
    let mut sum = y.clone();
    let mut alternating_sign = true; // next term is subtracted

    let max_iter = working_prec.saturating_mul(2).max(256);
    for n in 1u32..=max_iter {
        let (next_power, _) = x_power.mul(&y_sq, RoundingMode::NearestEven);
        x_power = next_power;
        let denom_val = i64::from(2 * n + 1);
        let denom = BigFloat::try_from_i64_exact(denom_val, working_prec).expect("precision >= 1");
        let (mut term, _) = x_power.div(&denom, RoundingMode::NearestEven);
        if alternating_sign {
            term = term.negated();
        }
        alternating_sign = !alternating_sign;
        let (next_sum, _) = sum.add(&term, RoundingMode::NearestEven);
        sum = next_sum;

        if let Class::Normal { exponent, .. } = &term.class {
            if *exponent < -i64::from(working_prec) - 4 {
                break;
            }
        } else {
            break;
        }
    }

    sum
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
    fn atan_zero() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.atan(RoundingMode::NearestEven);
            assert!(r.is_zero());
            assert_eq!(r.is_sign_negative(), matches!(s, Sign::Negative));
        }
    }

    #[test]
    fn atan_pos_inf_is_pi_over_2() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 113).unwrap();
        let (r, _) = pi.atan(RoundingMode::NearestEven);
        let pi_2 = super::super::pi_over_2_at(113);
        assert!(close_at(&r, &pi_2, 100));
    }

    #[test]
    fn atan_neg_inf_is_neg_pi_over_2() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 113).unwrap();
        let (r, _) = ni.atan(RoundingMode::NearestEven);
        let pi_2 = super::super::pi_over_2_at(113);
        let neg = pi_2.negated();
        assert!(close_at(&r, &neg, 100));
    }

    #[test]
    fn atan_one_is_pi_over_4() {
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.atan(RoundingMode::NearestEven);
        let pi_2 = super::super::pi_over_2_at(113);
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (pi_4, _) = pi_2.div(&two, RoundingMode::NearestEven);
        assert!(close_at(&r, &pi_4, 113 - 12));
    }

    #[test]
    fn atan_neg_one_is_neg_pi_over_4() {
        let neg_one = BigFloat::try_from_i64_exact(-1, 113).unwrap();
        let (r, _) = neg_one.atan(RoundingMode::NearestEven);
        let pi_2 = super::super::pi_over_2_at(113);
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (pi_4, _) = pi_2.div(&two, RoundingMode::NearestEven);
        let neg_pi_4 = pi_4.negated();
        assert!(close_at(&r, &neg_pi_4, 113 - 12));
    }

    #[test]
    fn atan_is_odd() {
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let neg_two = BigFloat::try_from_i64_exact(-2, 113).unwrap();
        let (a, _) = two.atan(RoundingMode::NearestEven);
        let (b, _) = neg_two.atan(RoundingMode::NearestEven);
        let neg_a = a.negated();
        assert!(close_at(&b, &neg_a, 113 - 12));
    }

    #[test]
    fn atan_tan_round_trip() {
        // atan(tan(x)) = x for x ∈ (−π/2, π/2).
        let p = 113u32;
        for n in &[-1i64, 0, 1] {
            let x = BigFloat::try_from_i64_exact(*n, p).unwrap();
            let (t, _) = x.tan(RoundingMode::NearestEven);
            let (back, _) = t.atan(RoundingMode::NearestEven);
            assert!(close_at(&back, &x, p - 16), "atan(tan({n})) = {back}");
        }
    }

    #[test]
    fn atan_large_argument() {
        // atan(1000) lies just below π/2. Round-trip through tan
        // gives back the original argument; that's a tighter test
        // than the π/2 − 1/x first-order approximation.
        let p = 113u32;
        let large = BigFloat::try_from_i64_exact(1000, p).unwrap();
        let (r, _) = large.atan(RoundingMode::NearestEven);
        let (back, _) = r.tan(RoundingMode::NearestEven);
        assert!(close_at(&back, &large, p - 24));
    }

    #[test]
    fn atan_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.atan(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn atan_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.atan(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
