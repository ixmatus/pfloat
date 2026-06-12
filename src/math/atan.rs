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
use super::ziv_calibration::ATAN_ERROR_GUARD;
use super::{pi_over_2_at, pi_over_2_at_round, signed_constant_at_round};

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
            // atan(±∞) = ±π/2, mode-aware. The negative case rounds π/2
            // under the mirrored mode before negating (Phase 4
            // directed-mode constant audit; atan(−∞, TowardNegative) used
            // to land above −π/2).
            let (signed, status) =
                signed_constant_at_round(pi_over_2_at_round, *sign, target_precision, mode);
            crate::status::auto_raise(status);
            return (signed, status);
        }
        Class::Normal { .. } => {}
    }

    // Tiny x: atan(x) = x − x³/3 + … shrinks toward zero (every
    // post-x term opposes the leading term in magnitude), and past
    // the representable band the correction can never reach any
    // working grid the Ziv driver visits: the eval collapses onto
    // the on-grid argument, the interval test never converges, and
    // the exhausted fall-through returned the argument itself — 1
    // ULP wrong under the inward modes (pf-e2ow's root, ADR-0102).
    // The correction c = |x| − |atan x| lies in (2^(3e−2), 2^(3e+2));
    // round_with_infinitesimal is exact when both c and its residue
    // stay strictly inside the boundary-free zone below |x|, whose
    // width is at least 2^(e − max(p, target) − 2) (one input ulp
    // when x sits off the target's rounding-change grid, half a
    // change step when on it, both with binade-crossing slack), so
    // `2|e| ≥ max(p, target) + 6` clears it with two bits to spare.
    // Unlike the ADR-0059 fast paths this trigger carries the
    // input-precision arm — the boundary-free zone shrinks with the
    // INPUT's grid, not only the target's.
    let e = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => unreachable!("specials dispatched above"),
    };
    // The rim guard mirrors hypot's (ADR-0102 verifier finding): within
    // max(p, target) + 5 of i64::MIN the infinitesimal's residue
    // placement saturates and the dispatch would certify a wrong value
    // with Status OK; refuse it there (pre-existing rim behavior, fixed
    // at the root by pf-a77o).
    let max_pt = i64::from(x.precision.max(target_precision));
    if e < 0
        && e.saturating_mul(-2) >= max_pt.saturating_add(6)
        && e >= i64::MIN.saturating_add(max_pt).saturating_add(5)
    {
        return crate::rounding::round_with_infinitesimal(
            x,
            x.sign(),
            true, // magnitude shrinks: the −x³/3 correction opposes x's sign
            target_precision,
            mode,
        );
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
        ATAN_ERROR_GUARD,
    );
    // atan(x) for finite normal x ≠ 0 is transcendental (Lindemann–
    // Weierstrass), hence irrational, hence INEXACT even where it rounds
    // onto a grid value (pf-uqd1, ADR-0063). atan(±0) = ±0 is dispatched
    // above; atan(±∞) = ±π/2 already carries INEXACT from
    // pi_over_2_at_round (the irrational constant rounds inexactly).
    let status = if matches!(result.class, Class::Normal { .. }) {
        status | Status::INEXACT
    } else {
        status
    };
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

    #[test]
    fn atan_inf_directed_rounding_is_sound() {
        // Regression (Phase 4 directed-mode constant audit): atan(−∞) used
        // to round on the wrong side of −π/2 (constant rounded under
        // `mode` then negated without mirroring TN↔TP).
        let hp = crate::math::pi_over_2_at(600);
        let nhp = hp.negated();
        for &p in &[24u32, 53, 113, 200] {
            let pi = BigFloat::try_new_infinity(Sign::Positive, p).unwrap();
            let ni = BigFloat::try_new_infinity(Sign::Negative, p).unwrap();
            assert_ne!(
                pi.atan(RoundingMode::TowardNegative).0.partial_cmp(&hp).0,
                Some(Ordering::Greater),
                "atan(+inf, TN) ≤ π/2 at p={p}"
            );
            assert_ne!(
                pi.atan(RoundingMode::TowardPositive).0.partial_cmp(&hp).0,
                Some(Ordering::Less),
                "atan(+inf, TP) ≥ π/2 at p={p}"
            );
            assert_ne!(
                ni.atan(RoundingMode::TowardNegative).0.partial_cmp(&nhp).0,
                Some(Ordering::Greater),
                "atan(-inf, TN) ≤ −π/2 at p={p}"
            );
            assert_ne!(
                ni.atan(RoundingMode::TowardPositive).0.partial_cmp(&nhp).0,
                Some(Ordering::Less),
                "atan(-inf, TP) ≥ −π/2 at p={p}"
            );
        }
    }

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
