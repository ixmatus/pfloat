//! `expm1(x) = exp(x) − 1`: exponential minus one.
//!
//! Naively `exp(x) − 1` loses precision near zero: for `x ≈ 2^−n`,
//! `exp(x) ≈ 1 + x`, so the subtraction cancels the leading bits of
//! the sum and the result has ~`n` bits of cancellation. The
//! kernel boosts the working precision by the cancellation amount
//! before calling `exp`, then subtracts.
//!
//! For `|x|` so small that `expm1(x)` rounds to `x` at the target
//! precision (specifically, `x.exponent ≤ −target − 8`, so the
//! `x²/2` term is below half a ULP of `x`), the kernel short-
//! circuits and returns `x` directly without invoking `exp` at all.
//!
//! Correctly rounded under every IEEE 754-2019 rounding mode via
//! the shared [`crate::math::ziv::ziv_round`] driver (slice p1.24,
//! ADR-0038). The cancellation-boost composition runs inside the
//! eval closure at each Ziv working precision; the outer envelope
//! certifies the rounding-mode interval test on the final round.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `expm1(±0) = ±0` (sign preserved).
//! - `expm1(+∞) = +∞`, `expm1(−∞) = −1`.
//! - `expm1(NaN) = NaN`; `sNaN` raises `INVALID`.

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

use super::ziv::ziv_round_with_depth;
use super::ziv_calibration::EXPM1_ERROR_GUARD;

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

impl BigFloat {
    /// `expm1(self) = exp(self) − 1` rounded under `mode` to
    /// `self.precision`.
    #[must_use]
    pub fn expm1(&self, mode: RoundingMode) -> (Self, Status) {
        self.expm1_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `expm1(self)` with explicit result precision.
    pub fn expm1_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(expm1_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `expm1(self)` for `FixedFloat`. Delegates to
    /// [`BigFloat::expm1`].
    #[must_use]
    pub fn expm1(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().expm1(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn expm1_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
            let neg_one =
                BigFloat::try_from_i64_exact(-1, target_precision).expect("precision >= 1");
            return (neg_one, Status::OK);
        }
        Class::Normal { .. } => {}
    }

    let e = match &x.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => unreachable!(),
    };

    // Tiny x: expm1(x) = x + x²/2 + … lies within a sub-ULP of x, so
    // round x with that positive infinitesimal directly. The Ziv path
    // below caps its cancellation boost at +1024 bits, so for very tiny
    // x the exp(x)−1 composition collapses to exactly 0 and the interval
    // test certifies the false 0 (half_width(0)=0); review 2026-05-29.
    // The depth must also clear the INPUT's grid (pf-fbjn, ADR-0104):
    // a high-precision x parked next to a rounding-change point puts
    // the quadratic series correction (position 2e) across a boundary
    // the residue (position e − p − 2) never reaches. Arm-failing
    // inputs go to the driver, whose ADR-0103 deep rung takes the
    // input at full precision and certifies the true boundary side.
    if e <= -(i64::from(target_precision) + 2) && e <= -(i64::from(x.precision) + 3) {
        return crate::rounding::round_with_infinitesimal(
            x,
            x.sign(),
            x.is_sign_negative(),
            target_precision,
            mode,
        );
    }

    // Large negative x: expm1(x) = −1 + e^x approaches −1 from above
    // (e^x → 0⁺), magnitude strictly < 1. Past the Ziv guard cap the
    // residual e^x underflows every working precision, the interval test
    // never converges, and the fallback returns the exact on-grid −1 —
    // correct under nearest and TowardNegative, but wrong under TowardZero
    // and TowardPositive (the modes rounding toward zero). Short-circuit
    // to the mode-aware rounding of −1 with a magnitude-shrinking
    // infinitesimal. (x → +∞ grows without bound; no saturation there.)
    if x.is_sign_negative() && e >= super::saturation_threshold_exponent(target_precision) {
        let one = BigFloat::try_from_i64_exact(1, target_precision).expect("precision >= 1");
        return crate::rounding::round_with_infinitesimal(
            &one,
            Sign::Negative,
            true,
            target_precision,
            mode,
        );
    }

    // Positive exponent-rim forwarding (pf-qm0h family, ADR-0101).
    // The Ziv closure below discards exp's Status, so exp's mode-aware
    // rim dispatch (ADR-0096) arrived as a bare +inf that
    // half_width(non-Normal) = 0 certified — and the Class::Normal
    // INEXACT force then skipped it, emitting Status::OK on a
    // transcendental result. Past the rim the −1 is unobservable at
    // ANY expressible target: e^x ≥ 2^(2^62) puts the ulp at
    // 2^(2^62 − target − 1) ≫ 1 for every target < 2^32, so
    // expm1(x) and exp(x) round identically under every mode, flags
    // included. Forward verbatim; exp's triage owns the rim.
    if !x.is_sign_negative() && e >= 62 {
        return x
            .exp_round(target_precision, mode)
            .expect("target_precision >= 1 (validated by expm1_round)");
    }

    // Cancellation boost: e^x − 1 loses ~|exponent(x)| leading
    // bits when x is small. The boost moves INSIDE the Ziv eval
    // closure so each working-precision retry inherits it.
    //
    // The pre-Phase-1f kernel had a short-circuit at
    // `e ≤ -target - 8` that returned `x rounded under mode`.
    // That shortcut was NE-correct (when x²/2 < ULP/2, expm1(x)
    // rounds to x under NE) but produced wrong directed-mode
    // results: for positive small x under TP, true expm1(x) > x
    // is strictly above x, so TP must round up to the next-up
    // f32 neighbour, but the short-circuit returned x. Pf-l6s5
    // precedent applies: trust the Ziv driver to certify the
    // correct rounding mode behaviour. The Ziv loop converges
    // in 1-2 iterations on tiny-x inputs.
    let cancellation: u32 = if e < 0 {
        u32::try_from(-e).unwrap_or(u32::MAX)
    } else {
        0
    };

    // Ziv-driven correct rounding under every IEEE mode. The eval
    // closure runs the existing composition (`exp(x_w) − 1`) at a
    // working precision boosted by `cancellation` above the Ziv
    // driver's requested working precision `w`, then rounds the
    // composition's result to `w` under NE so the Ziv interval
    // test sees a w-precision value with the cancellation absorbed.
    let (result, status) = ziv_round_with_depth(
        |w| {
            // The +1024 cap on the internal boost predates the
            // ADR-0059 fast path; with the ADR-0104 precision arm,
            // every driver-reached tiny input has |e| ≤ max(p + 2,
            // target + 2) (deeper inputs take the fast path), so the
            // uncapped boost is input-proportional — and the cap was
            // actively harmful in the arm-fail zone: it let the
            // composition collapse to exactly 0, which half_width(0)
            // = 0 certifies instantly with Status OK (the
            // review-2026-05-29 certified-zero class, resurrected).
            let inner_w = w.saturating_add(cancellation);
            let x_w = x
                .round_to_precision(inner_w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0;
            let (e_x, _) = x_w.exp(RoundingMode::NearestEven);
            let one = BigFloat::try_from_i64_exact(1, inner_w).expect("precision >= 1");
            let (diff, _) = e_x.sub(&one, RoundingMode::NearestEven);
            diff.round_to_precision(w, RoundingMode::NearestEven)
                .expect("precision >= 1")
                .0
        },
        target_precision,
        mode,
        EXPM1_ERROR_GUARD,
        // Parked-input certification depth (pf-fbjn, ADR-0104):
        // arm-rejected tiny inputs resolve at the deep rung, which
        // must reach both the input's precision and the series
        // correction's depth. Lazy: free unless the schedule exhausts.
        || {
            if e < 0 {
                u32::try_from(e.saturating_mul(-3))
                    .unwrap_or(u32::MAX)
                    .max(x.precision)
                    .saturating_add(64)
            } else {
                0
            }
        },
    );
    // expm1(x) = e^x − 1 for finite normal x ≠ 0 is transcendental
    // (e^x is transcendental for nonzero algebraic x by Lindemann–
    // Weierstrass, and subtracting 1 keeps it irrational), hence INEXACT
    // even where it rounds onto a grid value (pf-uqd1, ADR-0063).
    // expm1(±0) = ±0 and expm1(−∞) = −1 are dispatched above; the tiny-x
    // fast path sets INEXACT via round_with_infinitesimal.
    let status = super::force_transcendental_inexact(&result, status);
    auto_raise(status);
    (result, status)
}

#[cfg(test)]
mod tests {
    use super::*;

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
        use core::cmp::Ordering;
        matches!(
            abs_diff.partial_cmp(&bound).0,
            Some(Ordering::Less | Ordering::Equal)
        )
    }

    #[test]
    fn expm1_zero_is_zero() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.expm1(RoundingMode::NearestEven);
            assert!(r.is_zero());
            assert_eq!(r.is_sign_negative(), matches!(s, Sign::Negative));
        }
    }

    #[test]
    fn expm1_pos_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, _) = pi.expm1(RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
    }

    #[test]
    fn expm1_neg_inf_is_neg_one() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, _) = ni.expm1(RoundingMode::NearestEven);
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        use core::cmp::Ordering;
        assert_eq!(r.partial_cmp(&neg_one).0, Some(Ordering::Equal));
    }

    #[test]
    fn expm1_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.expm1(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn expm1_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.expm1(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn expm1_small_x_no_cancellation() {
        // The naïve `exp(x) − 1` at target precision loses ~|exp(x)|
        // worth of bits to cancellation against 1. expm1 boosts
        // working precision so the result keeps the leading bits of
        // `x`. Test: for x = 2^−50, expm1(x) agrees with x to about
        // 49 bits (relative error ~x/2 from the x²/2 second-order
        // term).
        let p = 113u32;
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let mut x = one;
        for _ in 0..50 {
            x = x.div(&two, RoundingMode::NearestEven).0;
        }
        let (r, _) = x.expm1(RoundingMode::NearestEven);
        assert!(close_at(&r, &x, 49));
    }

    #[test]
    fn expm1_at_one() {
        // expm1(1) = e - 1 ≈ 1.71828...
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.expm1(RoundingMode::NearestEven);
        let (e_one, _) = one.exp(RoundingMode::NearestEven);
        let (expected, _) = e_one.sub(&one, RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, 113 - 12));
    }

    #[test]
    fn expm1_large() {
        // expm1(20) = exp(20) − 1. No cancellation issue: the −1
        // is negligible relative to exp(20) ≈ 4.85e8.
        let twenty = BigFloat::try_from_i64_exact(20, 113).unwrap();
        let (r, _) = twenty.expm1(RoundingMode::NearestEven);
        let (ex, _) = twenty.exp(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (expected, _) = ex.sub(&one, RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, 113 - 12));
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixed_expm1() {
        let zero = FixedFloat::<53>::zero();
        let (r, _) = zero.expm1(RoundingMode::NearestEven);
        assert!(r.is_zero());
    }
}
