//! `sin(x)`: trigonometric sine.
//!
//! Argument reduction (see [`super::trig_reduce`]): factor `x` into
//! a quadrant `q ∈ {0, 1, 2, 3}` and a reduced argument
//! `r ∈ [−π/4, π/4]`. Then:
//!
//! - `q = 0`: `sin(x) = sin(r)`
//! - `q = 1`: `sin(x) = cos(r)`
//! - `q = 2`: `sin(x) = −sin(r)`
//! - `q = 3`: `sin(x) = −cos(r)`
//!
//! `sin(r)` and `cos(r)` are evaluated via Taylor series at working
//! precision: `sin(r) = r − r³/3! + r⁵/5! − …` and
//! `cos(r) = 1 − r²/2! + r⁴/4! − …`. With `|r| ≤ π/4 ≈ 0.785`, the
//! series converge geometrically; termination is when a term falls
//! below `2^−working_prec` relative to the running sum.
//!
//! Correctly rounded under every IEEE 754-2019 rounding mode via
//! the shared [`crate::math::ziv::ziv_round`] driver (slice p1.26,
//! ADR-0038). The Payne-Hanek reduction + quadrant-dispatched
//! Taylor composition runs inside the eval closure at each Ziv
//! working precision; the range-cap NaN check pre-empts Ziv at
//! the maximum working precision the driver could request, so
//! every Ziv iteration's `reduce` call succeeds.
//!
//! Special cases per IEEE 754-2019 §9.2:
//!
//! - `sin(±0) = ±0`.
//! - `sin(±∞) = qNaN + INVALID`.
//! - `sin(NaN) = NaN`; `sNaN` raises `INVALID`.
//! - `|x|` past the reduction table budget (`~2^3000`):
//!   `qNaN + INVALID` (range exhausted).

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::trig_reduce::{reduce, reduction_depth_hint, Reduction};
use super::ziv::{ziv_round_with_depth, ZIV_BASE_GUARD};
use super::ziv_calibration::SIN_ERROR_GUARD;

impl BigFloat {
    /// `sin(self)` rounded under `mode` to `self.precision`.
    #[must_use]
    pub fn sin(&self, mode: RoundingMode) -> (Self, Status) {
        self.sin_round(self.precision, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `sin(self)` with explicit result precision.
    pub fn sin_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(sin_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `sin(self)` for `FixedFloat`. Delegates to [`BigFloat::sin`].
    #[must_use]
    pub fn sin(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().sin(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn sin_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
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
        Class::Infinity { .. } => {
            // sin(±∞) is undefined.
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        Class::Normal { .. } => {}
    }

    // Range-cap check at the Ziv first-iteration working precision.
    // The reduction's range policy caps the supported `|x|` at
    // `2^4032` (`reduce`'s `e_x + 64 < 4096`; since ADR-0103 the cap
    // is on the input's exponent alone — the old working-coupled
    // form refused deep working precisions outright, which the deep
    // certification rung legitimately requests). Pre-check at
    // `target + ZIV_BASE_GUARD` (the first iteration the driver
    // runs): if reduce fails here the input is fundamentally out of
    // range and no iteration could recover; a mid-loop None (now
    // unreachable by construction, kept as defense) makes the
    // closure return NaN and the post-Ziv check below raises
    // INVALID.
    // The pre-pf-1axr pre-check at `target + 1024` (the Ziv ceiling)
    // was over-conservative: it fired spuriously at
    // `target_precision ≥ 3008 − e_x`, blocking any caller above
    // that cliff even when the first iteration would have succeeded.
    // pf-1axr surfaced this via `bessel_y_eval_normal_at_w`'s
    // working-precision boost pushing `bessel_y_asymptotic`'s
    // internal `cos`/`sin` past the cliff at `Y2(1025)` at `p = 53`
    // (recurrence working = 3061).
    let ziv_first_working = target_precision.saturating_add(ZIV_BASE_GUARD);
    if reduce(x, ziv_first_working).is_none() {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }

    // Ziv-driven correct rounding under every IEEE mode. The eval
    // closure runs Payne-Hanek reduction + quadrant-dispatched
    // Taylor at working precision `w` under NE; the outer envelope
    // certifies the rounding-mode interval test. If reduce returns
    // None at a Ziv-grown `w` (only possible when the pre-check
    // passed but a doubled-guard working exceeds the table for this
    // `|x|`), return a NaN at the working precision; Ziv treats the
    // NaN as "interval can't be certified" and the post-Ziv check
    // below raises INVALID on the final result.
    let (result, status) = ziv_round_with_depth(
        |w| match reduce(x, w) {
            Some(Reduction { quadrant, r }) => match quadrant {
                0 => sin_taylor(&r, w),
                1 => cos_taylor(&r, w),
                2 => sin_taylor(&r, w).negated(),
                _ => cos_taylor(&r, w).negated(),
            },
            None => BigFloat::try_new_quiet_nan(Sign::Positive, w, &[]).expect("precision >= 1"),
        },
        target_precision,
        mode,
        SIN_ERROR_GUARD,
        // Inputs near a multiple of π/2 park the truth ~2|e_r| bits
        // from the grid (pf-jl35, ADR-0103; the cos-shape quadrants
        // 1 and 3 are sin's exposed arms); resolved lazily on
        // schedule exhaustion.
        || reduction_depth_hint(x, ziv_first_working),
    );
    // Post-Ziv: a Ziv iteration's reduce hitting the table cap
    // propagated NaN through the driver; surface as INVALID for
    // shape-parity with the explicit pre-check path above.
    if matches!(result.class, Class::Nan { .. }) && !status.invalid() {
        let merged = status.merge(Status::INVALID);
        auto_raise(Status::INVALID);
        return (result, merged);
    }
    // sin(x) for finite normal x is transcendental (Lindemann–
    // Weierstrass: sin α is transcendental for nonzero algebraic α,
    // and sin(x) = 0 only at x = nπ, never at a nonzero dyadic x), so
    // the rounded result is INEXACT even where it lands on a grid
    // value — e.g. a huge argument whose reduced residual collapses
    // (pf-njs5 under-report, ADR-0060). sin(±0) = ±0 is the only exact
    // input and is special-cased above.
    let status = if matches!(result.class, Class::Normal { .. }) {
        status | Status::INEXACT
    } else {
        status
    };
    auto_raise(status);
    (result, status)
}

/// `sin(r) = r − r³/3! + r⁵/5! − …` for `|r| ≤ π/4`.
pub(super) fn sin_taylor(r: &BigFloat, working_prec: u32) -> BigFloat {
    if r.is_zero() {
        return BigFloat::try_new_zero(r.sign(), working_prec).expect("precision >= 1");
    }

    let r_sq = r.mul(r, RoundingMode::NearestEven).0;
    let mut term = r.clone();
    let mut sum = term.clone();

    // Iterate: term_{n+1} = -term_n · r² / ((2n)(2n+1)). Starting
    // n = 1 with term_0 = r.
    let max_iter = working_prec.saturating_mul(2).max(256);
    for n in 1u32..=max_iter {
        let denom_val = i64::from(2 * n) * i64::from(2 * n + 1);
        let denom = BigFloat::try_from_i64_exact(denom_val, working_prec).expect("precision >= 1");
        let (numer, _) = term.mul(&r_sq, RoundingMode::NearestEven);
        let (next, _) = numer.div(&denom, RoundingMode::NearestEven);
        term = next.negated();
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

/// `cos(r) = 1 − r²/2! + r⁴/4! − …` for `|r| ≤ π/4`.
pub(super) fn cos_taylor(r: &BigFloat, working_prec: u32) -> BigFloat {
    let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    if r.is_zero() {
        return one;
    }

    let r_sq = r.mul(r, RoundingMode::NearestEven).0;
    let mut term = one.clone();
    let mut sum = one;

    let max_iter = working_prec.saturating_mul(2).max(256);
    for n in 1u32..=max_iter {
        let denom_val = i64::from(2 * n - 1) * i64::from(2 * n);
        let denom = BigFloat::try_from_i64_exact(denom_val, working_prec).expect("precision >= 1");
        let (numer, _) = term.mul(&r_sq, RoundingMode::NearestEven);
        let (next, _) = numer.div(&denom, RoundingMode::NearestEven);
        term = next.negated();
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
    fn sin_zero() {
        for s in [Sign::Positive, Sign::Negative] {
            let z = BigFloat::try_new_zero(s, 53).unwrap();
            let (r, _) = z.sin(RoundingMode::NearestEven);
            assert!(r.is_zero());
            assert_eq!(r.is_sign_negative(), matches!(s, Sign::Negative));
        }
    }

    #[test]
    fn sin_pos_inf_is_invalid() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, status) = pi.sin(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn sin_neg_inf_is_invalid() {
        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r, status) = ni.sin(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn sin_pi_is_zero() {
        let pi = super::super::pi_at(113);
        let (r, _) = pi.sin(RoundingMode::NearestEven);
        let zero = BigFloat::try_new_zero(Sign::Positive, 113).unwrap();
        // sin(π) is exactly 0 in real math; numerically we get a tiny
        // value near 2^(-target).
        assert!(close_at(&r, &zero, 100));
    }

    #[test]
    fn sin_pi_over_2_is_one() {
        let pi_2 = super::super::pi_over_2_at(113);
        let (r, _) = pi_2.sin(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        assert!(close_at(&r, &one, 113 - 12));
    }

    #[test]
    fn sin_neg_pi_over_2_is_neg_one() {
        let pi_2 = super::super::pi_over_2_at(113);
        let neg = pi_2.negated();
        let (r, _) = neg.sin(RoundingMode::NearestEven);
        let neg_one = BigFloat::try_from_i64_exact(-1, 113).unwrap();
        assert!(close_at(&r, &neg_one, 113 - 12));
    }

    #[test]
    fn sin_one() {
        // sin(1) ≈ 0.8414709848078965
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (r, _) = one.sin(RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "0.84147098480789650665250232163029900",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 100));
    }

    #[test]
    fn sin_is_odd() {
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let neg_two = BigFloat::try_from_i64_exact(-2, 113).unwrap();
        let (a, _) = two.sin(RoundingMode::NearestEven);
        let (b, _) = neg_two.sin(RoundingMode::NearestEven);
        let neg_a = a.negated();
        assert!(close_at(&b, &neg_a, 113 - 12));
    }

    #[test]
    fn sin_huge_x_returns_invalid() {
        // 2^5000 exceeds the reduction table budget.
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let two = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let mut big = one;
        for _ in 0..5000 {
            big = big.mul(&two, RoundingMode::NearestEven).0;
        }
        let (r, status) = big.sin(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn sin_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, _) = q.sin(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn sin_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.sin(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
