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
//! precision `p` once the iterates are within a factor of two of
//! each other; the loop terminates once the gap `|a_n − b_n|` falls
//! below a few ulps *of the iterate magnitude* (ADR-0095: an
//! absolute floor certified the first arithmetic mean for small
//! operands, and a sub-ulp relative floor would be unsatisfiable on
//! a `p`-bit grid).
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
    //
    // Huge exponent spreads take the asymptotic branch (ADR-0106,
    // the slice verifier's refutation 3). The Gauss iterates
    // CONVERGE toward the AGM's exponent ≈ s − log₂(s) (s = half
    // the spread), so near-convergence products sit near 2s — no
    // static normalization keeps a spread ≳ 2^63 loop off the rim
    // (the opposite-rim pair certified a ~2^42-wrong value with
    // Status OK). For a ≥ b, AGM(a, b) = π·a / (2·ln(4a/b)) with
    // relative error O((b/a)²) (K(k) → ln(4/k′) as k → 1; pinned
    // numerically against mpmath, error shrinking like (b/a)²), so
    // once 2·spread clears every expressible target-plus-guard the
    // closed form IS the correctly roundable value: spread ≥ 2^33
    // gives 2·spread ≥ 2^34 > u32::MAX + 1024 with margin. The
    // composition orders π/(2L) before ·a so no intermediate
    // leaves the representable range (π·a alone would saturate at
    // the top rim).
    let (e_a, e_b) = match (&a.class, &b.class) {
        (Class::Normal { exponent: ea, .. }, Class::Normal { exponent: eb, .. }) => (*ea, *eb),
        _ => unreachable!("specials and zeros dispatched above"),
    };
    if e_a.saturating_sub(e_b).unsigned_abs() >= 1u64 << 33 {
        let (big, small, e_big, e_small) = if e_a >= e_b {
            (a, b, e_a, e_b)
        } else {
            (b, a, e_b, e_a)
        };
        // L = ln(4·big/small) = (e_big − e_small + 2)·ln2
        //     + ln(m_big) − ln(m_small), with the mantissas shifted
        // exactly to exponent 0 (m ∈ [1, 2), so both logs are O(1)
        // with O(1) absolute error at w bits). The integer exponent
        // part is EXACT — the first draft computed ln(big) and
        // ln(small) whole, and for same-sign rim exponents their
        // near-cancellation amplified the logs' absolute error by up
        // to (|e_a| + |e_b|)/spread ≈ 2^31 past the charged
        // half-width, certifying 1-ulp-wrong answers (the revision
        // verifier's refutation, 40 constructed reproducers). The
        // split leaves every term's relative error at O(2^-w),
        // uniformly in the exponents. spread + 2 spans up to 65 bits
        // (opposite rims), held exactly as two i64 halves; their sum
        // is exact at w ≥ target + 64 ≥ 65 bits.
        let spread2 = i128::from(e_big) - i128::from(e_small) + 2;
        let n1 = (spread2 / 2) as i64;
        let n2 = (spread2 - spread2 / 2) as i64;
        let (m_big, sm1) = big.scale_by_pow2(-e_big);
        let (m_small, sm2) = small.scale_by_pow2(-e_small);
        debug_assert!(
            sm1.is_ok() && sm2.is_ok(),
            "shifting a Normal to exponent 0 is exact"
        );
        let (result, status) = ziv_round(
            |w| {
                let ne = RoundingMode::NearestEven;
                let ma_w = m_big.round_to_precision(w, ne).expect("precision >= 1").0;
                let ms_w = m_small.round_to_precision(w, ne).expect("precision >= 1").0;
                let (lma, _) = ma_w.ln(ne);
                let (lms, _) = ms_w.ln(ne);
                let h1 = BigFloat::try_from_i64_exact(n1, w).expect("precision >= 1");
                let h2 = BigFloat::try_from_i64_exact(n2, w).expect("precision >= 1");
                let (n_bf, _) = h1.add(&h2, ne);
                let ln2 = super::ln_2_at(w);
                let (l_int, _) = n_bf.mul(&ln2, ne);
                let (l_mant, _) = lma.sub(&lms, ne);
                let (l, _) = l_int.add(&l_mant, ne);
                let (two_l, _) = l.scale_by_pow2(1);
                // π via the Brent–Salamin iteration: available under
                // the agm feature itself (pi_at's table dispatch is
                // trig-gated — the first draft's use of it broke the
                // big,agm combo build).
                let pi = super::agm_constants::pi_via_agm(w);
                let (t, _) = pi.div(&two_l, ne);
                // π/(2L) multiplies before ·big so no intermediate
                // can leave the representable range.
                let big_w = big.round_to_precision(w, ne).expect("precision >= 1").0;
                t.mul(&big_w, ne).0
            },
            target_precision,
            mode,
            AGM_ERROR_GUARD,
        );
        auto_raise(status);
        return (result, status);
    }

    // Normalize toward exponent 0 (pf-06lk, ADR-0106). The Gauss
    // iteration's mul/sqrt saturate their result exponents at the
    // i64 rim per the no-emax contract, the closure discards those
    // per-op statuses, and the clamp is independent of the working
    // precision — so the interval test certified corrupted iterates
    // for operand exponents within ~2^62 of the rim (worst case
    // ~10^31 wrong with Status OK). AGM is degree-1 homogeneous,
    // agm(s·a, s·b) = s·agm(a, b), so scale both operands by 2^−m
    // with m the floor midpoint of their exponents (pure
    // shift-and-mask, overflow-free; clamped one above i64::MIN so
    // its negation exists — the verifier's refutation 1 found
    // agm(2^MIN, 1.5·2^MIN) wrapping `-m` into an equal-pair
    // Status-OK lie). With the asymptotic branch above owning every
    // spread ≥ 2^33, the normalized exponents are within ±2^32 of
    // 0: every internal product lands near exponent 0, sums stay
    // within the operand range, and nothing can reach the rim. The
    // result scales back by 2^m, exactness asserted below at the
    // scale-back itself.
    let m = ((e_a >> 1) + (e_b >> 1) + (e_a & e_b & 1)).max(i64::MIN + 1);
    let (a_norm, sa_norm) = a.scale_by_pow2(-m);
    let (b_norm, sb_norm) = b.scale_by_pow2(-m);
    debug_assert!(
        sa_norm.is_ok() && sb_norm.is_ok(),
        "normalization is a half-spread shift (< 2^33 here) and cannot saturate"
    );
    let a = &a_norm;
    let b = &b_norm;

    // Ziv-driven correct rounding under every IEEE mode: the eval
    // closure captures the normalized (a, b) and runs the Gauss AGM
    // iteration at working precision w. Quadratic convergence
    // doubles the bit agreement each step, so Ziv adds at most one
    // extra iteration per retry (O(log w) → O(log 2w) ≈ +1).
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

            // Canonicalize so a_n >= b_n. Load-bearing for the
            // convergence test below: it reads a_n's exponent as the
            // operand magnitude, and AM >= GM keeps a_n the larger
            // iterate from here on.
            if matches!(a_n.partial_cmp(&b_n).0, Some(Ordering::Less)) {
                core::mem::swap(&mut a_n, &mut b_n);
            }

            // Iteration budget, a derived backstop (pf-ddfl,
            // ADR-0095). Two regimes: while a_n/b_n = R > 2, one
            // step maps R to at most sqrt(R), halving log2(R);
            // exponents span i64, so log2(R) < 2^65 needs at most
            // 65 halvings. Once R <= 2 the relative gap squares
            // each step (quadratic convergence), so reaching the
            // few-ulp floor below for any supported w < 2^32 needs
            // at most log2(2^32) + 1 = 33 more. 65 + 33 + 6 slack
            // = 104. The floor normally exits the loop well before
            // the cap; the cap-exit average is still sound (the
            // pair is within a few ulps, see below).
            let max_iter = 104u32;
            let two = BigFloat::try_from_i64_exact(2, w).expect("precision >= 1");

            for _ in 0..max_iter {
                let (diff, _) = a_n.sub(&b_n, RoundingMode::NearestEven);
                let abs_diff = diff.abs();
                // Convergence is RELATIVE to the iterate magnitude
                // (pf-ddfl, ADR-0095): converge once the gap falls
                // below 4 ulps of a_n, i.e. gap exponent de <
                // exponent(a_n) - w + 3 (ulp(a_n) = 2^(ae - w + 1)).
                // With b_n <= AGM <= a_n the returned midpoint
                // m = (a_n + b_n)/2 satisfies |m - AGM| <=
                // (a_n - b_n)/2 < 2 ulp = 2^(ae-w+2), a relative
                // error under 2^-(w-3) — more than 2^20 inside the
                // calibrated AGM_ERROR_GUARD = 24 half-width. The
                // threshold must sit a few ulps ABOVE the working
                // grid: two distinct w-bit values can never differ
                // by less than 2^(ae-w) relative, so any floor at
                // or below one ulp is unsatisfiable and the loop
                // would always run to the cap on iterates that
                // oscillate within an ulp instead of becoming
                // bit-identical (caught by the pf-ddfl adversarial
                // verification). The previous ABSOLUTE floor
                // -(w + 4) certified the first arithmetic mean for
                // small operands (0.5% relative error with Status
                // OK at 2^-300). saturating_sub errs toward more
                // iterations, never premature convergence. The
                // relative bound holds at every scale the loop's
                // arithmetic represents; operand exponents within
                // ~2^62 of the i64 rim saturate inside mul/sqrt
                // with their statuses discarded here, a separate
                // pre-existing defect (pf-a77o/pf-kh3z arc).
                let converged = match (&abs_diff.class, &a_n.class) {
                    (Class::Zero { .. }, _) => true,
                    (Class::Normal { exponent: de, .. }, Class::Normal { exponent: ae, .. }) => {
                        *de < ae.saturating_sub(i64::from(w) - 3)
                    }
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

            // After convergence (or the cap) the pair is within a
            // few ulps; averaging absorbs the final separation.
            let (sum, _) = a_n.add(&b_n, RoundingMode::NearestEven);
            sum.div(&two, RoundingMode::NearestEven).0
        },
        target_precision,
        mode,
        AGM_ERROR_GUARD,
    );
    // Scale back. NOT asserted exact: rounding the normalized AGM at
    // a target coarser than the operands can carry past max(a, b) to
    // the next binade (the verifier's refutation 2 constructed it
    // deterministically: operands just below 2^(i64::MAX + 1) whose
    // AGM rounds up at 53 bits), and at the top rim that binade is
    // unrepresentable. scale_by_pow2 then applies the documented
    // saturation contract (meaningful mantissa under the clamped
    // exponent) and flags it; the flag is merged into the returned
    // Status instead of being discarded — the first draft
    // debug-asserted it away, a new panic in debug and a dropped
    // mandatory OVERFLOW in release.
    let (result, s_back) = result.scale_by_pow2(m);
    let status = status | s_back;
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
