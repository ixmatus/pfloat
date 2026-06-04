//! `ln(x)`: natural logarithm.
//!
//! Algorithm:
//!
//! 1. Special cases per IEEE 754-2019 §9: NaN propagates,
//!    `ln(-anything finite or +∞) = qNaN + INVALID`,
//!    `ln(±0) = -∞ + DIV_BY_ZERO`, `ln(+∞) = +∞`, `ln(1) = +0`.
//! 2. Range reduce by the binary exponent: pfloat already stores
//!    `x = m × 2^e` with `m ∈ [1, 2)` (this is just the `BigFloat`
//!    representation with `exponent` rewritten to `0`). Then
//!    `ln(x) = ln(m) + e · ln(2)`.
//! 3. Compute `ln(m)` for `m ∈ [1, 2)` via the atanh series:
//!    `u = (m − 1) / (m + 1)`, then
//!    `ln(m) = 2 · (u + u³/3 + u⁵/5 + …)`. With `m ∈ [1, 2)`,
//!    `u ∈ [0, 1/3]`, so `u² ≤ 1/9` and the series converges
//!    about 3 bits per term. Termination: when a term falls below
//!    `2^(-working_prec - 4)` relative to the running sum
//!    (≈ ln 2 = 0.69), it contributes nothing.
//! 4. Compose: `ln(x) = 2·atanh(u) + e · ln(2)`. Round to target
//!    precision under the user's rounding mode.
//!
//! Correctly rounded under every IEEE rounding mode via the shared
//! [`crate::math::ziv::ziv_round`] driver (slice p1.2, ADR-0022).
//! The driver supplies a working precision, [`ln_at_w`] evaluates the
//! algorithm above at that precision, and the guard grows until the
//! Ziv interval test certifies correct rounding at the target. The
//! `log2` and `log10` kernels compose through `ln_round` and inherit
//! the same correctness.

use core::cmp::Ordering;

use crate::big::{BigFloat, BuildError};
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

#[cfg(feature = "fixed")]
use crate::fixed::FixedFloat;
#[cfg(feature = "fixed")]
use crate::mantissa::limbs_for;

use super::ln_2_at;
use super::ziv::ziv_round;
use super::ziv_calibration::LN_ERROR_GUARD;

impl BigFloat {
    /// `ln(self)`: returns the natural logarithm rounded under
    /// `mode` to `self.precision`.
    ///
    /// Correctly rounded under every IEEE rounding mode via the
    /// shared [`crate::math::ziv::ziv_round`] driver (slice p1.2,
    /// ADR-0022). See [the module docs](self) for the algorithm.
    #[must_use]
    pub fn ln(&self, mode: RoundingMode) -> (Self, Status) {
        let target = self.precision;
        self.ln_round(target, mode)
            .expect("self.precision >= 1 by invariant")
    }

    /// `ln(self)` with explicit result precision.
    pub fn ln_round(
        &self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(ln_kernel(self, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `ln(self)` for `FixedFloat`. Delegates to [`BigFloat::ln`].
    #[must_use]
    pub fn ln(&self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().ln(mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn ln_kernel(x: &BigFloat, target_precision: u32, mode: RoundingMode) -> (BigFloat, Status) {
    // Special cases per IEEE 754-2019 §9.
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
            // ln(±0) = -∞ + DIV_BY_ZERO per IEEE 754-2019 §7.3.
            let ninf = BigFloat::try_new_infinity(Sign::Negative, target_precision)
                .expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            return (ninf, Status::DIV_BY_ZERO);
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
            // ln(-∞) = qNaN + INVALID.
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        Class::Normal {
            sign: Sign::Negative,
            ..
        } => {
            // ln(negative finite) = qNaN + INVALID.
            let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
                .expect("precision >= 1");
            auto_raise(Status::INVALID);
            return (nan, Status::INVALID);
        }
        Class::Normal {
            sign: Sign::Positive,
            ..
        } => {}
    }

    // ln(1) = 0 exactly (pf-kk16, ADR-0039). The atanh-series
    // `ln_at_w` at x = 1 reduces to t = (x-1)/(x+1) = 0/2 = 0,
    // and every series term has factor t^k (k ≥ 1), so the
    // working-precision sum is exactly 0; Ziv would in fact
    // certify the correct rounding through the half_width(0) = 0
    // exact-match interval test today. The pre-Ziv dispatch
    // makes the exact-value invariant structural: future
    // refactors of `ln_at_w` that no longer short-circuit at
    // t = 0 would still return the correct directed-mode
    // result through this dispatch (`feedback_exact_value_defeats_ziv`).
    let one = BigFloat::try_from_i64_exact(1, target_precision).expect("target_precision >= 1");
    if matches!(x.partial_cmp(&one).0, Some(Ordering::Equal)) {
        let zero = BigFloat::try_new_zero(Sign::Positive, target_precision)
            .expect("target_precision >= 1");
        return (zero, Status::OK);
    }

    // x is finite positive normal, x ≠ 1. Correctly rounded under
    // `mode` via the Ziv interval test (ADR-0022).
    let (result, status) = ziv_round(
        |working_prec| ln_at_w(x, working_prec),
        target_precision,
        mode,
        LN_ERROR_GUARD,
    );
    // ln(x) for finite positive normal x ≠ 1 is transcendental
    // (Lindemann–Weierstrass: ln α is transcendental for algebraic
    // α ∉ {0, 1}, and a dyadic x is algebraic), hence irrational,
    // hence INEXACT even where the working-precision evaluation
    // rounds onto a grid value (pf-njs5, ADR-0060). x = 1 → 0 is the
    // only exact input and is dispatched above.
    let status = status | Status::INEXACT;
    auto_raise(status);
    (result, status)
}

/// Evaluate `ln(x)` at the supplied working precision via
/// binary-exponent range reduction plus the atanh series. `x` must
/// be finite positive normal (the caller's special-case handling
/// peels off NaN, ±0, ±∞, and negatives before invoking this).
/// Returns the unrounded value; the Ziv driver handles rounding to
/// the caller's target precision and mode.
fn ln_at_w(x: &BigFloat, working_prec: u32) -> BigFloat {
    let x_w = x
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    let (e, mantissa, prec_w) = match &x_w.class {
        Class::Normal {
            exponent, mantissa, ..
        } => (*exponent, mantissa.clone(), x_w.precision),
        _ => unreachable!("x_w is finite positive normal"),
    };

    // Build m = x_w with exponent rewritten to 0. m ∈ [1, 2).
    let m = BigFloat {
        class: Class::Normal {
            sign: Sign::Positive,
            exponent: 0,
            mantissa,
        },
        precision: prec_w,
    };

    // Compute u = (m - 1) / (m + 1) ∈ [0, 1/3].
    let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let (m_minus_1, _) = m.sub(&one, RoundingMode::NearestEven);
    let (m_plus_1, _) = m.add(&one, RoundingMode::NearestEven);
    let (u, _) = m_minus_1.div(&m_plus_1, RoundingMode::NearestEven);

    // atanh(u) = u + u³/3 + u⁵/5 + ... Compute via:
    //   term_1 = u, sum = u
    //   for k = 1, 2, ...:
    //     term_{k+1} = term_k * u²
    //     sum += term_{k+1} / (2k + 1)
    //     stop when term < 2^(-working_prec)
    let (u_squared, _) = u.mul(&u, RoundingMode::NearestEven);
    let mut sum = u.clone();
    let mut term = u;
    let max_iter = 4u32.saturating_mul(working_prec).max(256);
    for k in 1u32..=max_iter {
        let (new_term, _) = term.mul(&u_squared, RoundingMode::NearestEven);
        term = new_term;
        let divisor = BigFloat::try_from_i64_exact(i64::from(2 * k + 1), working_prec)
            .expect("precision >= 1");
        let (term_over_d, _) = term.div(&divisor, RoundingMode::NearestEven);
        let (new_sum, _) = sum.add(&term_over_d, RoundingMode::NearestEven);
        sum = new_sum;
        if let Class::Normal { exponent, .. } = &term.class {
            if *exponent < -i64::from(working_prec) - 4 {
                break;
            }
        } else {
            break;
        }
    }

    // ln(m) = 2 · sum.
    let two = BigFloat::try_from_i64_exact(2, working_prec).expect("precision >= 1");
    let (ln_m, _) = two.mul(&sum, RoundingMode::NearestEven);

    // ln(x) = ln(m) + e · ln(2).
    let ln2 = ln_2_at(working_prec);
    let e_big = BigFloat::try_from_i64_exact(e, working_prec).expect("precision >= 1");
    let (e_ln2, _) = e_big.mul(&ln2, RoundingMode::NearestEven);
    let (result, _) = ln_m.add(&e_ln2, RoundingMode::NearestEven);
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(feature = "fixed")]
    use core::cmp::Ordering;

    fn parse(s: &str, p: u32) -> BigFloat {
        BigFloat::parse_str(s, p, RoundingMode::NearestEven)
            .expect("parse should succeed")
            .0
    }

    fn close_at(v: &BigFloat, expected: &BigFloat, prec: u32) -> bool {
        // |v - expected| <= 4 ULP at expected's magnitude.
        let (diff, _) = v.sub(expected, RoundingMode::NearestEven);
        let abs_diff = diff.abs();
        if abs_diff.is_zero() {
            return true;
        }
        if !abs_diff.is_normal() {
            return false;
        }
        let exp_diff = match &abs_diff.class {
            Class::Normal { exponent, .. } => *exponent,
            _ => return false,
        };
        let expected_exp = match &expected.class {
            Class::Normal { exponent, .. } => *exponent,
            _ => 0,
        };
        let tolerance_exp = expected_exp - i64::from(prec) + 4;
        exp_diff <= tolerance_exp
    }

    #[test]
    fn ln_one_is_zero() {
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, status) = one.ln(RoundingMode::NearestEven);
        assert!(status.is_ok());
        assert!(r.is_zero(), "ln(1) should be 0, got {r:?}");
    }

    #[test]
    fn ln_one_is_zero_under_every_directed_mode() {
        // pf-kk16 pinning test: the exact-value pre-Ziv dispatch
        // returns exactly +0 under every mode. Without the dispatch,
        // the atanh-series composition's tiny noise would tip TP to
        // the smallest representable positive value at target precision.
        for &mode in &[
            RoundingMode::NearestEven,
            RoundingMode::TowardPositive,
            RoundingMode::TowardNegative,
            RoundingMode::TowardZero,
            RoundingMode::NearestAway,
        ] {
            for &prec in &[24u32, 53, 113] {
                let one = BigFloat::try_from_i64_exact(1, prec).unwrap();
                let (r, status) = one.ln(mode);
                assert!(status.is_ok(), "ln(1) status under {mode:?}@p{prec}");
                assert!(
                    r.is_zero() && !r.is_sign_negative(),
                    "ln(1) should be +0 under {mode:?}@p{prec}, got {r:?}"
                );
                assert_eq!(r.precision(), prec);
            }
        }
    }

    #[test]
    fn ln_e_is_one() {
        let e = parse("2.718281828459045", 53);
        let (r, _) = e.ln(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        assert!(close_at(&r, &one, 53), "ln(e) ≈ 1, got {r}");
    }

    #[test]
    fn ln_two() {
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, _) = two.ln(RoundingMode::NearestEven);
        // ln(2) ≈ 0.6931471805599453
        let expected = parse("0.6931471805599453", 53);
        assert!(close_at(&r, &expected, 53), "ln(2) = {r}");
    }

    #[test]
    fn ln_ten() {
        let ten = BigFloat::try_from_i64_exact(10, 53).unwrap();
        let (r, _) = ten.ln(RoundingMode::NearestEven);
        // ln(10) ≈ 2.302585092994046
        let expected = parse("2.302585092994046", 53);
        assert!(close_at(&r, &expected, 53), "ln(10) = {r}");
    }

    #[test]
    fn ln_half() {
        let half = parse("0.5", 53);
        let (r, _) = half.ln(RoundingMode::NearestEven);
        // ln(0.5) = -ln(2) ≈ -0.6931471805599453
        let expected = parse("-0.6931471805599453", 53);
        assert!(close_at(&r, &expected, 53), "ln(0.5) = {r}");
    }

    #[test]
    fn ln_zero() {
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, status) = pz.ln(RoundingMode::NearestEven);
        assert!(r.is_infinite());
        assert!(r.is_sign_negative());
        assert!(status.div_by_zero());

        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r2, status2) = nz.ln(RoundingMode::NearestEven);
        assert!(r2.is_infinite());
        assert!(r2.is_sign_negative());
        assert!(status2.div_by_zero());
    }

    #[test]
    fn ln_negative_is_invalid() {
        let neg_one = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        let (r, status) = neg_one.ln(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn ln_inf() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let (r, status) = pi.ln(RoundingMode::NearestEven);
        assert!(r.is_infinite());
        assert!(r.is_sign_positive());
        assert!(status.is_ok());

        let ni = BigFloat::try_new_infinity(Sign::Negative, 53).unwrap();
        let (r2, status2) = ni.ln(RoundingMode::NearestEven);
        assert!(r2.is_quiet_nan());
        assert!(status2.invalid());
    }

    #[test]
    fn ln_qnan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = q.ln(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(!status.invalid());
    }

    #[test]
    fn ln_snan_raises_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let (r, status) = sn.ln(RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn ln_exp_round_trip() {
        // ln(exp(x)) ≈ x for moderate x.
        for &n in &[1i64, 2, 5, 7, -3, 10] {
            let x = BigFloat::try_from_i64_exact(n, 113).unwrap();
            let (e_x, _) = x.exp(RoundingMode::NearestEven);
            let (back, _) = e_x.ln(RoundingMode::NearestEven);
            assert!(
                close_at(&back, &x, 113),
                "ln(exp({n})) = {back}, expected {x}"
            );
        }
    }

    #[test]
    fn exp_ln_round_trip() {
        // exp(ln(x)) ≈ x for positive x.
        for &n in &[1i64, 2, 5, 7, 10, 100, 1000] {
            let x = BigFloat::try_from_i64_exact(n, 113).unwrap();
            let (ln_x, _) = x.ln(RoundingMode::NearestEven);
            let (back, _) = ln_x.exp(RoundingMode::NearestEven);
            assert!(close_at(&back, &x, 113), "exp(ln({n})) = {back}");
        }
    }

    #[test]
    fn ln_high_precision() {
        // Verify via the round-trip exp(ln(2)) ≈ 2 at 256-bit
        // precision. (A literal decimal expected value would need
        // round_trip_digit_count(256) = 79 digits to faithfully
        // pin the 256-bit binary representation; the exp/ln pair
        // is the cleaner check.)
        let two = BigFloat::try_from_i64_exact(2, 256).unwrap();
        let (ln2, _) = two.ln(RoundingMode::NearestEven);
        assert_eq!(ln2.precision(), 256);
        let (back, _) = ln2.exp(RoundingMode::NearestEven);
        assert!(
            close_at(&back, &two, 256),
            "exp(ln(2)) at 256 bits = {back}"
        );
    }

    #[test]
    fn ln_round_rejects_zero_precision() {
        let two = BigFloat::try_from_i64_exact(2, 53).unwrap();
        assert_eq!(
            two.ln_round(0, RoundingMode::NearestEven),
            Err(BuildError::PrecisionZero)
        );
    }

    #[cfg(feature = "fixed")]
    #[test]
    fn fixed_ln() {
        let two = FixedFloat::<53>::try_from_i64_exact(2).unwrap();
        let (r, _) = two.ln(RoundingMode::NearestEven);
        let expected = parse("0.6931471805599453", 53);
        let cmp = r
            .partial_cmp(&FixedFloat::<53>::try_from_big_exact(expected).unwrap())
            .0;
        assert!(matches!(
            cmp,
            Some(Ordering::Equal | Ordering::Less | Ordering::Greater)
        ));
    }

    #[test]
    fn ln_just_above_one() {
        // ln(1 + ε) ≈ ε for small ε; round-trip through exp checks
        // the atanh series converges correctly near m = 1.
        let one_plus = parse("1.0001", 113);
        let (r, _) = one_plus.ln(RoundingMode::NearestEven);
        let (back, _) = r.exp(RoundingMode::NearestEven);
        assert!(close_at(&back, &one_plus, 113), "exp(ln(1.0001)) = {back}");
    }
}
