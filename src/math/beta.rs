//! `beta(a, b) = Γ(a) · Γ(b) / Γ(a + b)`: the Euler beta function.
//!
//! Implementation: route through `lgamma`:
//!
//! ```text
//! ln β(a, b) = lgamma(a) + lgamma(b) − lgamma(a + b)
//! β(a, b)    = sign · exp(ln β(a, b))
//! ```
//!
//! The combined sign is the product of the three Γ signs
//! (`gamma_sign_of`: `+` for a positive argument, `sign sin(πx)`
//! for a negative non-integer, derived from DLMF 5.5.3). For
//! positive `a, b` all three are positive and the sign is `+`.
//!
//! `beta` accepts the full real domain (ADR-0030, derived from
//! DLMF 5.12.1 / 5.5.3 / 5.2, every case pinned against mpmath).
//! With `Zle = {0, −1, −2, …}` the Γ poles:
//!
//! - `a, b, a+b` all off `Zle`: finite, magnitude via the `lgamma`
//!   composition, sign the product of the three Γ signs.
//! - `a+b ∈ Zle` with `a, b ∉ Zle`: `+0` (denominator pole).
//! - one operand a negative integer, the other a positive integer,
//!   `a+b ∈ Zle`: finite via pole cancellation,
//!   `(−1)^m (m−1)!(n−m)!/n!` evaluated through `lgamma` of the
//!   three positive-integer factorials (ADR-0030 case 4).
//! - one operand a negative integer otherwise: `qNaN + INVALID`
//!   (two-sided sign-ambiguous pole, mirrors `gamma`).
//! - one operand `±0` otherwise: `±∞ + DIV_BY_ZERO` (mirrors
//!   `gamma(±0)`).
//! - both operands in `Zle`: `qNaN + INVALID` (net pole).
//!
//! Special cases:
//!
//! - `beta(NaN, _) = beta(_, NaN) = NaN`; `sNaN` raises `INVALID`.
//! - `beta(+∞, b)` finite positive `b`: `+0`; symmetrically in `a`.
//! - any other infinite operand: `qNaN + INVALID` (the
//!   ∞-with-negative-operand edge is outside this extension; see
//!   ADR-0030 and DESIGN.md).

use super::gamma::gamma_sign_of;
use super::lgamma::is_integer_test;
use super::pow::{integer_parity, Parity};
use super::ziv::ziv_round;
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
    /// `beta(self, other)` rounded under `mode` to
    /// `max(self.precision, other.precision)`.
    #[must_use]
    pub fn beta(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let target = self.precision.max(other.precision);
        self.beta_round(other, target, mode)
            .expect("max of two valid precisions is valid")
    }

    /// `beta(self, other)` with explicit result precision.
    ///
    /// Correctly rounded under every IEEE 754-2019 rounding mode
    /// via the shared `crate::math::ziv::ziv_round` driver (slice
    /// p1.29, ADR-0038). The ADR-0030 case dispatch (poles,
    /// pole-cancellation, finite paths) stays before Ziv; only the
    /// finite-path lgamma compositions (cases 1/2 and the case-4
    /// closed form) run through the Ziv envelope.
    pub fn beta_round(
        &self,
        other: &Self,
        target_precision: u32,
        mode: RoundingMode,
    ) -> Result<(Self, Status), BuildError> {
        if target_precision == 0 {
            return Err(BuildError::PrecisionZero);
        }
        Ok(beta_kernel(self, other, target_precision, mode))
    }
}

#[cfg(feature = "fixed")]
impl<const PREC: u32> FixedFloat<PREC>
where
    [(); limbs_for(PREC)]:,
{
    /// `beta(self, other)` for `FixedFloat`. Delegates to
    /// [`BigFloat::beta`].
    #[must_use]
    pub fn beta(&self, other: &Self, mode: RoundingMode) -> (Self, Status) {
        let (big, status) = self.to_big().beta(&other.to_big(), mode);
        (
            Self::try_from_big_exact(big).expect("precision matches"),
            status,
        )
    }
}

fn beta_kernel(
    a: &BigFloat,
    b: &BigFloat,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    if a.is_signaling_nan() || b.is_signaling_nan() {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }
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

    // Infinity: beta(+∞, finite_positive) = +0.
    if (a.is_infinite() && matches!(a.sign(), Sign::Positive) && is_finite_positive(b))
        || (b.is_infinite() && matches!(b.sign(), Sign::Positive) && is_finite_positive(a))
    {
        return (
            BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1"),
            Status::OK,
        );
    }

    // Any infinite operand beyond the +∞ × finite-positive case
    // handled above is left at the conservative INVALID convention
    // (ADR-0030: the ∞-with-negative-operand edge is not part of the
    // negative-domain extension and is recorded in DESIGN.md).
    if a.is_infinite() || b.is_infinite() {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }

    // Both operands are finite. Classify against the Γ pole set
    // (the non-positive integers) per ADR-0030, derived from DLMF
    // 5.2 (poles only there, no zeros), 5.5.3 (reflection sign),
    // 5.12.1 (B = Γ(a)Γ(b)/Γ(a+b)). The decision table:
    //   6  a,b both non-positive integers      → qNaN+INVALID (pole)
    //   4  one a non-pos int, the other a pos
    //      int, a+b a non-pos int              → finite (cancellation)
    //   0  one operand ±0, not case 4          → ±∞+DIV_BY_ZERO
    //   3  one operand a negative integer,
    //      not case 4                          → qNaN+INVALID (pole)
    //   5  a+b a non-pos int, a,b not poles    → +0 (denominator pole)
    //   1/2 otherwise                          → signed lgamma path
    let a_npi = is_nonpos_integer(a);
    let b_npi = is_nonpos_integer(b);
    let (sum, _) = a.add(b, RoundingMode::NearestEven);
    let s_npi = is_nonpos_integer(&sum);

    // Case 6: both operands at Γ poles; the net is still a pole.
    if a_npi && b_npi {
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }

    // Case 4: pole cancellation to a finite value. One operand is a
    // negative integer −n, the other a positive integer m, and
    // a+b = m−n is a non-positive integer (so 1 ≤ m ≤ n). The closed
    // form (ADR-0030) is B(−n,m) = (−1)^m (m−1)!(n−m)!/n!; evaluated
    // through lgamma of the three positive-integer factorials (O(1),
    // no loop over the caller-supplied m).
    if a_npi && is_pos_integer(b) && s_npi {
        return beta_case4(a, b, target_precision, mode);
    }
    if b_npi && is_pos_integer(a) && s_npi {
        return beta_case4(b, a, target_precision, mode);
    }

    // Cases 0 and 3: exactly one operand sits at a Γ pole with no
    // compensating a+b pole. A ±0 operand mirrors gamma(±0) (signed
    // ∞ + DIV_BY_ZERO); a negative integer is a two-sided
    // sign-ambiguous pole (qNaN + INVALID), mirroring
    // gamma(negative integer).
    if a_npi || b_npi {
        let pole = if a_npi { a } else { b };
        if pole.is_zero() {
            let inf =
                BigFloat::try_new_infinity(pole.sign(), target_precision).expect("precision >= 1");
            auto_raise(Status::DIV_BY_ZERO);
            return (inf, Status::DIV_BY_ZERO);
        }
        let nan = BigFloat::try_new_quiet_nan(Sign::Positive, target_precision, &[])
            .expect("precision >= 1");
        auto_raise(Status::INVALID);
        return (nan, Status::INVALID);
    }

    // Case 5: a+b is a non-positive integer while a, b are not Γ
    // poles. Γ(a+b) is a denominator pole, 1/Γ(a+b) = 0, so B = +0
    // exactly (a genuine finite value, no exception).
    if s_npi {
        return (
            BigFloat::try_new_zero(Sign::Positive, target_precision).expect("precision >= 1"),
            Status::OK,
        );
    }

    // Cases 1 and 2: a, b, a+b all away from Γ poles. Magnitude via
    // the lgamma composition; sign is the product of the three Γ
    // signs (ADR-0030), reusing the single reflection derivation in
    // gamma_sign_of. The composition runs through ziv_round so the
    // directed-mode boundary lands correctly; the sign is binary
    // and applied after the magnitude is rounded.
    let (result, status) = ziv_round(
        |w| {
            let (lg_a, _) = a
                .lgamma_round(w, RoundingMode::NearestEven)
                .expect("precision >= 1");
            let (lg_b, _) = b
                .lgamma_round(w, RoundingMode::NearestEven)
                .expect("precision >= 1");
            let (lg_sum, _) = sum
                .lgamma_round(w, RoundingMode::NearestEven)
                .expect("precision >= 1");
            let (lg_a_plus_b, _) = lg_a.add(&lg_b, RoundingMode::NearestEven);
            let (lb, _) = lg_a_plus_b.sub(&lg_sum, RoundingMode::NearestEven);
            let (magnitude, _) = lb.exp(RoundingMode::NearestEven);
            let sign_a = gamma_sign_of(a, w);
            let sign_b = gamma_sign_of(b, w);
            let sign_sum = gamma_sign_of(&sum, w);
            if sign_product_is_negative(sign_a, sign_b, sign_sum) {
                magnitude.negated()
            } else {
                magnitude
            }
        },
        target_precision,
        mode,
    );
    auto_raise(status);
    (result, status)
}

/// `true` iff `x` is a non-positive integer, i.e. sits at a Γ pole
/// (DLMF 5.2: simple poles exactly at 0, −1, −2, …). Zero counts;
/// positive integers do not.
fn is_nonpos_integer(x: &BigFloat) -> bool {
    (x.is_zero() || matches!(x.sign(), Sign::Negative)) && is_integer_test(x)
}

/// `true` iff `x` is a strictly positive integer (1, 2, 3, …).
fn is_pos_integer(x: &BigFloat) -> bool {
    matches!(x.sign(), Sign::Positive) && !x.is_zero() && is_integer_test(x)
}

/// `true` iff the product of three Γ signs is negative.
fn sign_product_is_negative(a: Sign, b: Sign, c: Sign) -> bool {
    let neg = [a, b, c]
        .into_iter()
        .filter(|s| matches!(s, Sign::Negative))
        .count();
    neg % 2 == 1
}

/// Case 4 of ADR-0030: `B(neg, pos)` where `neg` is a negative
/// integer `−n`, `pos` a positive integer `m`, and `neg + pos` a
/// non-positive integer (so `1 ≤ m ≤ n`). The pole/pole cancellation
/// leaves the finite value
///
/// ```text
/// B(−n, m) = (−1)^m (m−1)! (n−m)! / n!,   n ≥ 1, 1 ≤ m ≤ n
/// ```
///
/// The three factorials are `Γ` of the *positive* integers `m`,
/// `n−m+1`, `n+1` (none a `Γ` pole), so the magnitude is
/// `exp(lgamma(m) + lgamma(n−m+1) − lgamma(n+1))` and the sign is
/// `(−1)^m`, read from the parity of `m`. This is `O(1)` (three
/// `lgamma` calls plus an `exp`); no huge binomial is formed and,
/// unlike the earlier reciprocal-product form, there is no loop over
/// the caller-supplied `m`. That form was exact but ran `m`
/// iterations — an unbounded cost on a caller-controlled integer
/// (a hang at large `m`); replacing it is the ADR-0030 robustness
/// fix. The `lgamma` composition trades the old form's exact
/// rational arithmetic for the same ~`p`-bit accuracy as the
/// negative-domain magnitude path (cases 1/2), the right call given
/// the alternative is non-termination.
fn beta_case4(
    neg: &BigFloat,
    pos: &BigFloat,
    target_precision: u32,
    mode: RoundingMode,
) -> (BigFloat, Status) {
    // Sign is (−1)^m, from the parity of m; binary and pinned
    // outside the Ziv envelope (`neg` and `pos` are exact integers
    // so the parity does not depend on the rounding mode or
    // working precision).
    let negative = matches!(integer_parity(pos), Some(Parity::Odd));
    let (result, status) = ziv_round(
        |w| {
            let ne = RoundingMode::NearestEven;
            let one = BigFloat::try_from_i64_exact(1, w).expect("precision >= 1");
            // n, m as exact integers at the working precision (neg,
            // pos are already exact integers, so widening cannot
            // round).
            let n = neg
                .negated()
                .round_to_precision(w, ne)
                .expect("precision >= 1")
                .0; // +n ≥ 1
            let m = pos.round_to_precision(w, ne).expect("precision >= 1").0; // m ≥ 1
            let (n_plus_1, _) = n.add(&one, ne); // n + 1 ≥ 2
            let (n_minus_m, _) = n.sub(&m, ne); // n − m ≥ 0
            let (n_minus_m_plus_1, _) = n_minus_m.add(&one, ne); // n−m+1 ≥ 1
            let (lg_m, _) = m.lgamma_round(w, ne).expect("precision >= 1"); // ln (m−1)!
            let (lg_nm1, _) = n_minus_m_plus_1
                .lgamma_round(w, ne)
                .expect("precision >= 1"); // ln (n−m)!
            let (lg_np1, _) = n_plus_1.lgamma_round(w, ne).expect("precision >= 1"); // ln n!
            let (sum, _) = lg_m.add(&lg_nm1, ne);
            let (ln_mag, _) = sum.sub(&lg_np1, ne);
            let (mag, _) = ln_mag.exp(ne);
            if negative {
                mag.negated()
            } else {
                mag
            }
        },
        target_precision,
        mode,
    );
    auto_raise(status);
    (result, status)
}

fn is_finite_positive(x: &BigFloat) -> bool {
    matches!(
        &x.class,
        Class::Normal {
            sign: Sign::Positive,
            ..
        }
    )
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
    fn beta_2_3_is_one_twelfth() {
        // β(2, 3) = 1/12.
        let a = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let b = BigFloat::try_from_i64_exact(3, 113).unwrap();
        let (r, _) = a.beta(&b, RoundingMode::NearestEven);
        let twelve = BigFloat::try_from_i64_exact(12, 113).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (expected, _) = one.div(&twelve, RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn beta_half_half_is_pi() {
        // β(1/2, 1/2) = π.
        let half = BigFloat::parse_str("0.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (r, _) = half.beta(&half, RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "3.1415926535897932384626433832795028841971693993751",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn beta_3_5() {
        // β(3, 5) = 1/105 ≈ 0.00952380952380952.
        let a = BigFloat::try_from_i64_exact(3, 113).unwrap();
        let b = BigFloat::try_from_i64_exact(5, 113).unwrap();
        let (r, _) = a.beta(&b, RoundingMode::NearestEven);
        let expected = BigFloat::parse_str(
            "0.0095238095238095238095238095238095238095238095238095",
            113,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, 80));
    }

    #[test]
    fn beta_is_symmetric() {
        // β(a, b) = β(b, a).
        let a = BigFloat::parse_str("2.5", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let b = BigFloat::parse_str("3.7", 113, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (ab, _) = a.beta(&b, RoundingMode::NearestEven);
        let (ba, _) = b.beta(&a, RoundingMode::NearestEven);
        assert!(close_at(&ab, &ba, 80));
    }

    fn p(s: &str) -> BigFloat {
        BigFloat::parse_str(s, 113, RoundingMode::NearestEven)
            .unwrap()
            .0
    }

    #[test]
    fn beta_negative_integer_no_cancellation_is_invalid() {
        // ADR-0030 case 3: a = −1 (negative integer), b = 2, a+b = 1
        // is positive so no pole cancellation → qNaN + INVALID.
        let a = BigFloat::try_from_i64_exact(-1, 53).unwrap();
        let b = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, status) = a.beta(&b, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
        // m > n (B(−3,5), mpmath +inf) is the same case.
        let a = BigFloat::try_from_i64_exact(-3, 53).unwrap();
        let b = BigFloat::try_from_i64_exact(5, 53).unwrap();
        let (r, status) = a.beta(&b, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan() && status.invalid());
    }

    #[test]
    fn beta_zero_operand_is_signed_pole() {
        // ADR-0030 case 0: a = +0, b = 2 → +∞ + DIV_BY_ZERO
        // (mirrors gamma(+0)); a = −0 → −∞ + DIV_BY_ZERO.
        let b = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let pz = BigFloat::try_new_zero(Sign::Positive, 53).unwrap();
        let (r, status) = pz.beta(&b, RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_positive());
        assert!(status.div_by_zero());
        let nz = BigFloat::try_new_zero(Sign::Negative, 53).unwrap();
        let (r, status) = nz.beta(&b, RoundingMode::NearestEven);
        assert!(r.is_infinite() && r.is_sign_negative());
        assert!(status.div_by_zero());
    }

    #[test]
    fn beta_case2_negative_non_integer_signed() {
        // ADR-0030 case 2, pinned vs mpmath (dps 55).
        let (r, _) = p("-0.5").beta(&p("0.25"), RoundingMode::NearestEven);
        assert!(close_at(
            &r,
            &p("2.622057554292119810464839589891119413682754951431623163"),
            80
        ));
        assert!(r.is_sign_positive());
        let (r, _) = p("-0.5").beta(&p("0.75"), RoundingMode::NearestEven);
        assert!(close_at(
            &r,
            &p("-1.198140234735592207439922492280323878227212663215651558"),
            80
        ));
        assert!(r.is_sign_negative());
        let (r, _) = p("-2.5").beta(&p("-1.25"), RoundingMode::NearestEven);
        assert!(close_at(
            &r,
            &p("-13.8385197111960899959311047858377407935243062601407755"),
            80
        ));
    }

    #[test]
    fn beta_case2_symmetric_on_negative_domain() {
        let (ab, _) = p("-0.5").beta(&p("0.75"), RoundingMode::NearestEven);
        let (ba, _) = p("0.75").beta(&p("-0.5"), RoundingMode::NearestEven);
        assert!(close_at(&ab, &ba, 80));
    }

    #[test]
    fn beta_case2_near_pole() {
        // a = −1 + 1/1024 (exact dyadic, non-integer), near the −1
        // Γ pole; finite, large magnitude. mpmath: −1025.0009775171…
        let a = p("-0.9990234375");
        let b = BigFloat::try_from_i64_exact(2, 113).unwrap();
        let (r, _) = a.beta(&b, RoundingMode::NearestEven);
        assert!(close_at(
            &r,
            &p("-1025.000977517106549364613880742913000977517106549364614"),
            70
        ));
    }

    #[test]
    fn beta_case5_sum_is_nonpos_integer_is_zero() {
        // ADR-0030 case 5: a,b not poles, a+b ∈ {0,−1,…} → +0.
        for (x, y) in [("-0.5", "0.5"), ("-1.5", "0.5"), ("-0.5", "-0.5")] {
            let (r, status) = p(x).beta(&p(y), RoundingMode::NearestEven);
            assert!(r.is_zero() && r.is_sign_positive(), "B({x},{y})");
            assert!(!status.invalid() && !status.div_by_zero());
        }
    }

    #[test]
    fn beta_case4_pole_cancellation_closed_form() {
        // ADR-0030 case 4: B(−n,m) = (−1)^m (m−1)!(n−m)!/n!,
        // pinned vs mpmath; symmetric.
        let cases: &[(i64, i64, &str)] = &[
            (
                -3,
                2,
                "0.16666666666666666666666666666666666666666666666667",
            ),
            (-1, 1, "-1.0"),
            (-5, 5, "-0.2"),
            (
                -5,
                3,
                "-0.03333333333333333333333333333333333333333333333333",
            ),
            (
                -3,
                1,
                "-0.33333333333333333333333333333333333333333333333333",
            ),
        ];
        for &(na, mb, exp) in cases {
            let a = BigFloat::try_from_i64_exact(na, 113).unwrap();
            let b = BigFloat::try_from_i64_exact(mb, 113).unwrap();
            let (r, status) = a.beta(&b, RoundingMode::NearestEven);
            assert!(!status.invalid(), "B({na},{mb}) should be finite");
            assert!(close_at(&r, &p(exp), 90), "B({na},{mb})");
            // symmetric: B(m,−n) == B(−n,m)
            let (rs, _) = b.beta(&a, RoundingMode::NearestEven);
            assert!(close_at(&rs, &p(exp), 90), "B({mb},{na}) symmetry");
        }
    }

    #[test]
    fn beta_case4_factorial_exact_rational() {
        // ADR-0030 case 4 cross-checked against a hand-derived exact
        // rational (no mpmath): B(−10, 4) = (−1)^4·3!·6!/10!
        // = (6·720)/3628800 = 4320/3628800 = 1/840.
        let a = BigFloat::try_from_i64_exact(-10, 113).unwrap();
        let b = BigFloat::try_from_i64_exact(4, 113).unwrap();
        let (r, status) = a.beta(&b, RoundingMode::NearestEven);
        assert!(!status.invalid() && !status.div_by_zero());
        assert!(r.is_sign_positive());
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let n840 = BigFloat::try_from_i64_exact(840, 113).unwrap();
        let (expected, _) = one.div(&n840, RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, 90));
        // symmetric
        let (rs, _) = b.beta(&a, RoundingMode::NearestEven);
        assert!(close_at(&rs, &expected, 90));
    }

    #[test]
    fn beta_case4_large_m_terminates() {
        // Robustness regression (ADR-0030): the earlier reciprocal-
        // product form looped `m` times, so this case — m = 10^12 —
        // did not terminate. The closed form is O(1); the test
        // returning at all is the property. With n = m+1,
        // B(−(m+1), m) = (−1)^m (m−1)!·1!/(m+1)! = (−1)^m /((m+1)·m);
        // m even ⇒ +1/(m·(m+1)).
        let m: i64 = 1_000_000_000_000;
        let neg = BigFloat::try_from_i64_exact(-(m + 1), 113).unwrap();
        let pos = BigFloat::try_from_i64_exact(m, 113).unwrap();
        let (r, status) = neg.beta(&pos, RoundingMode::NearestEven);
        assert!(!r.is_nan() && !r.is_infinite(), "got {r}");
        assert!(!status.invalid() && !status.div_by_zero());
        assert!(r.is_sign_positive(), "m even ⇒ positive");
        // Value cross-check: 1/(m·(m+1)).
        let m_big = BigFloat::try_from_i64_exact(m, 113).unwrap();
        let m_plus_1 = BigFloat::try_from_i64_exact(m + 1, 113).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 113).unwrap();
        let (denom, _) = m_big.mul(&m_plus_1, RoundingMode::NearestEven);
        let (expected, _) = one.div(&denom, RoundingMode::NearestEven);
        assert!(close_at(&r, &expected, 50), "B(−(m+1), m) = {r}");
    }

    #[test]
    fn beta_case6_both_nonpos_integers_is_invalid() {
        // ADR-0030 case 6: B(−2,−3) net pole (mpmath +inf) →
        // qNaN + INVALID.
        let a = BigFloat::try_from_i64_exact(-2, 53).unwrap();
        let b = BigFloat::try_from_i64_exact(-3, 53).unwrap();
        let (r, status) = a.beta(&b, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }

    #[test]
    fn beta_pos_inf_is_zero() {
        let pi = BigFloat::try_new_infinity(Sign::Positive, 53).unwrap();
        let b = BigFloat::try_from_i64_exact(2, 53).unwrap();
        let (r, _) = pi.beta(&b, RoundingMode::NearestEven);
        assert!(r.is_zero() && r.is_sign_positive());
    }

    #[test]
    fn beta_nan_propagates() {
        let q = BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, _) = q.beta(&one, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        let (r, _) = one.beta(&q, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
    }

    #[test]
    fn beta_snan_invalid() {
        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, 53, &[]).unwrap();
        let one = BigFloat::try_from_i64_exact(1, 53).unwrap();
        let (r, status) = sn.beta(&one, RoundingMode::NearestEven);
        assert!(r.is_quiet_nan());
        assert!(status.invalid());
    }
}
