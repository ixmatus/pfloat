//! Property tests for `Ci` (slice 6m): the small-`x` envelope
//! `Ci(x) ≈ γ + ln x − x²/4` and the domain rules.

#![cfg(all(feature = "big", feature = "integrals"))]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Sign};
use proptest::prelude::*;

/// Euler–Mascheroni γ, a mathematical fact to ~40 digits (safe as a
/// literal; the in-repo value is pinned separately by slice 6m0).
const GAMMA: &str = "0.5772156649015328606065120900824024310421593359399235988058";

proptest! {
    /// Ci(x) = γ + ln x − x²/4 + O(x⁴) for small x > 0. Using
    /// x = 1/d with d ≥ 8, the O(x⁴) remainder is ≤ ~2⁻¹² of the
    /// leading γ term, so agreement to ~p−40 bits.
    #[test]
    fn ci_small_x_envelope(d in 8i64..=400) {
        let p = 160u32;
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let db = BigFloat::try_from_i64_exact(d, p).unwrap();
        let (x, _) = one.div(&db, RoundingMode::NearestEven);
        let (ci, _) = x.ci(RoundingMode::NearestEven);

        let gamma = BigFloat::parse_str(GAMMA, p, RoundingMode::NearestEven).unwrap().0;
        let (ln_x, _) = x.ln(RoundingMode::NearestEven);
        let (x_sq, _) = x.mul(&x, RoundingMode::NearestEven);
        let four = BigFloat::try_from_i64_exact(4, p).unwrap();
        let (x_sq_4, _) = x_sq.div(&four, RoundingMode::NearestEven);
        let (g_ln, _) = gamma.add(&ln_x, RoundingMode::NearestEven);
        let (approx, _) = g_ln.sub(&x_sq_4, RoundingMode::NearestEven);

        // Ci(x) − (γ + ln x − x²/4) = Σ_{k≥2} (−1)ᵏ x^{2k}/((2k)(2k)!),
        // an alternating decreasing tail, so its magnitude is ≤ the
        // first omitted term x⁴/96. Use x⁴/48 to absorb kernel
        // rounding while keeping the bound tight (a real O(x⁴)
        // property, not a hand-tuned bit count).
        let (diff, _) = ci.sub(&approx, RoundingMode::NearestEven);
        let err = diff.abs();
        let (x4, _) = x_sq.mul(&x_sq, RoundingMode::NearestEven);
        let f48 = BigFloat::try_from_i64_exact(48, p).unwrap();
        let (bound, _) = x4.div(&f48, RoundingMode::NearestEven);
        prop_assert!(
            matches!(
                err.partial_cmp(&bound).0,
                Some(Ordering::Less | Ordering::Equal)
            ),
            "Ci(1/{d})={ci} vs γ+ln x−x²/4={approx}, err={err}, bound x⁴/48={bound}"
        );
    }

    /// Ci(+0) = −∞ + DIV_BY_ZERO (pole).
    #[test]
    fn ci_zero_is_pole(p in prop_oneof![Just(53u32), Just(113u32)]) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, status) = z.ci(RoundingMode::NearestEven);
        prop_assert!(r.is_infinite() && r.is_sign_negative());
        prop_assert!(status.div_by_zero());
    }

    /// Ci(x < 0) = NaN + INVALID (complex in the reals).
    #[test]
    fn ci_negative_is_nan_invalid(num in 1i64..=1000) {
        let p = 64u32;
        let x = BigFloat::try_from_i64_exact(-num, p).unwrap();
        let (r, status) = x.ci(RoundingMode::NearestEven);
        prop_assert!(r.is_nan());
        prop_assert!(status.invalid());
    }
}
