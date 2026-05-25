//! Property tests for Airy `Ai` (slice 6n). Airy has no parity, so
//! the load-bearing identity is the Wronskian
//! `Ai·Bi′ − Ai′·Bi = 1/π` (DLMF 9.2.7); plus precision
//! self-consistency and the sign/monotonicity of `Ai` on `[0, ∞)`.
//!
//! Each `Ai` evaluation drives the gamma kernel twice, so the case
//! counts are kept small (these are correctness properties, not a
//! fuzz sweep).

#![cfg(all(feature = "big", feature = "airy"))]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Sign};
use proptest::prelude::*;

/// 1/π to 50 digits (mpmath; treated as a fact) — the crate exposes
/// no public π accessor.
const INV_PI: &str = "0.31830988618379067153776752674502872406891929148091";

fn close_within(a: &BigFloat, b: &BigFloat, bits: u32) -> bool {
    let (diff, _) = a.sub(b, RoundingMode::NearestEven);
    let abs_diff = diff.abs();
    if abs_diff.is_zero() {
        return true;
    }
    let p = a.precision().max(b.precision());
    let two = BigFloat::try_from_i64_exact(2, p).unwrap();
    let one = BigFloat::try_from_i64_exact(1, p).unwrap();
    let abs_b = b.abs();
    let mut bound = if abs_b.is_zero() { one } else { abs_b };
    for _ in 0..bits {
        bound = bound.div(&two, RoundingMode::NearestEven).0;
    }
    matches!(
        abs_diff.partial_cmp(&bound).0,
        Some(Ordering::Less | Ordering::Equal)
    )
}

fn rat(num: i64, den: i64, p: u32) -> BigFloat {
    let n = BigFloat::try_from_i64_exact(num, p).unwrap();
    let d = BigFloat::try_from_i64_exact(den, p).unwrap();
    n.div(&d, RoundingMode::NearestEven).0
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]

    /// Ai(0) is a finite positive normal value (≈ 0.355).
    #[test]
    fn ai_zero_is_finite_positive(p in prop_oneof![Just(53u32), Just(113u32)]) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (r, _) = z.ai(RoundingMode::NearestEven);
        prop_assert!(!r.is_nan() && !r.is_infinite() && !r.is_zero());
        prop_assert!(!r.is_sign_negative());
    }

    /// Wronskian: Ai·Bi′ − Ai′·Bi = 1/π at any x (DLMF 9.2.7).
    #[test]
    fn wronskian(num in -12i64..=12, den in 1i64..=4) {
        let p = 96u32;
        let x = rat(num, den, p);
        let (ai, _) = x.ai(RoundingMode::NearestEven);
        let (bi, _) = x.bi(RoundingMode::NearestEven);
        let (aip, _) = x.ai_prime(RoundingMode::NearestEven);
        let (bip, _) = x.bi_prime(RoundingMode::NearestEven);
        let (t1, _) = ai.mul(&bip, RoundingMode::NearestEven);
        let (t2, _) = aip.mul(&bi, RoundingMode::NearestEven);
        let (w, _) = t1.sub(&t2, RoundingMode::NearestEven);
        let inv_pi = BigFloat::parse_str(INV_PI, p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        prop_assert!(close_within(&w, &inv_pi, p - 16));
    }

    /// Precision self-consistency: Ai at p agrees with Ai at p+48
    /// rounded to p. The argument is an exact dyadic rational
    /// (denominator a power of two) so both precisions evaluate the
    /// *same* real point — a non-dyadic argument would differ in its
    /// low bits between the two precisions, and near a zero of Ai
    /// (the first at x ≈ −2.338) that argument mismatch is amplified
    /// by |Ai′/Ai| into a spurious failure even at the loose
    /// `p − 12` tolerance (the property of the test construction the
    /// pf-ok9 lesson identifies, also encoded in
    /// `property_yn::self_consistent`, `property_ik::self_consistent`,
    /// `property_zeta::self_consistent`, and
    /// `property_jn::self_consistent` per ADR-0036).
    #[test]
    fn ai_self_consistent(
        num in -10i64..=10,
        den in prop_oneof![Just(1i64), Just(2), Just(4)],
    ) {
        let p = 80u32;
        let x_lo = rat(num, den, p);
        let x_hi = rat(num, den, p + 48);
        let (lo, _) = x_lo.ai(RoundingMode::NearestEven);
        let (hi_raw, _) = x_hi.ai(RoundingMode::NearestEven);
        let hi = hi_raw
            .round_to_precision(p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        prop_assert!(close_within(&lo, &hi, p - 12));
    }

    /// On [0, ∞) Ai is strictly positive and strictly decreasing.
    #[test]
    fn ai_positive_decreasing(a in 0i64..=8, b in 9i64..=20) {
        let p = 80u32;
        let xa = BigFloat::try_from_i64_exact(a, p).unwrap();
        let xb = BigFloat::try_from_i64_exact(b, p).unwrap();
        let (fa, _) = xa.ai(RoundingMode::NearestEven);
        let (fb, _) = xb.ai(RoundingMode::NearestEven);
        prop_assert_eq!(fa.partial_cmp(&BigFloat::try_new_zero(Sign::Positive, p).unwrap()).0, Some(Ordering::Greater));
        prop_assert_eq!(fb.partial_cmp(&BigFloat::try_new_zero(Sign::Positive, p).unwrap()).0, Some(Ordering::Greater));
        prop_assert_eq!(fb.partial_cmp(&fa).0, Some(Ordering::Less));
    }
}
