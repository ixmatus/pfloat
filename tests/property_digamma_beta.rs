//! Property-based tests for the second half of Phase 4's
//! tier-1 specials shipped in slice 4c: `digamma` and `beta`.

#![cfg(all(feature = "big", feature = "specials"))]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode};
use proptest::prelude::*;

fn arb_precision() -> impl Strategy<Value = u32> {
    prop_oneof![Just(53u32), Just(64u32), Just(113u32)]
}

/// `num/den` as an exact pfloat dyadic (`den` a power of two) at
/// precision `p`.
fn dyadic(num: i64, den: i64, p: u32) -> BigFloat {
    let n = BigFloat::try_from_i64_exact(num, p).unwrap();
    let d = BigFloat::try_from_i64_exact(den, p).unwrap();
    n.div(&d, RoundingMode::NearestEven).0
}

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

proptest! {
    // ---- digamma ----

    /// digamma(x + 1) = digamma(x) + 1/x. Recurrence identity.
    #[test]
    fn digamma_recurrence(n in 1i64..=20) {
        let p = 113u32;
        let x = BigFloat::try_from_i64_exact(n, p).unwrap();
        let x_plus_1 = BigFloat::try_from_i64_exact(n + 1, p).unwrap();
        let (psi_x, _) = x.digamma(RoundingMode::NearestEven);
        let (psi_x_plus_1, _) = x_plus_1.digamma(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (recip, _) = one.div(&x, RoundingMode::NearestEven);
        let (lhs, _) = psi_x.add(&recip, RoundingMode::NearestEven);
        prop_assert!(
            close_within(&lhs, &psi_x_plus_1, p.saturating_sub(20)),
            "ψ({}) + 1/{} = {}, ψ({}+1) = {}",
            n, n, lhs, n, psi_x_plus_1,
        );
    }

    /// digamma is monotonic on positive integers (strictly
    /// increasing).
    #[test]
    fn digamma_monotonic(n in 2i64..=30) {
        let p = 113u32;
        let a = BigFloat::try_from_i64_exact(n, p).unwrap();
        let b = BigFloat::try_from_i64_exact(n + 1, p).unwrap();
        let (pa, _) = a.digamma(RoundingMode::NearestEven);
        let (pb, _) = b.digamma(RoundingMode::NearestEven);
        prop_assert_eq!(pb.partial_cmp(&pa).0, Some(Ordering::Greater));
    }

    /// digamma(1/2) reflection: ψ(1/2) = -γ - 2·ln(2).
    /// Equivalent: ψ(1) - ψ(1/2) = 2·ln(2).
    #[test]
    fn digamma_half_relates_to_ln_two(_x in 0u32..1) {
        let p = 113u32;
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let half = BigFloat::parse_str("0.5", p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (psi_one, _) = one.digamma(RoundingMode::NearestEven);
        let (psi_half, _) = half.digamma(RoundingMode::NearestEven);
        let (diff, _) = psi_one.sub(&psi_half, RoundingMode::NearestEven);
        // Should equal 2·ln(2).
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let (ln_2, _) = two.ln(RoundingMode::NearestEven);
        let (expected, _) = two.mul(&ln_2, RoundingMode::NearestEven);
        prop_assert!(close_within(&diff, &expected, p.saturating_sub(20)));
    }

    // ---- beta ----

    /// β(a, b) = β(b, a). Symmetry.
    #[test]
    fn beta_is_symmetric(a in 1i64..=10, b in 1i64..=10) {
        let p = 113u32;
        let av = BigFloat::try_from_i64_exact(a, p).unwrap();
        let bv = BigFloat::try_from_i64_exact(b, p).unwrap();
        let (ab, _) = av.beta(&bv, RoundingMode::NearestEven);
        let (ba, _) = bv.beta(&av, RoundingMode::NearestEven);
        prop_assert!(close_within(&ab, &ba, p.saturating_sub(16)));
    }

    /// For positive integers a, b: β(a, b) = (a-1)! · (b-1)! / (a + b - 1)!.
    #[test]
    fn beta_integer_factorial(a in 1i64..=8, b in 1i64..=8) {
        let p = 113u32;
        let av = BigFloat::try_from_i64_exact(a, p).unwrap();
        let bv = BigFloat::try_from_i64_exact(b, p).unwrap();
        let (r, _) = av.beta(&bv, RoundingMode::NearestEven);

        // Build the expected factorial ratio.
        let mut numer = BigFloat::try_from_i64_exact(1, p).unwrap();
        for k in 1..a {
            let k_big = BigFloat::try_from_i64_exact(k, p).unwrap();
            numer = numer.mul(&k_big, RoundingMode::NearestEven).0;
        }
        for k in 1..b {
            let k_big = BigFloat::try_from_i64_exact(k, p).unwrap();
            numer = numer.mul(&k_big, RoundingMode::NearestEven).0;
        }
        let mut denom = BigFloat::try_from_i64_exact(1, p).unwrap();
        for k in 1..(a + b) {
            let k_big = BigFloat::try_from_i64_exact(k, p).unwrap();
            denom = denom.mul(&k_big, RoundingMode::NearestEven).0;
        }
        let (expected, _) = numer.div(&denom, RoundingMode::NearestEven);
        prop_assert!(
            close_within(&r, &expected, p.saturating_sub(20)),
            "β({a}, {b}) = {r}, expected = {expected}",
        );
    }

    /// β(1, b) = 1/b.
    #[test]
    fn beta_one_b_is_reciprocal(b in 1i64..=20) {
        let p = 113u32;
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let bv = BigFloat::try_from_i64_exact(b, p).unwrap();
        let (r, _) = one.beta(&bv, RoundingMode::NearestEven);
        let (expected, _) = one.div(&bv, RoundingMode::NearestEven);
        prop_assert!(close_within(&r, &expected, p.saturating_sub(16)));
    }

    /// β(a, b) > 0 for a, b > 0.
    #[test]
    fn beta_positive(a in 1i64..=10, b in 1i64..=10) {
        let p = 113u32;
        let av = BigFloat::try_from_i64_exact(a, p).unwrap();
        let bv = BigFloat::try_from_i64_exact(b, p).unwrap();
        let (r, _) = av.beta(&bv, RoundingMode::NearestEven);
        prop_assert!(r.is_sign_positive());
        prop_assert!(!r.is_zero());
    }

    /// ADR-0030 row 3: a *negative integer* with no a+b pole
    /// cancellation is a two-sided sign-ambiguous Γ pole, so β
    /// returns qNaN/INVALID at any precision. (β is no longer
    /// invalid for negative *non-integers* — see the properties
    /// below; only the integer-pole subcase keeps this behavior.)
    #[test]
    fn beta_non_positive_invalid(p in arb_precision()) {
        let neg = BigFloat::try_from_i64_exact(-1, p).unwrap();
        let pos = BigFloat::try_from_i64_exact(2, p).unwrap();
        let (r, status) = neg.beta(&pos, RoundingMode::NearestEven);
        prop_assert!(r.is_quiet_nan());
        prop_assert!(status.invalid());
    }

    /// β(a, b) = β(b, a) on the negative domain (ADR-0030 case 2):
    /// a = −(2i+1)/2 (negative non-integer), b = ±(2j+1)/4, so a+b
    /// has an odd numerator over 4 and is never a pole.
    #[test]
    fn beta_symmetric_negative_non_integer(i in 1i64..=8, j in 0i64..=8, bneg in any::<bool>()) {
        let p = 113u32;
        let a = dyadic(-(2 * i + 1), 2, p);
        let bn = if bneg { -(2 * j + 1) } else { 2 * j + 1 };
        let b = dyadic(bn, 4, p);
        let (ab, sab) = a.beta(&b, RoundingMode::NearestEven);
        let (ba, sba) = b.beta(&a, RoundingMode::NearestEven);
        prop_assert!(ab.is_finite() && !ab.is_nan());
        prop_assert!(!sab.invalid() && !sba.invalid());
        prop_assert!(
            close_within(&ab, &ba, p.saturating_sub(16)),
            "β({a}, {b}) = {ab} vs β({b}, {a}) = {ba}",
        );
    }

    /// β(a, b) = Γ(a)·Γ(b)/Γ(a+b) on the negative domain (DLMF
    /// 5.12.1, the defining Γ-quotient continuation; ADR-0030
    /// case 2). An in-crate reflection-consistency check that needs
    /// no external oracle and carries the sign through Γ. Modest
    /// ranges keep the Γ magnitudes tame.
    #[test]
    fn beta_equals_gamma_quotient_negative(i in 1i64..=4, j in 0i64..=4, bneg in any::<bool>()) {
        let p = 113u32;
        let a = dyadic(-(2 * i + 1), 2, p);
        let bn = if bneg { -(2 * j + 1) } else { 2 * j + 1 };
        let b = dyadic(bn, 4, p);
        let (sum, _) = a.add(&b, RoundingMode::NearestEven);
        let (r, status) = a.beta(&b, RoundingMode::NearestEven);
        prop_assert!(r.is_finite() && !status.invalid());
        let (ga, _) = a.gamma(RoundingMode::NearestEven);
        let (gb, _) = b.gamma(RoundingMode::NearestEven);
        let (gab, _) = sum.gamma(RoundingMode::NearestEven);
        let (num, _) = ga.mul(&gb, RoundingMode::NearestEven);
        let (expected, _) = num.div(&gab, RoundingMode::NearestEven);
        prop_assert!(
            close_within(&r, &expected, p.saturating_sub(30)),
            "β({a}, {b}) = {r}, Γ-quotient = {expected}",
        );
    }
}
