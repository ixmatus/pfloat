//! Property tests for `BigFloat::zeta` (Riemann zeta, real
//! argument; slice 6r).
//!
//! The load-bearing cross-ties are the **π-free even-integer
//! identities** derived from `ζ(2n) = (2π)^{2n}|B_{2n}|/(2·(2n)!)`
//! (DLMF 25.6.2) by eliminating the common `π` power:
//!
//! - `ζ(4) = (2/5)·ζ(2)²`   (from `(π⁴/90) = (2/5)(π²/6)²`)
//! - `ζ(6) = (4/7)·ζ(2)·ζ(4)` (from `(π⁶/945) = (4/7)(π²/6)(π⁴/90)`)
//!
//! These bind several independent Borwein-core evaluations with no
//! transcendental constant on either side (the constant-free
//! Wronskian-analog form, the `property_ik` precedent; only the
//! `zeta` feature is needed). Plus exact-rational value pins that
//! bind the functional-equation lane and the special-value dispatch
//! to integer arithmetic (`ζ(0)=−1/2`, `ζ(−1)=−1/12`, `ζ(−3)=1/120`,
//! the trivial zeros `ζ(−2)=ζ(−4)=0`), the pole / NaN / ±∞ domain
//! conventions, and precision self-consistency on exact dyadic
//! arguments **from the start** (the pf-ok9 lesson) spanning all
//! three code paths.

#![cfg(all(feature = "big", feature = "zeta"))]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Sign};
use proptest::prelude::*;

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
    if den == 1 {
        return n;
    }
    let d = BigFloat::try_from_i64_exact(den, p).unwrap();
    n.div(&d, RoundingMode::NearestEven).0
}

fn ne(x: &BigFloat) -> BigFloat {
    let (r, _) = x.zeta(RoundingMode::NearestEven);
    r
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(10))]

    /// π-free even-integer cross-ties (DLMF 25.6.2 with the common
    /// `π` power eliminated), binding independent Borwein-core
    /// evaluations: `ζ(4) = (2/5)·ζ(2)²` and
    /// `ζ(6) = (4/7)·ζ(2)·ζ(4)`. No transcendental constant appears;
    /// a coefficient or sign error in the accelerator fails here.
    #[test]
    fn even_integer_pi_free_identities(p in 64u32..=200) {
        let z2 = ne(&rat(2, 1, p));
        let z4 = ne(&rat(4, 1, p));
        let z6 = ne(&rat(6, 1, p));

        // ζ(4) = (2/5)·ζ(2)².
        let (z2_sq, _) = z2.mul(&z2, RoundingMode::NearestEven);
        let (lhs4, _) = z2_sq.mul(&rat(2, 5, p), RoundingMode::NearestEven);
        prop_assert!(close_within(&z4, &lhs4, p - 8), "ζ(4)=(2/5)ζ(2)²");

        // ζ(6) = (4/7)·ζ(2)·ζ(4).
        let (z2z4, _) = z2.mul(&z4, RoundingMode::NearestEven);
        let (lhs6, _) = z2z4.mul(&rat(4, 7, p), RoundingMode::NearestEven);
        prop_assert!(close_within(&z6, &lhs6, p - 8), "ζ(6)=(4/7)ζ(2)ζ(4)");
    }

    /// Exact-rational value pins binding the functional-equation
    /// lane and the special-value dispatch to integer arithmetic:
    /// ζ(0)=−1/2 (DLMF 25.6.1), ζ(−1)=−1/12, ζ(−3)=1/120
    /// (DLMF 25.6.3 via the FE), and the trivial zeros ζ(−2)=ζ(−4)=0
    /// (DLMF 25.6.4, special-cased exactly).
    #[test]
    fn exact_rational_pins(p in 64u32..=200) {
        let zero = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (z0, _) = BigFloat::try_new_zero(Sign::Positive, p)
            .unwrap()
            .zeta(RoundingMode::NearestEven);
        prop_assert!(close_within(&z0, &rat(-1, 2, p), p - 4), "ζ(0)=−1/2");

        let zm1 = ne(&rat(-1, 1, p));
        prop_assert!(close_within(&zm1, &rat(-1, 12, p), p - 8), "ζ(−1)=−1/12");

        let zm3 = ne(&rat(-3, 1, p));
        prop_assert!(close_within(&zm3, &rat(1, 120, p), p - 8), "ζ(−3)=1/120");

        for k in [-2i64, -4] {
            let zt = ne(&rat(k, 1, p));
            prop_assert_eq!(
                zt.partial_cmp(&zero).0,
                Some(Ordering::Equal),
                "ζ({}) trivial zero",
                k
            );
        }
    }

    /// Precision self-consistency: value at `p` agrees with value
    /// at `p + 96` rounded back. The argument is an exact dyadic
    /// rational (power-of-two denominator) **from the start** (the
    /// pf-ok9 lesson) so both precisions evaluate the same real
    /// point; spans s > 1 (Borwein), 0 < s < 1 (Borwein, critical
    /// strip), and s < 0 (functional equation).
    #[test]
    fn self_consistent(
        num in -20i64..=20,
        den in prop_oneof![Just(1i64), Just(2), Just(4), Just(8)],
    ) {
        // Skip the pole s = 1 and the s = 0 special-case (covered
        // elsewhere); keep genuine finite-regime points.
        prop_assume!(num != den && num != 0);
        let p = 96u32;
        let x_lo = rat(num, den, p);
        let x_hi = rat(num, den, p + 96);
        let lo = ne(&x_lo);
        let (hi_raw, _) = x_hi.zeta(RoundingMode::NearestEven);
        let hi = hi_raw
            .round_to_precision(p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        prop_assert!(close_within(&lo, &hi, p - 2), "ζ({num}/{den}) p vs p+96");
    }

    /// Domain conventions: the pole ζ(1)=+∞ + DIV_BY_ZERO; quiet
    /// NaN propagates; signaling NaN raises INVALID; ζ(+∞)=1;
    /// ζ(−∞)=NaN + INVALID (the unbounded-non-converging convention).
    #[test]
    fn domain_conventions(p in 24u32..=160) {
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let (zp, sp) = one.zeta(RoundingMode::NearestEven);
        prop_assert!(zp.is_infinite() && zp.is_sign_positive());
        prop_assert!(sp.div_by_zero());

        let q = BigFloat::try_new_quiet_nan(Sign::Positive, p, &[]).unwrap();
        let (zq, sq) = q.zeta(RoundingMode::NearestEven);
        prop_assert!(zq.is_quiet_nan());
        prop_assert!(!sq.invalid());

        let sn = BigFloat::try_new_signaling_nan(Sign::Positive, p, &[]).unwrap();
        let (zsn, ssn) = sn.zeta(RoundingMode::NearestEven);
        prop_assert!(zsn.is_quiet_nan());
        prop_assert!(ssn.invalid());

        let pinf = BigFloat::try_new_infinity(Sign::Positive, p).unwrap();
        let (zi, _) = pinf.zeta(RoundingMode::NearestEven);
        prop_assert_eq!(zi.partial_cmp(&one).0, Some(Ordering::Equal), "ζ(+∞)=1");

        let ninf = BigFloat::try_new_infinity(Sign::Negative, p).unwrap();
        let (zni, sni) = ninf.zeta(RoundingMode::NearestEven);
        prop_assert!(zni.is_quiet_nan());
        prop_assert!(sni.invalid());
    }
}
