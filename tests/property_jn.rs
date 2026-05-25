//! Property tests for Bessel `J0`/`J1`/`Jn` (slice 6o). The
//! J/Y Wronskian is unavailable until `Y` ships (slice 6p), so the
//! load-bearing cross-tie here is the three-term recurrence
//! `J_{n−1}(x) + J_{n+1}(x) = (2n/x)·J_n(x)` (DLMF 10.6.1), which
//! binds three independently descended orders; plus argument
//! parity, the exact boundary values, and precision
//! self-consistency.
//!
//! Each `Jn` evaluation drives the Miller recurrence (whose seed
//! index grows with precision), so case counts are kept small
//! (these are correctness properties, not a fuzz sweep).

#![cfg(all(feature = "big", feature = "bessel"))]

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
    let d = BigFloat::try_from_i64_exact(den, p).unwrap();
    n.div(&d, RoundingMode::NearestEven).0
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(12))]

    /// Recurrence cross-tie `J_{n−1}(x)+J_{n+1}(x) = (2n/x)·J_n(x)`
    /// (DLMF 10.6.1), binding three independently descended orders.
    #[test]
    fn recurrence(n in 2i32..=6, num in 1i64..=30, den in 1i64..=4) {
        let p = 96u32;
        let x = rat(num, den, p);
        let (jm1, _) = x.jn(n - 1, RoundingMode::NearestEven);
        let (jp1, _) = x.jn(n + 1, RoundingMode::NearestEven);
        let (jn, _) = x.jn(n, RoundingMode::NearestEven);
        let (lhs, _) = jm1.add(&jp1, RoundingMode::NearestEven);
        let two_n = BigFloat::try_from_i64_exact(2 * i64::from(n), p).unwrap();
        let (r1, _) = two_n.mul(&jn, RoundingMode::NearestEven);
        let (rhs, _) = r1.div(&x, RoundingMode::NearestEven);
        prop_assert!(close_within(&lhs, &rhs, p - 8));
    }

    /// Argument parity `J_n(−x) = (−1)^n J_n(x)` (DLMF 10.11.1):
    /// the kernel reduces to `|x|` then flips the sign, so this is
    /// bit-exact.
    #[test]
    fn parity(n in 0i32..=6, num in 1i64..=30, den in 1i64..=4) {
        let p = 80u32;
        let xp = rat(num, den, p);
        let xn = rat(-num, den, p);
        let (fp, _) = xp.jn(n, RoundingMode::NearestEven);
        let (fn_, _) = xn.jn(n, RoundingMode::NearestEven);
        let expected = if n % 2 == 0 { fp.clone() } else { fp.negated() };
        prop_assert_eq!(fn_.partial_cmp(&expected).0, Some(Ordering::Equal));
    }

    /// Boundary: `J_0(0) = 1` exactly, `J_n(0) = +0` for `n ≥ 1`
    /// (DLMF 10.2.2).
    #[test]
    fn boundary(n in 1i32..=8, p in prop_oneof![Just(53u32), Just(113u32)]) {
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (j0, _) = z.j0(RoundingMode::NearestEven);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        prop_assert_eq!(j0.partial_cmp(&one).0, Some(Ordering::Equal));
        let (jn, _) = z.jn(n, RoundingMode::NearestEven);
        prop_assert!(jn.is_zero() && !jn.is_sign_negative());
    }

    /// Precision self-consistency: `Jn` at `p` agrees with `Jn` at
    /// `p + 96` rounded back to `p`. The argument is an exact dyadic
    /// rational (denominator a power of two) so both precisions
    /// evaluate the *same* real point — a non-dyadic argument would
    /// differ in its low bits between the two precisions, and near a
    /// zero of `Jn` that argument mismatch is amplified by
    /// `|Jn′/Jn|` into a spurious failure (this is a property of the
    /// test construction, not the kernel). Matches the pf-ok9 lesson
    /// already encoded in `property_yn::self_consistent`.
    #[test]
    fn self_consistent(
        n in 0i32..=5,
        num in 1i64..=20,
        den in prop_oneof![Just(1i64), Just(2), Just(4)],
    ) {
        let p = 96u32;
        let x_lo = rat(num, den, p);
        let x_hi = rat(num, den, p + 96);
        let (lo, _) = x_lo.jn(n, RoundingMode::NearestEven);
        let (hi_raw, _) = x_hi.jn(n, RoundingMode::NearestEven);
        let hi = hi_raw
            .round_to_precision(p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        prop_assert!(close_within(&lo, &hi, p - 2));
    }
}
