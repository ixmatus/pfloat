//! Property tests for modified Bessel `I0`/`I1`/`In` and
//! `K0`/`K1`/`Kn` (slice 6q). The load-bearing cross-tie is the
//! **DLMF 10.28.2 I/K identity**
//! `I_ν(x)·K_{ν+1}(x) + I_{ν+1}(x)·K_ν(x) = 1/x`, the modified-Bessel
//! analog of 6p's J/Y Wronskian, binding the `bessel_i` and
//! `bessel_k` kernels. It is checked in the constant-free invariance
//! form: `(I_ν K_{ν+1} + I_{ν+1} K_ν)·x` is independent of both the
//! order `ν` and the argument `x` — it equals the exact constant
//! `1`, so it is π-free and two independent `(ν, x)` evaluations must
//! agree (and equal `1`); only the `bessel` feature is needed. Plus
//! the order / argument parities (which differ from `J`/`Y`: `I`/`K`
//! are **even in order with no sign**), the pole / domain
//! conventions, and precision self-consistency.
//!
//! Each `I`/`K` evaluation composes a recurrence or log series (and
//! `K` composes `bessel_i` plus `γ`/`ln`), so case counts are kept
//! small (these are correctness properties, not a fuzz sweep). All
//! self-consistency arguments are exact dyadic rationals from the
//! start (the pf-ok9 lesson; a non-dyadic argument differs in its
//! low bits between `p` and `p+96` and spuriously fails near a
//! steep region).

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

/// `(I_ν(x)·K_{ν+1}(x) + I_{ν+1}(x)·K_ν(x)) · x` (DLMF 10.28.2 in
/// constant-free invariance form: the value is exactly `1` for every
/// `ν`, `x > 0`).
fn cross_tie_times_x(nu: i32, x: &BigFloat) -> BigFloat {
    let (i_nu, _) = x.in_(nu, RoundingMode::NearestEven);
    let (i_nu1, _) = x.in_(nu + 1, RoundingMode::NearestEven);
    let (k_nu, _) = x.kn(nu, RoundingMode::NearestEven);
    let (k_nu1, _) = x.kn(nu + 1, RoundingMode::NearestEven);
    let (a, _) = i_nu.mul(&k_nu1, RoundingMode::NearestEven);
    let (b, _) = i_nu1.mul(&k_nu, RoundingMode::NearestEven);
    let (w, _) = a.add(&b, RoundingMode::NearestEven);
    w.mul(x, RoundingMode::NearestEven).0
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(12))]

    /// DLMF 10.28.2: `(I_ν K_{ν+1} + I_{ν+1} K_ν)·x` is independent
    /// of `ν` and `x` and equals the exact constant `1`, binding the
    /// `bessel_i` and `bessel_k` kernels (π-free, so no constant
    /// import is needed; the `1` is checked directly too).
    #[test]
    fn ik_cross_tie_invariant(
        n1 in 0i32..=5,
        n2 in 0i32..=5,
        a1 in 1i64..=30, d1 in 1i64..=4,
        a2 in 1i64..=30, d2 in 1i64..=4,
    ) {
        let p = 128u32;
        let x1 = rat(a1, d1, p);
        let x2 = rat(a2, d2, p);
        let w1 = cross_tie_times_x(n1, &x1);
        let w2 = cross_tie_times_x(n2, &x2);
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        prop_assert!(close_within(&w1, &w2, p - 12));
        prop_assert!(close_within(&w1, &one, p - 12));
    }

    /// Order parity `I_{−n}(x) = I_n(x)` (DLMF 10.27.1): **even in
    /// order, no sign** (unlike `J`/`Y`'s `(−1)^n`). The kernel
    /// reduces on `m = |n|` and applies no order sign, so bit-exact.
    #[test]
    fn i_order_parity(n in 0i32..=6, num in 1i64..=30, den in 1i64..=4) {
        let p = 96u32;
        let x = rat(num, den, p);
        let (pos, _) = x.in_(n, RoundingMode::NearestEven);
        let (neg, _) = x.in_(-n, RoundingMode::NearestEven);
        prop_assert_eq!(neg.partial_cmp(&pos).0, Some(Ordering::Equal));
    }

    /// Argument parity `I_n(−x) = (−1)^n I_n(x)` (the `(x/2)^n`
    /// prefactor of DLMF 10.25.2). Bit-exact (the kernel negates
    /// exactly when `m` is odd and `x < 0`).
    #[test]
    fn i_argument_parity(n in 0i32..=6, num in 1i64..=30, den in 1i64..=4) {
        let p = 96u32;
        let xp = rat(num, den, p);
        let xn = rat(-num, den, p);
        let (pos, _) = xp.in_(n, RoundingMode::NearestEven);
        let (neg, _) = xn.in_(n, RoundingMode::NearestEven);
        let expected = if n % 2 == 0 { pos.clone() } else { pos.negated() };
        prop_assert_eq!(neg.partial_cmp(&expected).0, Some(Ordering::Equal));
    }

    /// Order parity `K_{−n}(x) = K_n(x)` (DLMF 10.27.3): **even in
    /// order, no sign** (unlike `Y`'s `(−1)^n`). Bit-exact.
    #[test]
    fn k_order_parity(n in 0i32..=6, num in 1i64..=30, den in 1i64..=4) {
        let p = 96u32;
        let x = rat(num, den, p);
        let (pos, _) = x.kn(n, RoundingMode::NearestEven);
        let (neg, _) = x.kn(-n, RoundingMode::NearestEven);
        prop_assert_eq!(neg.partial_cmp(&pos).0, Some(Ordering::Equal));
    }

    /// Domain conventions. `I` is entire: `I_n(−x)` is finite (the
    /// argument parity) and never `INVALID`. `K` is `x > 0` only:
    /// `K_n(+0) = +∞` raising `DIV_BY_ZERO` (a pole, positive — the
    /// opposite of `Y_n(+0) = −∞`), and a negative argument is
    /// `NaN` + `INVALID` (`K` complex in ℝ).
    #[test]
    fn domain_conventions(n in 0i32..=8, num in 1i64..=30, den in 1i64..=4) {
        let p = 80u32;

        // I: negative argument is finite, OK (entire).
        let xneg = rat(-num, den, p);
        let (iv, is) = xneg.in_(n, RoundingMode::NearestEven);
        prop_assert!(!iv.is_nan());
        prop_assert!(!is.invalid());

        // K: positive-zero pole, +∞ + DIV_BY_ZERO.
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (kp, ks) = z.kn(n, RoundingMode::NearestEven);
        prop_assert!(kp.is_infinite() && kp.is_sign_positive());
        prop_assert!(ks.div_by_zero());

        // K: negative argument is NaN + INVALID (complex in ℝ).
        let (kn, kns) = xneg.kn(n, RoundingMode::NearestEven);
        prop_assert!(kn.is_quiet_nan());
        prop_assert!(kns.invalid());
    }

    /// Precision self-consistency for both kernels: value at `p`
    /// agrees with value at `p + 96` rounded back to `p`. The
    /// argument is an exact dyadic rational (denominator a power of
    /// two) **from the start** (the pf-ok9 lesson) so both
    /// precisions evaluate the same real point.
    #[test]
    fn self_consistent(
        n in 0i32..=5,
        num in 1i64..=20,
        den in prop_oneof![Just(1i64), Just(2), Just(4), Just(8)],
    ) {
        let p = 96u32;
        let x_lo = rat(num, den, p);
        let x_hi = rat(num, den, p + 96);

        let (i_lo, _) = x_lo.in_(n, RoundingMode::NearestEven);
        let (i_hi_raw, _) = x_hi.in_(n, RoundingMode::NearestEven);
        let i_hi = i_hi_raw
            .round_to_precision(p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        prop_assert!(close_within(&i_lo, &i_hi, p - 2));

        let (k_lo, _) = x_lo.kn(n, RoundingMode::NearestEven);
        let (k_hi_raw, _) = x_hi.kn(n, RoundingMode::NearestEven);
        let k_hi = k_hi_raw
            .round_to_precision(p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        prop_assert!(close_within(&k_lo, &k_hi, p - 2));
    }
}
