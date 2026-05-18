//! Property tests for Bessel `Y0`/`Y1`/`Yn` (slice 6p). The
//! load-bearing cross-tie is the **J/Y Wronskian** (DLMF 10.5.2)
//! `J_{n+1}(x)·Y_n(x) − J_n(x)·Y_{n+1}(x) = 2/(πx)`, the binding
//! identity slice 6o could not yet exercise (it shipped only the
//! `J` recurrence). It is checked in the π-free invariance form:
//! `(J_{n+1}Y_n − J_n Y_{n+1})·x` is independent of both the order
//! `n` and the argument `x` (it equals the constant `2/π`), so two
//! independent `(n, x)` evaluations must agree — no π constant is
//! needed, only the `bessel` feature. Plus order parity, the pole /
//! domain conventions, and precision self-consistency.
//!
//! Each `Yn` evaluation composes the 6o `J` Miller kernel plus the
//! log series or asymptotic (and, for `n ≥ 2`, the upward
//! recurrence), so case counts are kept small (these are
//! correctness properties, not a fuzz sweep).

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

/// `(J_{n+1}(x)·Y_n(x) − J_n(x)·Y_{n+1}(x)) · x` (DLMF 10.5.2 in
/// π-free invariance form: the value is `2/π` for every `n`, `x`).
fn wronskian_times_x(n: i32, x: &BigFloat) -> BigFloat {
    let (jn, _) = x.jn(n, RoundingMode::NearestEven);
    let (jn1, _) = x.jn(n + 1, RoundingMode::NearestEven);
    let (yn, _) = x.yn(n, RoundingMode::NearestEven);
    let (yn1, _) = x.yn(n + 1, RoundingMode::NearestEven);
    let (a, _) = jn1.mul(&yn, RoundingMode::NearestEven);
    let (b, _) = jn.mul(&yn1, RoundingMode::NearestEven);
    let (w, _) = a.sub(&b, RoundingMode::NearestEven);
    w.mul(x, RoundingMode::NearestEven).0
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(12))]

    /// J/Y Wronskian (DLMF 10.5.2): `(J_{n+1}Y_n − J_n Y_{n+1})·x`
    /// is independent of `n` and `x` (both evaluations equal `2/π`),
    /// binding the 6o `J` kernel to the new `Y` kernel.
    #[test]
    fn jy_wronskian_invariant(
        n1 in 0i32..=5,
        n2 in 0i32..=5,
        a1 in 1i64..=30, d1 in 1i64..=4,
        a2 in 1i64..=30, d2 in 1i64..=4,
    ) {
        let p = 128u32;
        let x1 = rat(a1, d1, p);
        let x2 = rat(a2, d2, p);
        let w1 = wronskian_times_x(n1, &x1);
        let w2 = wronskian_times_x(n2, &x2);
        prop_assert!(close_within(&w1, &w2, p - 12));
    }

    /// Order parity `Y_{−n}(x) = (−1)^n Y_n(x)` (DLMF 10.4.1): the
    /// kernel reduces on `m = |n|` then flips the sign, so this is
    /// bit-exact.
    #[test]
    fn order_parity(n in 0i32..=6, num in 1i64..=30, den in 1i64..=4) {
        let p = 96u32;
        let x = rat(num, den, p);
        let (pos, _) = x.yn(n, RoundingMode::NearestEven);
        let (neg, _) = x.yn(-n, RoundingMode::NearestEven);
        let expected = if n % 2 == 0 { pos.clone() } else { pos.negated() };
        prop_assert_eq!(neg.partial_cmp(&expected).0, Some(Ordering::Equal));
    }

    /// Domain: `Y_n(+0) = −∞` raising `DIV_BY_ZERO` (a pole), and a
    /// negative argument is `NaN` + `INVALID` (`Y` complex in ℝ).
    #[test]
    fn pole_and_domain(n in 0i32..=8, num in 1i64..=30, den in 1i64..=4) {
        let p = 80u32;
        let z = BigFloat::try_new_zero(Sign::Positive, p).unwrap();
        let (yp, sp) = z.yn(n, RoundingMode::NearestEven);
        prop_assert!(yp.is_infinite() && yp.is_sign_negative());
        prop_assert!(sp.div_by_zero());

        let xneg = rat(-num, den, p);
        let (yn, sn) = xneg.yn(n, RoundingMode::NearestEven);
        prop_assert!(yn.is_quiet_nan());
        prop_assert!(sn.invalid());
    }

    /// Precision self-consistency: `Yn` at `p` agrees with `Yn` at
    /// `p + 96` rounded back to `p`. The argument is an exact dyadic
    /// rational (denominator a power of two) so both precisions
    /// evaluate the *same* real point — a non-dyadic argument would
    /// differ in its low bits between the two precisions, and near a
    /// zero of `Yn` that argument mismatch is amplified by
    /// `|Yn′/Yn|` into a spurious failure (this is a property of the
    /// test construction, not the kernel).
    #[test]
    fn self_consistent(
        n in 0i32..=5,
        num in 1i64..=20,
        den in prop_oneof![Just(1i64), Just(2), Just(4)],
    ) {
        let p = 96u32;
        let x_lo = rat(num, den, p);
        let x_hi = rat(num, den, p + 96);
        let (lo, _) = x_lo.yn(n, RoundingMode::NearestEven);
        let (hi_raw, _) = x_hi.yn(n, RoundingMode::NearestEven);
        let hi = hi_raw
            .round_to_precision(p, RoundingMode::NearestEven)
            .unwrap()
            .0;
        prop_assert!(close_within(&lo, &hi, p - 2));
    }
}
