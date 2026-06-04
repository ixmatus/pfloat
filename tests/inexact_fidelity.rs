//! IEEE 754-2019 §7.6 `INEXACT`-flag fidelity for the transcendental
//! exp/log family and sin/cos (pf-njs5, ADR-0060).
//!
//! Two regimes the directed-pair libm shell surfaced, both with a
//! correctly rounded *value* but a wrong *flag* before the fix:
//!
//! - **Over-report**: a composed transcendental whose true result is
//!   exactly representable (`exp2(10) = 1024`, `exp10(2) = 100`,
//!   `log10(1000) = 3`, `ln(1) = 0`, `log2(2^k) = k`) set `INEXACT`
//!   because the kernel rounds internally. The fix dispatches the
//!   decidable exact-input set before the Ziv loop and returns
//!   `Status::OK`.
//!
//! - **Under-report**: a transcendental whose true result is
//!   irrational but rounds onto a grid value because the residual fell
//!   below the kernel's working precision (`exp(2^-1074) = 1.0` at low
//!   target precision) cleared `INEXACT`. The fix forces `INEXACT` on
//!   the transcendental fall-through: by Lindemann–Weierstrass /
//!   Gelfond–Schneider the result is irrational, so the flag is
//!   unconditionally correct there.
//!
//! Each assertion holds under every IEEE rounding mode.
//!
//! Run: `cargo test --test inexact_fidelity \
//!        --features std,big,exp-log,trig`

#![cfg(all(feature = "big", feature = "exp-log", feature = "trig"))]

use core::cmp::Ordering;
use pfloat::{BigFloat, RoundingMode};

const MODES: [RoundingMode; 5] = [
    RoundingMode::NearestEven,
    RoundingMode::NearestAway,
    RoundingMode::TowardZero,
    RoundingMode::TowardPositive,
    RoundingMode::TowardNegative,
];

fn from_i(n: i64, p: u32) -> BigFloat {
    BigFloat::try_from_i64_exact(n, p).expect("precision >= 1")
}

fn eq(a: &BigFloat, b: &BigFloat) -> bool {
    matches!(a.partial_cmp(b).0, Some(Ordering::Equal))
}

/// `2^k` (exact, k ≥ 0) by repeated doubling at precision `p`.
fn two_pow(k: u32, p: u32) -> BigFloat {
    let two = from_i(2, p);
    let mut x = from_i(1, p);
    for _ in 0..k {
        x = x.mul(&two, RoundingMode::NearestEven).0;
    }
    x
}

/// `2^-k` (exact, a power of two) by repeated halving at precision `p`.
fn two_pow_neg(k: u32, p: u32) -> BigFloat {
    let two = from_i(2, p);
    let mut x = from_i(1, p);
    for _ in 0..k {
        x = x.div(&two, RoundingMode::NearestEven).0;
    }
    x
}

/// `10^k` (exact for `p` ≥ sig-bits) by repeated multiply at `p`.
fn ten_pow(k: u32, p: u32) -> BigFloat {
    let ten = from_i(10, p);
    let mut x = from_i(1, p);
    for _ in 0..k {
        x = x.mul(&ten, RoundingMode::NearestEven).0;
    }
    x
}

// ---------------------------------------------------------------------
// Over-report: exact results must report INEXACT *clear*.
// ---------------------------------------------------------------------

#[test]
fn exp2_of_integer_is_exact_every_mode() {
    // 2^k is a power of two, exact at every precision and sign of k.
    for k in [-130_i64, -2, 0, 1, 3, 10, 53, 200] {
        let expected = if k >= 0 {
            two_pow(k as u32, 256)
        } else {
            two_pow_neg((-k) as u32, 256)
        };
        for m in MODES {
            let (v, s) = from_i(k, 64).exp2(m);
            assert!(
                !s.inexact(),
                "exp2({k}) is exactly 2^{k}; INEXACT must be clear (mode {m:?})"
            );
            assert!(eq(&v, &expected), "exp2({k}) value wrong (mode {m:?})");
        }
    }
}

#[test]
fn exp10_of_small_nonneg_integer_is_exact_every_mode() {
    // 10^k = 5^k·2^k is exact while 5^k fits the target precision.
    for k in [0_i64, 1, 2, 5, 10] {
        let expected = ten_pow(k as u32, 256);
        for m in MODES {
            let (v, s) = from_i(k, 64).exp10_round(256, m).expect("precision >= 1");
            assert!(
                !s.inexact(),
                "exp10({k}) = 10^{k} fits 256 bits; INEXACT must be clear (mode {m:?})"
            );
            assert!(eq(&v, &expected), "exp10({k}) value wrong (mode {m:?})");
        }
    }
}

#[test]
fn log10_of_power_of_ten_is_exact_every_mode() {
    // log10(10^k) = k, the only dyadic inputs with a rational log10
    // (k = 0 is x = 1 → 0).
    for k in [0_i64, 1, 2, 3, 6, 10] {
        let input = ten_pow(k as u32, 256);
        let expected = from_i(k, 53);
        for m in MODES {
            let (v, s) = input.log10_round(53, m).expect("precision >= 1");
            assert!(
                !s.inexact(),
                "log10(10^{k}) = {k} exactly; INEXACT must be clear (mode {m:?})"
            );
            assert!(eq(&v, &expected), "log10(10^{k}) value wrong (mode {m:?})");
        }
    }
}

#[test]
fn log2_of_power_of_two_is_exact_every_mode() {
    for k in [0_i64, 1, 3, 10, 64] {
        let input = if k >= 0 {
            two_pow(k as u32, 256)
        } else {
            two_pow_neg((-k) as u32, 256)
        };
        let expected = from_i(k, 53);
        for m in MODES {
            let (v, s) = input.log2_round(53, m).expect("precision >= 1");
            assert!(
                !s.inexact(),
                "log2(2^{k}) = {k} exactly; INEXACT must be clear (mode {m:?})"
            );
            assert!(eq(&v, &expected), "log2(2^{k}) value wrong (mode {m:?})");
        }
    }
}

#[test]
fn ln_of_one_is_exact_zero_every_mode() {
    let one = from_i(1, 256);
    let zero = from_i(0, 53);
    for m in MODES {
        let (v, s) = one.ln_round(53, m).expect("precision >= 1");
        assert!(
            !s.inexact(),
            "ln(1) = 0 exactly; INEXACT must be clear (mode {m:?})"
        );
        assert!(eq(&v, &zero), "ln(1) value wrong (mode {m:?})");
    }
}

// ---------------------------------------------------------------------
// Under-report: a sub-working-precision residual must still report
// INEXACT *set* (the result is transcendental, hence irrational).
// ---------------------------------------------------------------------

#[test]
fn exp_of_subnormal_tiny_is_inexact_every_mode() {
    // exp(2^-1074) = 1 + 2^-1074 + … rounds to 1.0 at target 53; the
    // residual 2^-1074 is far below the Ziv working precision, so the
    // kernel never observes a rounding — yet the true value ≠ 1.0.
    let tiny = two_pow_neg(1074, 64);
    for m in MODES {
        let (v, s) = tiny.exp_round(53, m);
        assert!(
            s.inexact(),
            "exp(2^-1074) ≠ 1 exactly; INEXACT must be set (mode {m:?})"
        );
        // Under directed-up modes the correctly rounded value is the
        // successor of 1.0; under the others it is 1.0. Either way it
        // must be finite and ≈ 1.
        assert!(
            v.is_normal(),
            "exp(2^-1074) must be finite normal (mode {m:?})"
        );
    }
}

#[test]
fn exp2_exp10_of_tiny_are_inexact_every_mode() {
    let tiny = two_pow_neg(1074, 64);
    for m in MODES {
        let (_, s2) = tiny.exp2_round(53, m).expect("precision >= 1");
        assert!(
            s2.inexact(),
            "exp2(2^-1074) is transcendental; INEXACT set (mode {m:?})"
        );
        let (_, s10) = tiny.exp10_round(53, m).expect("precision >= 1");
        assert!(
            s10.inexact(),
            "exp10(2^-1074) is transcendental; INEXACT set (mode {m:?})"
        );
    }
}

#[test]
fn sin_cos_of_huge_argument_are_inexact_every_mode() {
    // sin/cos of a large representable 2^k: the reduced residual can
    // collapse the result onto a grid value, but the true value is
    // irrational (Lindemann–Weierstrass), so INEXACT must be set.
    let big = two_pow(1000, 64);
    for m in MODES {
        let (vs, ss) = big.sin_round(53, m).expect("precision >= 1");
        let (vc, sc) = big.cos_round(53, m).expect("precision >= 1");
        assert!(
            ss.inexact(),
            "sin(2^1000) cannot be exact; INEXACT set (mode {m:?})"
        );
        assert!(
            sc.inexact(),
            "cos(2^1000) cannot be exact; INEXACT set (mode {m:?})"
        );
        assert!(
            !vs.is_nan() && !vc.is_nan(),
            "2^1000 is within the reduction table (mode {m:?})"
        );
    }
}

#[test]
fn ordinary_transcendentals_are_inexact_every_mode() {
    // A representable, non-exact input on each kernel: the result is
    // irrational and must report INEXACT under every mode.
    let two = from_i(2, 64);
    let three = from_i(3, 64);
    let half = two_pow_neg(1, 64); // 0.5, exact, but exp2(0.5)=√2 irrational
    let five_halves = {
        let (q, _) = from_i(5, 64).div(&two, RoundingMode::NearestEven); // 2.5
        q
    };
    for m in MODES {
        assert!(
            two.exp_round(53, m).1.inexact(),
            "exp(2) inexact (mode {m:?})"
        );
        assert!(
            two.ln_round(53, m).expect("p>=1").1.inexact(),
            "ln(2) inexact (mode {m:?})"
        );
        assert!(
            two.log10_round(53, m).expect("p>=1").1.inexact(),
            "log10(2) inexact (mode {m:?})"
        );
        // 3 is not a power of two, so log2(3) is irrational.
        assert!(
            three.log2_round(53, m).expect("p>=1").1.inexact(),
            "log2(3) inexact (mode {m:?})"
        );
        assert!(
            half.exp2_round(53, m).expect("p>=1").1.inexact(),
            "exp2(0.5)=√2 inexact (mode {m:?})"
        );
        assert!(
            five_halves.exp10_round(53, m).expect("p>=1").1.inexact(),
            "exp10(2.5) inexact (mode {m:?})"
        );
    }
}
