//! Correctness regression for the beta exact-value-defeats-Ziv defect
//! (pf-umlm). For positive-integer and case-4 (pole-cancellation) inputs,
//! `β` is rational; the old exp(lgamma) Ziv path returned the
//! exactly-dyadic ones (β(1,2ᵏ)=2⁻ᵏ, β(1,1)=1, B(−1,1)=−1) off by a ULP
//! under directed modes and over-reported INEXACT. The construct-and-check
//! dispatch builds the exact rational and divides once, so dyadic outputs
//! are exact with INEXACT clear and non-dyadic ones are correctly rounded
//! with INEXACT set. Every assertion holds under all five rounding modes.
//!
//! The breadth check is an oracle sweep: the p=53 (value, INEXACT) must
//! match the kernel's own p=200 result rounded down to 53.
//!
//! Run: `cargo test --test beta_exact_fidelity --features std,big,specials`
#![cfg(all(feature = "big", feature = "specials"))]

use core::cmp::Ordering;
use pfloat::{BigFloat, RoundingMode};

const MODES: [RoundingMode; 5] = [
    RoundingMode::NearestEven,
    RoundingMode::NearestAway,
    RoundingMode::TowardZero,
    RoundingMode::TowardPositive,
    RoundingMode::TowardNegative,
];

fn fi(n: i64) -> BigFloat {
    BigFloat::try_from_i64_exact(n, 64).expect("p>=1")
}
fn two_pow_neg(k: u32) -> BigFloat {
    let two = fi(2);
    let mut x = fi(1);
    for _ in 0..k {
        x = x.div(&two, RoundingMode::NearestEven).0;
    }
    x
}
fn eqv(a: &BigFloat, b: &BigFloat) -> bool {
    matches!(a.partial_cmp(b).0, Some(Ordering::Equal))
}

#[test]
fn beta_dyadic_integer_outputs_are_exact_every_mode() {
    // β(1,1) = 1 and β(1,2ᵏ) = β(2ᵏ,1) = 2⁻ᵏ are exactly representable:
    // correct value, INEXACT clear, under every mode.
    for m in MODES {
        let (v, s) = fi(1).beta_round(&fi(1), 53, m).expect("p>=1");
        assert!(!s.inexact(), "β(1,1)=1 INEXACT must be clear (mode {m:?})");
        assert!(eqv(&v, &fi(1)), "β(1,1) value (mode {m:?})");
        for k in 1..=12u32 {
            let n = 1i64 << k;
            let expected = two_pow_neg(k);
            let (v1, s1) = fi(1).beta_round(&fi(n), 53, m).expect("p>=1");
            assert!(
                !s1.inexact(),
                "β(1,{n})=2^-{k} INEXACT must be clear (mode {m:?})"
            );
            assert!(eqv(&v1, &expected), "β(1,{n}) value (mode {m:?})");
            // symmetric
            let (v2, s2) = fi(n).beta_round(&fi(1), 53, m).expect("p>=1");
            assert!(
                !s2.inexact(),
                "β({n},1)=2^-{k} INEXACT must be clear (mode {m:?})"
            );
            assert!(eqv(&v2, &expected), "β({n},1) value (mode {m:?})");
        }
    }
}

#[test]
fn beta_nondyadic_integer_outputs_set_inexact_every_mode() {
    // Rational but non-dyadic: INEXACT must be set (and the value
    // correctly rounded, checked by the oracle sweep below).
    let cases = [(2i64, 2i64), (2, 3), (3, 5), (1, 3), (1, 6), (4, 6)];
    for m in MODES {
        for (a, b) in cases {
            let (_, s) = fi(a).beta_round(&fi(b), 53, m).expect("p>=1");
            assert!(
                s.inexact(),
                "β({a},{b}) non-dyadic INEXACT set (mode {m:?})"
            );
        }
    }
}

#[test]
fn beta_case4_dyadic_is_exact_every_mode() {
    // Case 4 (pole cancellation): B(−1,1) = −1 is dyadic → exact, INEXACT
    // clear; B(−5,5) = −1/5 is non-dyadic → INEXACT set.
    let neg_one = fi(-1);
    for m in MODES {
        let (v, s) = fi(-1).beta_round(&fi(1), 53, m).expect("p>=1");
        assert!(!s.inexact(), "B(-1,1)=-1 INEXACT clear (mode {m:?})");
        assert!(eqv(&v, &neg_one), "B(-1,1) value (mode {m:?})");
        let (_, s5) = fi(-5).beta_round(&fi(5), 53, m).expect("p>=1");
        assert!(
            s5.inexact(),
            "B(-5,5)=-1/5 non-dyadic INEXACT set (mode {m:?})"
        );
    }
}

/// Oracle sweep: p=53 (value, INEXACT) must match the kernel's own p=200
/// result rounded down to 53, across a grid of integer / case-4 inputs.
#[test]
fn beta_oracle_sweep_no_under_over_or_value_bugs() {
    let mut pairs: Vec<(BigFloat, BigFloat)> = vec![];
    for a in 1..=16i64 {
        for b in 1..=16i64 {
            pairs.push((fi(a), fi(b)));
        }
    }
    for k in 1..=12u32 {
        let n = 1i64 << k;
        pairs.push((fi(1), fi(n)));
        pairs.push((fi(n), fi(1)));
    }
    // case 4: B(−n, m), 1 ≤ m ≤ n
    for n in 1..=14i64 {
        for m in 1..=n {
            pairs.push((fi(-n), fi(m)));
        }
    }

    let mut under = 0usize;
    let mut over = 0usize;
    let mut valbug = 0usize;
    let mut first: Option<String> = None;
    for (a, b) in &pairs {
        let vo = a
            .beta_round(b, 200, RoundingMode::NearestEven)
            .expect("p>=1")
            .0;
        for m in MODES {
            let (vt, st) = a.beta_round(b, 53, m).expect("p>=1");
            let (vor, sor) = vo.round_to_precision(53, m).expect("p>=1");
            let flag_bad = st.inexact() != sor.inexact();
            let val_bad = !eqv(&vt, &vor);
            if flag_bad {
                if st.inexact() {
                    over += 1;
                } else {
                    under += 1;
                }
            }
            if val_bad {
                valbug += 1;
            }
            if (flag_bad || val_bad) && first.is_none() {
                first = Some(format!(
                    "β({},{}) [{m:?}] kernel=({}, INEXACT={}) oracle=({}, INEXACT={})",
                    a.to_f64(),
                    b.to_f64(),
                    vt.to_f64(),
                    st.inexact(),
                    vor.to_f64(),
                    sor.inexact()
                ));
            }
        }
    }
    assert!(
        under == 0 && over == 0 && valbug == 0,
        "beta oracle sweep: under={under} over={over} valbug={valbug}; first={first:?}"
    );
}
