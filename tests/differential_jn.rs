//! MPFR differential for `BigFloat::j0` / `j1` / `jn` (Bessel
//! functions of the first kind, integer order, the
//! `differential_bi` "one file, several functions" precedent).
//!
//! rug 1.30 exposes MPFR `mpfr_j0` / `mpfr_j1` / `mpfr_jn`
//! (`j0_ref` / `j1_ref` / `jn_ref`), a genuine external oracle, so
//! this is a **bit-exact** lane (`assert_eq!`, the `differential_ei`
//! idiom) under `NearestEven`.
//!
//! Cost note (`feedback_differential_lane_cost`): the integer sweep
//! lands almost entirely in the Miller recurrence regime, whose seed
//! index `M` grows with target precision, so `p = 1024` is the cost
//! driver. The sweeps are therefore a small bounded representative
//! set, not a wide random sweep; this is a deliberately slow CI tier
//! (the `differential_si` / `differential_ci` posture). Validate a
//! fast subset locally (one function at `p = 53`); the full lane is
//! the CI slow tier.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, rug_from_i64, sweep_size,
    NEAREST_EVEN_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};
use pfloat::RoundingMode;

/// Exact dyadic-rational inputs (denominator a power of two, so the
/// constructed argument is bit-identical in pfloat and rug at every
/// precision). Spans the tiny `|x| < 1` Maclaurin regime and the
/// moderate Miller regime, both signs.
const DYADIC: &[(i64, i64)] = &[
    (1, 2),
    (1, 4),
    (3, 8),
    (1, 16),
    (5, 4),
    (9, 2),
    (-3, 4),
    (-7, 8),
    (-11, 2),
];

#[test]
fn j0_matches_mpfr() {
    let mut state: u64 = u64::from_le_bytes(*b"pf6oj0__");
    let cases = sweep_size().min(10);
    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            let a = next_i64_in(&mut state, 1, 40);
            for &mode in NEAREST_EVEN_ROUNDING_MODES {
                let bf_r = {
                    let (r, _s) = bigfloat_from_i64(a, p).j0(mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r =
                    rug::Float::with_val_round(p, rug_from_i64(a, p).j0_ref(), mpfr_round_of(mode))
                        .0;
                assert_eq!(bf_r, rug_r, "J0({a}) at p={p}, mode={mode:?}");
            }
        }
    }
}

#[test]
fn j1_matches_mpfr() {
    let mut state: u64 = u64::from_le_bytes(*b"pf6oj1__");
    let cases = sweep_size().min(10);
    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            let a = next_i64_in(&mut state, 1, 40);
            for &mode in NEAREST_EVEN_ROUNDING_MODES {
                let bf_r = {
                    let (r, _s) = bigfloat_from_i64(a, p).j1(mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r =
                    rug::Float::with_val_round(p, rug_from_i64(a, p).j1_ref(), mpfr_round_of(mode))
                        .0;
                assert_eq!(bf_r, rug_r, "J1({a}) at p={p}, mode={mode:?}");
            }
        }
    }
}

#[test]
fn jn_matches_mpfr() {
    let mut state: u64 = u64::from_le_bytes(*b"pf6ojn__");
    let cases = sweep_size().min(8);
    for &p in TRANSCENDENTAL_PRECISIONS {
        for n in [2i32, 3, 5] {
            for _ in 0..cases {
                let a = next_i64_in(&mut state, 1, 40);
                for &mode in NEAREST_EVEN_ROUNDING_MODES {
                    let bf_r = {
                        let (r, _s) = bigfloat_from_i64(a, p).jn(n, mode);
                        bigfloat_to_rug(&r)
                    };
                    let rug_r = rug::Float::with_val_round(
                        p,
                        rug_from_i64(a, p).jn_ref(n),
                        mpfr_round_of(mode),
                    )
                    .0;
                    assert_eq!(bf_r, rug_r, "J{n}({a}) at p={p}, mode={mode:?}");
                }
            }
        }
    }
}

/// Negative `x` (the `(−1)^n` argument-parity path) and exact dyadic
/// rationals (the tiny-`|x|` regime), cross-checked bit-exact for
/// `J0`/`J1`/`J2`/`J3` against MPFR.
#[test]
fn jn_negative_and_dyadic_matches_mpfr() {
    let mut state: u64 = u64::from_le_bytes(*b"pf6oneg_");
    let cases = sweep_size().min(8);
    for &p in TRANSCENDENTAL_PRECISIONS {
        // Negative integers: parity reduction vs MPFR directly.
        for _ in 0..cases {
            let a = next_i64_in(&mut state, -30, -1);
            let mode = RoundingMode::NearestEven;
            for n in [0i32, 1, 2, 3] {
                let bf_r = {
                    let (r, _s) = bigfloat_from_i64(a, p).jn(n, mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = rug::Float::with_val_round(
                    p,
                    rug_from_i64(a, p).jn_ref(n),
                    mpfr_round_of(mode),
                )
                .0;
                assert_eq!(bf_r, rug_r, "J{n}({a}) at p={p}");
            }
        }
        // Exact dyadic rationals: identical inputs both sides.
        for &(num, den) in DYADIC {
            let mode = RoundingMode::NearestEven;
            let x_bf = {
                let (q, _) = bigfloat_from_i64(num, p)
                    .div(&bigfloat_from_i64(den, p), RoundingMode::NearestEven);
                q
            };
            let x_rg = rug_from_i64(num, p) / rug_from_i64(den, p);
            for n in [0i32, 1, 2, 3] {
                let bf_r = {
                    let (r, _s) = x_bf.jn(n, mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = rug::Float::with_val_round(p, x_rg.jn_ref(n), mpfr_round_of(mode)).0;
                assert_eq!(bf_r, rug_r, "J{n}({num}/{den}) at p={p}");
            }
        }
    }
}
