//! MPFR differential for `BigFloat::y0` / `y1` / `yn` (Bessel
//! functions of the second kind, integer order, the
//! `differential_jn` "one file, several functions" precedent).
//!
//! rug 1.30 exposes MPFR `mpfr_y0` / `mpfr_y1` / `mpfr_yn`
//! (`y0_ref` / `y1_ref` / `yn_ref`), a genuine external oracle, so
//! this is a **bit-exact** lane (`assert_eq!`, the `differential_ei`
//! / `differential_jn` idiom) under `NearestEven`.
//!
//! Domain note: `Y` is real-valued only for `x > 0` (the Ci/li
//! convention; `Y` is complex off the positive axis, with a pole at
//! the origin), so every sweep argument is strictly positive — there
//! is no negative-argument arm (unlike `differential_jn`, where `J`
//! is entire).
//!
//! Cost note (`feedback_differential_lane_cost`): each `Y`
//! evaluation composes the 6o `J` kernel (Miller, whose seed index
//! grows with precision), the log series or the asymptotic, and the
//! `γ`/`ln` constants; the `n ≥ 2` orders add an upward recurrence.
//! `p = 1024` is therefore the cost driver and this lane is markedly
//! heavier than `differential_jn`. It runs the full
//! `TRANSCENDENTAL_PRECISIONS` (including `p = 1024`) as a
//! deliberately slow CI tier (the `differential_si` /
//! `differential_ci` posture, confirmed for 6p). The sweeps are a
//! small bounded representative set, not a wide random sweep.
//! Validate a fast subset locally (one function at `p = 53`); the
//! full lane is the CI slow tier.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, rug_from_i64, sweep_size,
    NEAREST_EVEN_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};
use pfloat::RoundingMode;

/// Exact positive dyadic-rational inputs (denominator a power of
/// two, so the constructed argument is bit-identical in pfloat and
/// rug at every precision). Spans the small-`x` log-series regime
/// near the pole and the moderate regime; all strictly positive (the
/// `Y` domain).
const DYADIC: &[(i64, i64)] = &[
    (1, 2),
    (1, 4),
    (3, 8),
    (1, 16),
    (5, 4),
    (9, 2),
    (15, 2),
    (21, 4),
];

#[test]
fn y0_matches_mpfr() {
    let mut state: u64 = u64::from_le_bytes(*b"pf6py0__");
    let cases = sweep_size().min(10);
    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            let a = next_i64_in(&mut state, 1, 40);
            for &mode in NEAREST_EVEN_ROUNDING_MODES {
                let bf_r = {
                    let (r, _s) = bigfloat_from_i64(a, p).y0(mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r =
                    rug::Float::with_val_round(p, rug_from_i64(a, p).y0_ref(), mpfr_round_of(mode))
                        .0;
                assert_eq!(bf_r, rug_r, "Y0({a}) at p={p}, mode={mode:?}");
            }
        }
    }
}

#[test]
fn y1_matches_mpfr() {
    let mut state: u64 = u64::from_le_bytes(*b"pf6py1__");
    let cases = sweep_size().min(10);
    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            let a = next_i64_in(&mut state, 1, 40);
            for &mode in NEAREST_EVEN_ROUNDING_MODES {
                let bf_r = {
                    let (r, _s) = bigfloat_from_i64(a, p).y1(mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r =
                    rug::Float::with_val_round(p, rug_from_i64(a, p).y1_ref(), mpfr_round_of(mode))
                        .0;
                assert_eq!(bf_r, rug_r, "Y1({a}) at p={p}, mode={mode:?}");
            }
        }
    }
}

#[test]
fn yn_matches_mpfr() {
    let mut state: u64 = u64::from_le_bytes(*b"pf6pyn__");
    let cases = sweep_size().min(8);
    for &p in TRANSCENDENTAL_PRECISIONS {
        for n in [2i32, 3, 5] {
            for _ in 0..cases {
                let a = next_i64_in(&mut state, 1, 40);
                for &mode in NEAREST_EVEN_ROUNDING_MODES {
                    let bf_r = {
                        let (r, _s) = bigfloat_from_i64(a, p).yn(n, mode);
                        bigfloat_to_rug(&r)
                    };
                    let rug_r = rug::Float::with_val_round(
                        p,
                        rug_from_i64(a, p).yn_ref(n),
                        mpfr_round_of(mode),
                    )
                    .0;
                    assert_eq!(bf_r, rug_r, "Y{n}({a}) at p={p}, mode={mode:?}");
                }
            }
        }
    }
}

/// Negative order (the `(−1)^n` order-parity path) and exact
/// positive dyadic rationals (the small-`x` log-series regime near
/// the pole), cross-checked bit-exact for `Y0`/`Y1`/`Y2`/`Y3`
/// against MPFR.
#[test]
fn yn_negative_order_and_dyadic_matches_mpfr() {
    for &p in TRANSCENDENTAL_PRECISIONS {
        let mode = RoundingMode::NearestEven;
        // Negative orders: parity reduction vs MPFR mpfr_yn directly
        // (MPFR's yn_ref accepts negative n and applies the same
        // (−1)^n parity).
        for &(num, den) in DYADIC {
            let x_bf = {
                let (q, _) = bigfloat_from_i64(num, p)
                    .div(&bigfloat_from_i64(den, p), RoundingMode::NearestEven);
                q
            };
            let x_rg = rug_from_i64(num, p) / rug_from_i64(den, p);
            for n in [0i32, 1, 2, 3] {
                let bf_r = {
                    let (r, _s) = x_bf.yn(n, mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = rug::Float::with_val_round(p, x_rg.yn_ref(n), mpfr_round_of(mode)).0;
                assert_eq!(bf_r, rug_r, "Y{n}({num}/{den}) at p={p}");

                let bf_neg = {
                    let (r, _s) = x_bf.yn(-n, mode);
                    bigfloat_to_rug(&r)
                };
                let rug_neg = rug::Float::with_val_round(p, x_rg.yn_ref(-n), mpfr_round_of(mode)).0;
                assert_eq!(bf_neg, rug_neg, "Y(-{n})({num}/{den}) at p={p}");
            }
        }
    }
}
