//! MPFR differential: `BigFloat::rootn` matches MPFR's `rootn`
//! bit-for-bit across all five IEEE rounding modes.
//!
//! rug 1.30.0 (MPFR 4.2.2) exposes the IEEE 754-2019 §9.2.1 n-th root
//! directly: `Float::root_ref(u32)` is `mpfr_rootn_ui` (n > 0) and
//! `Float::root_i_ref(i32)` is `mpfr_rootn_si` (signed n, including the
//! negative-order reciprocal). `mpfr_oracle_for_mode` synthesizes the
//! `NearestAway` oracle from a high-precision result, like the other
//! Ziv-routed lanes.
//!
//! A positive integer base's n-th root is an integer (perfect power) or
//! irrational, so it never lands on a tie. This lane sweeps positive
//! bases and small nonzero orders; the full §9.2.1 table (order zero,
//! signed zeros, negative-even domain errors, infinities) is covered by
//! the in-crate unit tests, since those produce NaN/∞ results that
//! `assert_eq!` cannot compare (NaN ≠ NaN).

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_oracle_for_mode, next_i64_in, rug_from_i64,
    sweep_size, BIT_EXACT_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};

#[test]
fn rootn_matches_mpfr_on_positive_base_small_order() {
    let mut state: u64 = u64::from_le_bytes(*b"pflm1brt");
    let cases = sweep_size().min(500);

    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            // Positive base, small nonzero order (n = 0 is the §9.2.1
            // domain error, unit-tested). Negative orders exercise the
            // reciprocal branch on both sides.
            let base = next_i64_in(&mut state, 1, 1000);
            let mut n = next_i64_in(&mut state, -8, 8) as i32;
            if n == 0 {
                n = 1;
            }
            for &mode in BIT_EXACT_ROUNDING_MODES {
                let bf_r = {
                    let b_bf = bigfloat_from_i64(base, p);
                    let (r, _status) = b_bf.rootn(n, mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = mpfr_oracle_for_mode(
                    |prec, round| {
                        let b_rg = rug_from_i64(base, prec);
                        if n > 0 {
                            rug::Float::with_val_round(prec, b_rg.root_ref(n as u32), round).0
                        } else {
                            rug::Float::with_val_round(prec, b_rg.root_i_ref(n), round).0
                        }
                    },
                    mode,
                    p,
                );
                assert_eq!(bf_r, rug_r, "rootn({base}, {n}) at p={p}, mode={mode:?}");
            }
        }
    }
}
