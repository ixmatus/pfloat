//! MPFR differential: `BigFloat::hypot` matches `rug::Float::hypot`
//! bit-for-bit across all five IEEE rounding modes.
//!
//! hypot evaluates `sqrt(x² + y²)` at an inflated Ziv working precision
//! (ADR-0056); the result is irrational except for Pythagorean pairs, so
//! `mpfr_oracle_for_mode` synthesizes the `NearestAway` oracle from a
//! high-precision MPFR result like the other Ziv-routed lanes. The
//! §9.2.1 special cases (infinity dominating NaN, the sNaN signal) are
//! covered by the in-crate unit tests; this lane confines itself to
//! finite signed integer pairs, including the near-equal operands that
//! exercise the sum-of-squares path.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_oracle_for_mode, next_i64_in, rug_from_i64,
    sweep_size, BIT_EXACT_ROUNDING_MODES, SWEEP_PRECISIONS,
};

#[test]
fn hypot_matches_mpfr_on_signed_i64_pairs() {
    let mut state: u64 = u64::from_le_bytes(*b"pflm1bhy");
    let cases = sweep_size();

    for &p in SWEEP_PRECISIONS {
        let cap = if p >= 64 {
            i64::MAX
        } else {
            (1_i64 << (p as i64 - 1)) - 1
        };
        for i in 0..cases {
            let x = next_i64_in(&mut state, -cap, cap);
            // Every fourth case forces near-equal operands (x ≈ y) to
            // stress the sum-of-squares path with no cancellation.
            let y = if i % 4 == 0 {
                x.saturating_add(1).min(cap)
            } else {
                next_i64_in(&mut state, -cap, cap)
            };
            for &mode in BIT_EXACT_ROUNDING_MODES {
                let bf_r = {
                    let x_bf = bigfloat_from_i64(x, p);
                    let y_bf = bigfloat_from_i64(y, p);
                    let (r, _status) = x_bf.hypot(&y_bf, mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = mpfr_oracle_for_mode(
                    |prec, round| {
                        let x_rg = rug_from_i64(x, prec);
                        let y_rg = rug_from_i64(y, prec);
                        rug::Float::with_val_round(prec, x_rg.hypot_ref(&y_rg), round).0
                    },
                    mode,
                    p,
                );
                assert_eq!(bf_r, rug_r, "hypot({x}, {y}) at p={p}, mode={mode:?}");
            }
        }
    }
}
