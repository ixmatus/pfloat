//! MPFR differential: `BigFloat::div` matches `rug::Float`
//! division bit-for-bit at every tested precision and IEEE
//! rounding mode (the arithmetic core is correctly rounded under
//! every mode by construction; `differential_pow.rs` is the
//! template for the `NearestAway` synthesis pattern this lane
//! now consumes through `mpfr_oracle_for_mode`).
//!
//! Phase 1f slice p1.23 widened this lane from NE-only
//! (`ALL_ROUNDING_MODES`) to all five IEEE modes
//! (`BIT_EXACT_ROUNDING_MODES`) via the shared
//! `mpfr_oracle_for_mode` helper (ADR-0038). The `NearestAway` arm
//! routes through the p+128 synthesis pattern that
//! `differential_pow` already used.
//!
//! Divide-by-zero is exercised by the Kani harness in
//! `src/verify/div.rs`; this lane skips the zero-divisor case to
//! avoid the `Display(±inf) → rug parse` round-trip on infinity.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_oracle_for_mode, next_i64_in, rug_from_i64,
    sweep_size, BIT_EXACT_ROUNDING_MODES, SWEEP_PRECISIONS,
};

#[test]
fn div_matches_mpfr_on_i64_pairs() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6b");
    let cases = sweep_size();

    for &p in SWEEP_PRECISIONS {
        let cap = if p >= 64 {
            i64::MAX
        } else {
            (1_i64 << (p as i64 - 2)) - 1
        };
        for _ in 0..cases {
            let a = next_i64_in(&mut state, -cap, cap);
            let mut b = next_i64_in(&mut state, -cap, cap);
            while b == 0 {
                b = next_i64_in(&mut state, -cap, cap);
            }
            for &mode in BIT_EXACT_ROUNDING_MODES {
                let bf_quot = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let b_bf = bigfloat_from_i64(b, p);
                    let (quot, _status) = a_bf.div(&b_bf, mode);
                    bigfloat_to_rug(&quot)
                };
                let rug_quot = mpfr_oracle_for_mode(
                    |prec, round| {
                        let a_rg = rug_from_i64(a, prec);
                        let b_rg = rug_from_i64(b, prec);
                        rug::Float::with_val_round(prec, &a_rg / &b_rg, round).0
                    },
                    mode,
                    p,
                );
                assert_eq!(bf_quot, rug_quot, "div({a}, {b}) at p={p}, mode={mode:?}");
            }
        }
    }
}
