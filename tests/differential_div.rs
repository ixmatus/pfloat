//! MPFR differential: `BigFloat::div` matches `rug::Float`
//! division bit-for-bit at every tested precision and rounding mode.
//!
//! Divide-by-zero is exercised by the Kani harness in
//! `src/verify/div.rs`; this lane skips the zero-divisor case to
//! avoid the `Display(±inf) → rug parse` round-trip on infinity.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, rug_from_i64, sweep_size,
    ALL_ROUNDING_MODES, SWEEP_PRECISIONS,
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
            for &mode in ALL_ROUNDING_MODES {
                let bf_quot = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let b_bf = bigfloat_from_i64(b, p);
                    let (quot, _status) = a_bf.div(&b_bf, mode);
                    bigfloat_to_rug(&quot)
                };
                let rug_quot = {
                    let a_rg = rug_from_i64(a, p);
                    let b_rg = rug_from_i64(b, p);
                    let (quot, _ord) = rug::Float::with_val_round(
                        p,
                        &a_rg / &b_rg,
                        mpfr_round_of(mode)
                            .expect("NE-only lane: NearestEven has an MPFR equivalent (pf-suo)"),
                    );
                    quot
                };
                assert_eq!(bf_quot, rug_quot, "div({a}, {b}) at p={p}, mode={mode:?}");
            }
        }
    }
}
