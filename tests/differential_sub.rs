//! MPFR differential: `BigFloat::sub` matches `rug::Float` subtraction
//! bit-for-bit at every tested precision and rounding mode.
//!
//! Mirrors the integer-operand pattern from `differential_add.rs`.
//! ADR-0014 records the lane's gating.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, rug_from_i64, sweep_size,
    BIT_EXACT_ROUNDING_MODES, SWEEP_PRECISIONS,
};

#[test]
fn sub_matches_mpfr_on_i64_pairs() {
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
            let b = next_i64_in(&mut state, -cap, cap);
            for &mode in BIT_EXACT_ROUNDING_MODES {
                let bf_diff = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let b_bf = bigfloat_from_i64(b, p);
                    let (diff, _status) = a_bf.sub(&b_bf, mode);
                    bigfloat_to_rug(&diff)
                };
                let rug_diff = {
                    let a_rg = rug_from_i64(a, p);
                    let b_rg = rug_from_i64(b, p);
                    let (diff, _ord) =
                        rug::Float::with_val_round(p, &a_rg - &b_rg, mpfr_round_of(mode));
                    diff
                };
                assert_eq!(bf_diff, rug_diff, "sub({a}, {b}) at p={p}, mode={mode:?}");
            }
        }
    }
}
