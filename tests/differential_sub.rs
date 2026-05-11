//! MPFR differential: `BigFloat::sub` matches `rug::Float` subtraction
//! bit-for-bit at every tested precision and rounding mode.
//!
//! Mirrors the integer-operand pattern from `differential_add.rs`.
//! ADR-0014 records the lane's gating.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, rug_from_i64, sweep_size,
    ALL_ROUNDING_MODES, SWEEP_PRECISIONS,
};

fn next_u64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn next_i64_in(state: &mut u64, lo: i64, hi: i64) -> i64 {
    debug_assert!(lo <= hi);
    let span = (hi - lo) as u64 + 1;
    lo + (next_u64(state) % span) as i64
}

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
            for &mode in ALL_ROUNDING_MODES {
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
