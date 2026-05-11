//! MPFR differential: `BigFloat::pow` matches `rug::Float::pow`
//! bit-for-bit at every tested precision and rounding mode.
//!
//! Restricted to positive bases and small finite exponents. The
//! full §9.2.1 table (zero base, infinity base, negative base
//! with integer / non-integer exponent) is covered by the Kani
//! harnesses in `src/verify/pow.rs`.

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
fn pow_matches_mpfr_on_positive_base_small_exponent() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6c");
    let cases = sweep_size().min(500);

    for &p in SWEEP_PRECISIONS {
        for _ in 0..cases {
            // Base in [1, 100], exponent in [-10, 10] keeps the
            // result inside i64 range without testing the
            // overflow / underflow paths.
            let base = next_i64_in(&mut state, 1, 100);
            let exp = next_i64_in(&mut state, -10, 10);
            for &mode in ALL_ROUNDING_MODES {
                let bf_r = {
                    let b_bf = bigfloat_from_i64(base, p);
                    let e_bf = bigfloat_from_i64(exp, p);
                    let (r, _status) = b_bf.pow(&e_bf, mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = {
                    let b_rg = rug_from_i64(base, p);
                    let e_rg = rug_from_i64(exp, p);
                    let (r, _ord) =
                        rug::Float::with_val_round(p, b_rg.pow_ref(&e_rg), mpfr_round_of(mode));
                    r
                };
                assert_eq!(bf_r, rug_r, "pow({base}, {exp}) at p={p}, mode={mode:?}");
            }
        }
    }
}
