//! MPFR differential: `BigFloat::mul` matches `rug::Float`
//! multiplication bit-for-bit at every tested precision and rounding
//! mode.

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
fn mul_matches_mpfr_on_i64_pairs() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6b");
    let cases = sweep_size();

    for &p in SWEEP_PRECISIONS {
        // Cap operands so a*b fits in i64 (the conversion source).
        // The actual rounding to p bits happens inside both pfloat
        // and rug; we only need the inputs to be representable as
        // i64 and the product to fit in p bits for the exact
        // pre-rounding form.
        let per_operand_bits = (p / 2).min(31);
        let cap = (1_i64 << per_operand_bits) - 1;
        for _ in 0..cases {
            let a = next_i64_in(&mut state, -cap, cap);
            let b = next_i64_in(&mut state, -cap, cap);
            for &mode in ALL_ROUNDING_MODES {
                let bf_prod = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let b_bf = bigfloat_from_i64(b, p);
                    let (prod, _status) = a_bf.mul(&b_bf, mode);
                    bigfloat_to_rug(&prod)
                };
                let rug_prod = {
                    let a_rg = rug_from_i64(a, p);
                    let b_rg = rug_from_i64(b, p);
                    let (prod, _ord) =
                        rug::Float::with_val_round(p, &a_rg * &b_rg, mpfr_round_of(mode));
                    prod
                };
                assert_eq!(bf_prod, rug_prod, "mul({a}, {b}) at p={p}, mode={mode:?}");
            }
        }
    }
}
