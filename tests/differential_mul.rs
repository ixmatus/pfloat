//! MPFR differential: `BigFloat::mul` matches `rug::Float`
//! multiplication bit-for-bit at every tested precision and rounding
//! mode.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, round_ties_to_away,
    rug_from_i64, sweep_size, BIT_EXACT_ROUNDING_MODES, SWEEP_PRECISIONS,
};
use pfloat::RoundingMode;

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
            for &mode in BIT_EXACT_ROUNDING_MODES {
                let bf_prod = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let b_bf = bigfloat_from_i64(b, p);
                    let (prod, _status) = a_bf.mul(&b_bf, mode);
                    bigfloat_to_rug(&prod)
                };
                let rug_prod = {
                    let a_rg = rug_from_i64(a, p);
                    let b_rg = rug_from_i64(b, p);
                    if matches!(mode, RoundingMode::NearestAway) {
                        // MPFR has no roundTiesToAway; synthesize it
                        // from an exact high-precision product (pf-suo).
                        let hp = rug::Float::with_val(p + 128, &a_rg * &b_rg);
                        round_ties_to_away(&hp, p)
                    } else {
                        rug::Float::with_val_round(
                            p,
                            &a_rg * &b_rg,
                            mpfr_round_of(mode).expect("non-NearestAway has an MPFR equivalent"),
                        )
                        .0
                    }
                };
                assert_eq!(bf_prod, rug_prod, "mul({a}, {b}) at p={p}, mode={mode:?}");
            }
        }
    }
}
