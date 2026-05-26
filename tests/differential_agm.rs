//! MPFR differential: `BigFloat::agm` matches `rug::Float::agm`
//! bit-for-bit under every IEEE rounding mode (agm is correctly
//! rounded via Ziv per slice p1.36, ADR-0038).

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_oracle_for_mode, next_i64_in, rug_from_i64,
    sweep_size, BIT_EXACT_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};

#[test]
fn agm_matches_mpfr_on_positive_integer_inputs() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6l");
    let cases = sweep_size().min(500);

    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            // Strictly positive inputs only; agm requires non-negative
            // operands and the zero case is exercised separately in
            // the property and unit suites.
            let a = next_i64_in(&mut state, 1, 1_000);
            let b = next_i64_in(&mut state, 1, 1_000);
            for &mode in BIT_EXACT_ROUNDING_MODES {
                let bf_r = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let b_bf = bigfloat_from_i64(b, p);
                    let (r, _status) = a_bf.agm(&b_bf, mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = mpfr_oracle_for_mode(
                    |prec, round| {
                        let a_rg = rug_from_i64(a, prec);
                        let b_rg = rug_from_i64(b, prec);
                        rug::Float::with_val_round(prec, a_rg.agm_ref(&b_rg), round).0
                    },
                    mode,
                    p,
                );
                assert_eq!(bf_r, rug_r, "agm({a}, {b}) at p={p}, mode={mode:?}");
            }
        }
    }
}
