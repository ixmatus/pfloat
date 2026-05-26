//! MPFR differential: `BigFloat::cos` matches `rug::Float::cos`
//! bit-for-bit under every IEEE rounding mode (cos is correctly
//! rounded via Ziv per slice p1.26, ADR-0038).

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_oracle_for_mode, next_i64_in, rug_from_i64,
    sweep_size, BIT_EXACT_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};

#[test]
fn cos_matches_mpfr_on_small_integer_inputs() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6d");
    let cases = sweep_size().min(1_000);

    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            let a = next_i64_in(&mut state, -(1_i64 << 31), 1_i64 << 31);
            for &mode in BIT_EXACT_ROUNDING_MODES {
                let bf_r = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let (r, _status) = a_bf.cos(mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = mpfr_oracle_for_mode(
                    |prec, round| {
                        let a_rg = rug_from_i64(a, prec);
                        rug::Float::with_val_round(prec, a_rg.cos_ref(), round).0
                    },
                    mode,
                    p,
                );
                assert_eq!(bf_r, rug_r, "cos({a}) at p={p}, mode={mode:?}");
            }
        }
    }
}
