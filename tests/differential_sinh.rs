//! MPFR differential: `BigFloat::sinh` matches `rug::Float::sinh`
//! bit-for-bit under every IEEE rounding mode (sinh is correctly
//! rounded via Ziv per slice p1.27, ADR-0038).

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_oracle_for_mode, next_i64_in, rug_from_i64,
    sweep_size, BIT_EXACT_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};

#[test]
fn sinh_matches_mpfr_on_small_integer_inputs() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6d");
    let cases = sweep_size().min(1_000);

    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            // sinh grows exponentially; cap at |x| ≤ 30 to keep the
            // result inside f64 range and the working precision boost
            // modest.
            let a = next_i64_in(&mut state, -30, 30);
            for &mode in BIT_EXACT_ROUNDING_MODES {
                let bf_r = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let (r, _status) = a_bf.sinh(mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = mpfr_oracle_for_mode(
                    |prec, round| {
                        let a_rg = rug_from_i64(a, prec);
                        rug::Float::with_val_round(prec, a_rg.sinh_ref(), round).0
                    },
                    mode,
                    p,
                );
                assert_eq!(bf_r, rug_r, "sinh({a}) at p={p}, mode={mode:?}");
            }
        }
    }
}
