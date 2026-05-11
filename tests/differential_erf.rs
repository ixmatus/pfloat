//! MPFR differential: `BigFloat::erf` matches `rug::Float::erf`.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, rug_from_i64, sweep_size,
    ALL_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};

#[test]
fn erf_matches_mpfr_on_small_integer_inputs() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6e");
    let cases = sweep_size().min(500);

    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            // erf saturates quickly; |x| ≤ 10 covers the interesting
            // range without hitting the asymptotic 1 plateau.
            let a = next_i64_in(&mut state, -10, 10);
            for &mode in ALL_ROUNDING_MODES {
                let bf_r = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let (r, _status) = a_bf.erf(mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = {
                    let a_rg = rug_from_i64(a, p);
                    let (r, _ord) =
                        rug::Float::with_val_round(p, a_rg.erf_ref(), mpfr_round_of(mode));
                    r
                };
                assert_eq!(bf_r, rug_r, "erf({a}) at p={p}, mode={mode:?}");
            }
        }
    }
}
