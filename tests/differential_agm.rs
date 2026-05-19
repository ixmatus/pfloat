//! MPFR differential: `BigFloat::agm` matches `rug::Float::agm`.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, rug_from_i64, sweep_size,
    ALL_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
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
            for &mode in ALL_ROUNDING_MODES {
                let bf_r = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let b_bf = bigfloat_from_i64(b, p);
                    let (r, _status) = a_bf.agm(&b_bf, mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = {
                    let a_rg = rug_from_i64(a, p);
                    let b_rg = rug_from_i64(b, p);
                    let (r, _ord) = rug::Float::with_val_round(
                        p,
                        a_rg.agm_ref(&b_rg),
                        mpfr_round_of(mode)
                            .expect("NE-only lane: NearestEven has an MPFR equivalent (pf-suo)"),
                    );
                    r
                };
                assert_eq!(bf_r, rug_r, "agm({a}, {b}) at p={p}, mode={mode:?}");
            }
        }
    }
}
