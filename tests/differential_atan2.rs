//! MPFR differential: `BigFloat::atan2` matches `rug::Float::atan2`.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, rug_from_i64, sweep_size,
    ALL_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};

#[test]
fn atan2_matches_mpfr_on_integer_pairs() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6d");
    let cases = sweep_size().min(1_000);

    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            let y = next_i64_in(&mut state, -(1_i64 << 30), 1_i64 << 30);
            let x = next_i64_in(&mut state, -(1_i64 << 30), 1_i64 << 30);
            // Skip (0, 0) — that is the dispatch boundary; Kani
            // harness covers it.
            if x == 0 && y == 0 {
                continue;
            }
            for &mode in ALL_ROUNDING_MODES {
                let bf_r = {
                    let y_bf = bigfloat_from_i64(y, p);
                    let x_bf = bigfloat_from_i64(x, p);
                    let (r, _status) = y_bf.atan2(&x_bf, mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = {
                    let y_rg = rug_from_i64(y, p);
                    let x_rg = rug_from_i64(x, p);
                    let (r, _ord) = rug::Float::with_val_round(
                        p,
                        y_rg.atan2_ref(&x_rg),
                        mpfr_round_of(mode)
                            .expect("NE-only lane: NearestEven has an MPFR equivalent (pf-suo)"),
                    );
                    r
                };
                assert_eq!(bf_r, rug_r, "atan2({y}, {x}) at p={p}, mode={mode:?}");
            }
        }
    }
}
