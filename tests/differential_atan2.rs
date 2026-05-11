//! MPFR differential: `BigFloat::atan2` matches `rug::Float::atan2`.

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
fn atan2_matches_mpfr_on_integer_pairs() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6d");
    let cases = sweep_size().min(1_000);

    for &p in SWEEP_PRECISIONS {
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
                    let (r, _ord) =
                        rug::Float::with_val_round(p, y_rg.atan2_ref(&x_rg), mpfr_round_of(mode));
                    r
                };
                assert_eq!(bf_r, rug_r, "atan2({y}, {x}) at p={p}, mode={mode:?}");
            }
        }
    }
}
