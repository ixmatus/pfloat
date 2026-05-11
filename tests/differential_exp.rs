//! MPFR differential: `BigFloat::exp` matches `rug::Float::exp`
//! bit-for-bit at every tested precision and rounding mode.
//!
//! Argument magnitude is capped at a small absolute range so the
//! result stays representable in i64 dynamic range and does not
//! overflow to ±∞ (which would force a Display round-trip through
//! "inf"; covered by the Kani harness).

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
fn exp_matches_mpfr_on_small_integer_inputs() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6c");
    let cases = sweep_size().min(1_000); // exp is expensive; cap CI sweep

    for &p in SWEEP_PRECISIONS {
        for _ in 0..cases {
            // |x| <= 30 keeps exp(x) inside f64 range and small
            // enough that the working-precision boost stays modest.
            let a = next_i64_in(&mut state, -30, 30);
            for &mode in ALL_ROUNDING_MODES {
                let bf_r = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let (r, _status) = a_bf.exp(mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = {
                    let a_rg = rug_from_i64(a, p);
                    let (r, _ord) =
                        rug::Float::with_val_round(p, a_rg.exp_ref(), mpfr_round_of(mode));
                    r
                };
                assert_eq!(bf_r, rug_r, "exp({a}) at p={p}, mode={mode:?}");
            }
        }
    }
}
