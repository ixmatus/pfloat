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
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, rug_from_i64, sweep_size,
    ALL_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};

#[test]
fn exp_matches_mpfr_on_small_integer_inputs() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6c");
    let cases = sweep_size().min(1_000); // exp is expensive; cap CI sweep

    for &p in TRANSCENDENTAL_PRECISIONS {
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
