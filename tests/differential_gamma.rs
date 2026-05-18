//! MPFR differential: `BigFloat::gamma` matches `rug::Float::gamma`.
//!
//! Capped at precisions ≤ 512 bits per the asymptotic `z_min`
//! memory: pfloat's Stirling implementation uses a 17-pair
//! Bernoulli table that caps the target precision around ~600 bits.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, rug_from_i64, sweep_size,
    ALL_ROUNDING_MODES,
};

#[test]
fn gamma_matches_mpfr_on_small_positive_integers() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6e");
    let cases = sweep_size().min(200);

    // Cap at p ≤ 256 to stay well below the 600-bit Stirling
    // ceiling; cap a ≤ 15 so gamma(a) stays in i64 dynamic range
    // and the asymptotic shift is modest.
    let precisions: &[u32] = &[53, 113, 256];
    for &p in precisions {
        for _ in 0..cases {
            let a = next_i64_in(&mut state, 1, 15);
            for &mode in ALL_ROUNDING_MODES {
                let bf_r = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let (r, _status) = a_bf.gamma(mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = {
                    let a_rg = rug_from_i64(a, p);
                    let (r, _ord) =
                        rug::Float::with_val_round(p, a_rg.gamma_ref(), mpfr_round_of(mode));
                    r
                };
                assert_eq!(bf_r, rug_r, "gamma({a}) at p={p}, mode={mode:?}");
            }
        }
    }
}
