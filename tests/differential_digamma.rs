//! MPFR differential: `BigFloat::digamma` matches `rug::Float::digamma`.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, rug_from_i64, sweep_size,
    ALL_ROUNDING_MODES,
};

#[test]
fn digamma_matches_mpfr_on_positive_integers() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6e");
    let cases = sweep_size().min(200);

    let precisions: &[u32] = &[53, 113, 256];
    for &p in precisions {
        for _ in 0..cases {
            let a = next_i64_in(&mut state, 1, 1000);
            for &mode in ALL_ROUNDING_MODES {
                let bf_r = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let (r, _status) = a_bf.digamma(mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = {
                    let a_rg = rug_from_i64(a, p);
                    let (r, _ord) =
                        rug::Float::with_val_round(p, a_rg.digamma_ref(), mpfr_round_of(mode));
                    r
                };
                assert_eq!(bf_r, rug_r, "digamma({a}) at p={p}, mode={mode:?}");
            }
        }
    }
}
