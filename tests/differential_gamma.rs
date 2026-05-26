//! MPFR differential: `BigFloat::gamma` matches `rug::Float::gamma`
//! bit-for-bit under every IEEE rounding mode (gamma is correctly
//! rounded via Ziv per slice p1.29, ADR-0038).
//!
//! Capped at precisions ≤ 256 bits per the lgamma Stirling memory:
//! pfloat's Stirling+Spouge dispatch (pf-l6s5) keeps lgamma
//! correctly rounded at all working precisions, but
//! `TRANSCENDENTAL_PRECISIONS` exercises the regime that p1.29
//! actually migrates here.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_oracle_for_mode, next_i64_in, rug_from_i64,
    sweep_size, BIT_EXACT_ROUNDING_MODES,
};

#[test]
fn gamma_matches_mpfr_on_small_positive_integers() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6e");
    let cases = sweep_size().min(200);

    // Cap a ≤ 15 so gamma(a) stays in i64 dynamic range and the
    // asymptotic shift is modest.
    let precisions: &[u32] = &[53, 113, 256];
    for &p in precisions {
        for _ in 0..cases {
            let a = next_i64_in(&mut state, 1, 15);
            for &mode in BIT_EXACT_ROUNDING_MODES {
                let bf_r = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let (r, _status) = a_bf.gamma(mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = mpfr_oracle_for_mode(
                    |prec, round| {
                        let a_rg = rug_from_i64(a, prec);
                        rug::Float::with_val_round(prec, a_rg.gamma_ref(), round).0
                    },
                    mode,
                    p,
                );
                assert_eq!(bf_r, rug_r, "gamma({a}) at p={p}, mode={mode:?}");
            }
        }
    }
}
