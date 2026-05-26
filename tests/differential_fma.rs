//! MPFR differential: `BigFloat::fma` matches `rug::Float::mul_add`
//! bit-for-bit at every tested precision and IEEE rounding mode
//! (FMA is correctly rounded under every mode by construction).
//!
//! Phase 1f slice p1.23 widened this lane from NE-only to all five
//! IEEE modes via `mpfr_oracle_for_mode` (ADR-0038).
//!
//! Integer operands sized so `|a × b|` fits in `i64`. The rounding
//! to p bits happens inside both pfloat and rug.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_oracle_for_mode, next_i64_in, rug_from_i64,
    sweep_size, BIT_EXACT_ROUNDING_MODES, SWEEP_PRECISIONS,
};

#[test]
fn fma_matches_mpfr_on_i64_triples() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6b");
    let cases = sweep_size();

    for &p in SWEEP_PRECISIONS {
        // Cap `a` and `b` so |a * b| stays inside i64. `c` cap
        // matches add/sub.
        let ab_bits = 31_u32;
        let ab_cap = (1_i64 << ab_bits) - 1;
        let c_cap = if p >= 64 {
            i64::MAX
        } else {
            (1_i64 << (p as i64 - 2)) - 1
        };
        for _ in 0..cases {
            let a = next_i64_in(&mut state, -ab_cap, ab_cap);
            let b = next_i64_in(&mut state, -ab_cap, ab_cap);
            let c = next_i64_in(&mut state, -c_cap, c_cap);
            for &mode in BIT_EXACT_ROUNDING_MODES {
                let bf_r = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let b_bf = bigfloat_from_i64(b, p);
                    let c_bf = bigfloat_from_i64(c, p);
                    let (r, _status) = a_bf.fma(&b_bf, &c_bf, mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = mpfr_oracle_for_mode(
                    |prec, round| {
                        let a_rg = rug_from_i64(a, prec);
                        let b_rg = rug_from_i64(b, prec);
                        let c_rg = rug_from_i64(c, prec);
                        rug::Float::with_val_round(prec, a_rg.mul_add_ref(&b_rg, &c_rg), round).0
                    },
                    mode,
                    p,
                );
                assert_eq!(bf_r, rug_r, "fma({a}, {b}, {c}) at p={p}, mode={mode:?}");
            }
        }
    }
}
