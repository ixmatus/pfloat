//! MPFR differential: `BigFloat::cbrt` matches `rug::Float::cbrt`
//! bit-for-bit at every tested precision and IEEE rounding mode.
//!
//! cbrt is an exact-integer root (ADR-0056): a perfect cube is exact and
//! every other integer's cube root is irrational, so it can never land on
//! a half-way tie, and the result is correctly rounded under every mode by
//! construction. Unlike `sqrt`, cbrt is defined for negatives, so this
//! lane sweeps **signed** integer inputs — the negative branch (the real
//! cube root, sign carried through `round_finite_to_precision`) is the
//! novel part rug pins down.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_oracle_for_mode, next_i64_in, rug_from_i64,
    sweep_size, BIT_EXACT_ROUNDING_MODES, SWEEP_PRECISIONS,
};

#[test]
fn cbrt_matches_mpfr_on_signed_i64() {
    let mut state: u64 = u64::from_le_bytes(*b"pflm1bcb");
    let cases = sweep_size();

    for &p in SWEEP_PRECISIONS {
        let cap = if p >= 64 {
            i64::MAX
        } else {
            (1_i64 << (p as i64 - 1)) - 1
        };
        for _ in 0..cases {
            let a = next_i64_in(&mut state, -cap, cap);
            for &mode in BIT_EXACT_ROUNDING_MODES {
                let bf_root = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let (root, _status) = a_bf.cbrt(mode);
                    bigfloat_to_rug(&root)
                };
                let rug_root = mpfr_oracle_for_mode(
                    |prec, round| {
                        let a_rg = rug_from_i64(a, prec);
                        rug::Float::with_val_round(prec, a_rg.cbrt_ref(), round).0
                    },
                    mode,
                    p,
                );
                assert_eq!(bf_root, rug_root, "cbrt({a}) at p={p}, mode={mode:?}");
            }
        }
    }
}
