//! MPFR differential: `BigFloat::log10` matches `rug::Float::log10`
//! bit-for-bit at every tested precision and IEEE rounding mode.
//!
//! Phase 4 (pf-3rtr.6) adds this lane after the `log10` kernel moved onto
//! the Ziv interval test (ADR-0081); the prior fixed-guard kernel had no
//! dedicated five-mode differential. The directed-mode bug-hunt lane
//! (`directed_mode_bug_hunt.rs`) covers the off-grid hard-to-round cases;
//! this lane is the standing integer-input regression against MPFR.
//!
//! `log10` domain is `x > 0`; non-positive inputs are covered by the Kani
//! harnesses in `src/verify/log10.rs`.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_oracle_for_mode, next_i64_in, rug_from_i64,
    sweep_size, BIT_EXACT_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};

#[test]
fn log10_matches_mpfr_on_positive_integer_inputs() {
    let mut state: u64 = u64::from_le_bytes(*b"pf3rlg10");
    let cases = sweep_size().min(1_000);

    for &p in TRANSCENDENTAL_PRECISIONS {
        let cap = if p >= 64 {
            i64::MAX
        } else {
            (1_i64 << (p as i64 - 1)) - 1
        };
        for _ in 0..cases {
            let a = next_i64_in(&mut state, 1, cap);
            for &mode in BIT_EXACT_ROUNDING_MODES {
                let bf_r = {
                    let a_bf = bigfloat_from_i64(a, p);
                    let (r, _status) = a_bf.log10(mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = mpfr_oracle_for_mode(
                    |prec, round| {
                        let a_rg = rug_from_i64(a, prec);
                        rug::Float::with_val_round(prec, a_rg.log10_ref(), round).0
                    },
                    mode,
                    p,
                );
                assert_eq!(bf_r, rug_r, "log10({a}) at p={p}, mode={mode:?}");
            }
        }
    }
}
