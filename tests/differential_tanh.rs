//! MPFR differential: `BigFloat::tanh` matches `rug::Float::tanh`
//! bit-for-bit at every tested precision and IEEE rounding mode.
//!
//! Phase 4 (pf-3rtr.6) adds the standing five-mode differential `tanh`
//! lacked. The directed-mode bug-hunt lane (`directed_mode_bug_hunt.rs`)
//! carries the heavier coverage: a bounded all-mode control, the
//! saturation boundary-and-tail sweep (the ADR-0080 fix), and the tiny-x
//! short-circuit. This lane is the standing integer-and-fractional
//! regression against MPFR.
//!
//! `tanh` is odd and defined on the whole real line; both arms exercise
//! both signs. The integer arm stays in the Ziv-converging band (the
//! saturation tail toward ±1 is the bug-hunt sweep's job); the fractional
//! arm reaches tiny-x and off-grid hard-to-round inputs.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_oracle_for_mode, next_i64_in, rug_from_i64,
    sweep_size, BIT_EXACT_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};
use pfloat::RoundingMode;

#[test]
fn tanh_matches_mpfr_on_integer_and_fractional_inputs() {
    let mut state: u64 = u64::from_le_bytes(*b"pf3rtanh");
    let cases = sweep_size().min(1_000);

    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            // Integer arm: moderate magnitude (both signs).
            let a = next_i64_in(&mut state, -256, 256);
            // Fractional arm: a/b, off-grid, reaching tiny-x for large b.
            let num = next_i64_in(&mut state, -(1 << 20), 1 << 20);
            let den = next_i64_in(&mut state, 1, 1 << 20);
            let x_frac = {
                let n = bigfloat_from_i64(num, p);
                let d = bigfloat_from_i64(den, p);
                n.div(&d, RoundingMode::NearestEven).0
            };
            let x_frac_rug = bigfloat_to_rug(&x_frac);

            for &mode in BIT_EXACT_ROUNDING_MODES {
                let bf_i = bigfloat_to_rug(&bigfloat_from_i64(a, p).tanh(mode).0);
                let rug_i = mpfr_oracle_for_mode(
                    |prec, round| {
                        rug::Float::with_val_round(prec, rug_from_i64(a, prec).tanh_ref(), round).0
                    },
                    mode,
                    p,
                );
                assert_eq!(bf_i, rug_i, "tanh({a}) at p={p}, mode={mode:?}");

                let bf_f = bigfloat_to_rug(&x_frac.tanh(mode).0);
                let rug_f = mpfr_oracle_for_mode(
                    |prec, round| rug::Float::with_val_round(prec, x_frac_rug.tanh_ref(), round).0,
                    mode,
                    p,
                );
                assert_eq!(bf_f, rug_f, "tanh({num}/{den}) at p={p}, mode={mode:?}");
            }
        }
    }
}
