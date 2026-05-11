//! MPFR differential: `BigFloat::pow` against `rug::Float::pow`
//! with **2 ULP tolerance**.
//!
//! Slice 3c's pow ships `exp(y · ln(x))` at working precision and
//! rounds back to target. The composition accumulates rounding
//! from ln, the multiplication, and exp, so the result is not
//! correctly-rounded — a 1 ULP difference from MPFR's
//! `mpfr_pow` is expected on integer-exponent inputs (where
//! MPFR has a fast path that avoids the exp/ln composition).
//! Closing the gap to bit-exact is a Phase 5 / Phase 7 follow-up
//! that needs either a Ziv-strategy retry or an integer-exponent
//! fast path in pfloat's pow.
//!
//! Restricted to positive bases and small finite exponents. The
//! full IEEE 754-2019 §9.2.1 table (zero base, infinity base,
//! negative base with integer / non-integer exponent) is covered
//! by the Kani harnesses in `src/verify/pow.rs`.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, rug_from_i64, sweep_size,
    ALL_ROUNDING_MODES, TRANSCENDENTAL_PRECISIONS,
};

#[test]
fn pow_matches_mpfr_on_positive_base_small_exponent() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6c");
    let cases = sweep_size().min(500);

    for &p in TRANSCENDENTAL_PRECISIONS {
        for _ in 0..cases {
            // Base in [1, 100], exponent in [-10, 10] keeps the
            // result inside i64 range without testing the
            // overflow / underflow paths.
            let base = next_i64_in(&mut state, 1, 100);
            let exp = next_i64_in(&mut state, -10, 10);
            for &mode in ALL_ROUNDING_MODES {
                let bf_r = {
                    let b_bf = bigfloat_from_i64(base, p);
                    let e_bf = bigfloat_from_i64(exp, p);
                    let (r, _status) = b_bf.pow(&e_bf, mode);
                    bigfloat_to_rug(&r)
                };
                let rug_r = {
                    use rug::ops::Pow;
                    let b_rg = rug_from_i64(base, p);
                    let e_rg = rug_from_i64(exp, p);
                    let (r, _ord) =
                        rug::Float::with_val_round(p, Pow::pow(&b_rg, &e_rg), mpfr_round_of(mode));
                    r
                };
                // Compare with 2 ULP tolerance. ULP at the result
                // magnitude is `|r| * 2^-(p-1)`; we accept twice
                // that for the exp/ln rounding composition.
                let diff = rug::Float::with_val(p + 32, &bf_r - &rug_r).abs();
                let ulp_scale = 2.0_f64.powi(-(p as i32 - 2));
                let tol = rug::Float::with_val(p + 32, &rug_r).abs()
                    * rug::Float::with_val(p + 32, ulp_scale);
                assert!(
                    diff <= tol || rug_r.is_zero(),
                    "pow({base}, {exp}) at p={p}, mode={mode:?}: pfloat={bf_r}, mpfr={rug_r}, diff={diff}"
                );
            }
        }
    }
}
