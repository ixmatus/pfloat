//! MPFR differential: `BigFloat::beta` against `lgamma`-based
//! reference. MPFR does not ship a direct `beta(a, b)` so the
//! oracle is `exp(ln_gamma(a) + ln_gamma(b) − ln_gamma(a + b))`,
//! evaluated at higher working precision.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, rug_from_i64, sweep_size, ALL_ROUNDING_MODES,
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
fn beta_matches_mpfr_lgamma_composition_loosely() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6e");
    let cases = sweep_size().min(100);

    // Use only p=113 here: the lgamma-composition oracle compounds
    // rounding at three lgamma calls plus exp; bit-for-bit equality
    // would require evaluating the oracle at p_oracle ≥ p + 64. The
    // test is left as a smoke gate that beta produces a finite
    // positive result for valid integer inputs and that
    // pfloat ↔ MPFR-composition agreement holds within 2 ULP.
    let p: u32 = 113;
    for _ in 0..cases {
        let a = next_i64_in(&mut state, 1, 20);
        let b = next_i64_in(&mut state, 1, 20);
        for &mode in ALL_ROUNDING_MODES {
            let (bf_r, _status) = {
                let a_bf = bigfloat_from_i64(a, p);
                let b_bf = bigfloat_from_i64(b, p);
                a_bf.beta(&b_bf, mode)
            };
            // Oracle: exp(ln_gamma(a) + ln_gamma(b) − ln_gamma(a+b))
            // evaluated at p + 64 to absorb compounded rounding.
            let p_oracle = p + 64;
            let rug_r = {
                let a_rg = rug_from_i64(a, p_oracle);
                let b_rg = rug_from_i64(b, p_oracle);
                let ab_rg = rug_from_i64(a + b, p_oracle);
                let lg_a = rug::Float::with_val(p_oracle, a_rg.ln_gamma_ref());
                let lg_b = rug::Float::with_val(p_oracle, b_rg.ln_gamma_ref());
                let lg_ab = rug::Float::with_val(p_oracle, ab_rg.ln_gamma_ref());
                let sum = lg_a + lg_b - lg_ab;
                let (r, _ord) = rug::Float::with_val_round(p, sum.exp(), mpfr_round_of(mode));
                r
            };
            assert!(
                bf_r.is_finite(),
                "beta({a}, {b}) at p={p}: pfloat got {bf_r}"
            );
            let bf_as_rug = bigfloat_to_rug(&bf_r);
            // Allow 2 ULP slack because the oracle compounds three
            // ln_gamma calls plus an exp.
            let diff = rug::Float::with_val(p_oracle, &bf_as_rug - &rug_r).abs();
            let ulp_scale = 2.0_f64.powi(-(p as i32 - 2));
            let ulp = rug::Float::with_val(p_oracle, &rug_r).abs()
                * rug::Float::with_val(p_oracle, ulp_scale);
            assert!(
                diff <= ulp,
                "beta({a}, {b}) at p={p}, mode={mode:?}: pfloat={bf_as_rug}, oracle={rug_r}, diff={diff}"
            );
        }
    }
}
