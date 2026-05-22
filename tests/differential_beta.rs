//! MPFR differential: `BigFloat::beta` against an `lgamma`-based
//! reference. MPFR ships no direct `beta(a, b)`, so the oracle is
//! `sign · exp(ln|Γ(a)| + ln|Γ(b)| − ln|Γ(a + b)|)` evaluated at a
//! higher working precision.
//!
//! The sign is taken from MPFR's own `mpfr_lgamma` sign output (rug
//! [`Float::ln_abs_gamma`], which returns
//! `Ordering::Less`/`Greater` for `Γ < 0`/`Γ > 0`). It is therefore
//! an *independent* oracle for the negative domain: it does not
//! reuse pfloat's `gamma_sign_of`, so the negative-domain extension
//! (ADR-0030, derived from DLMF 5.12.1 / 5.5.3 / 5.2) is cross-
//! checked against MPFR rather than against the same reflection rule
//! it implements. `Γ` has no zeros (DLMF 5.2), so on every non-pole
//! argument the sign is strictly `Less` or `Greater`.
//!
//! Two arms share the one oracle: the original positive-integer
//! smoke arm, and a negative-non-integer arm (ADR-0030 case 2:
//! `a, b, a+b` all off the Γ poles, finite signed result). Both keep
//! the 2-ULP slack — the oracle compounds three `lgamma` calls plus
//! an `exp`, and `beta` is not Ziv-corrected (the NE-only loose
//! tier, `feedback_differential_lane_cost`).

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use core::cmp::Ordering;

use differential::{
    bigfloat_from_i64, bigfloat_to_rug, mpfr_round_of, next_i64_in, rug_from_i64, sweep_size,
    ALL_ROUNDING_MODES,
};
use pfloat::RoundingMode;
use rug::Float;

/// `sign · exp(ln|Γ(a)| + ln|Γ(b)| − ln|Γ(a+b)|)`, the
/// sign-tracking `lgamma`-composition oracle for `β(a, b)`. `a` and
/// `b` are already at `p_oracle`; `a + b` is formed exactly there.
/// The sign is the product of MPFR's three `Γ`-sign Orderings.
fn beta_via_lgamma(a: &Float, b: &Float, p: u32, p_oracle: u32, mode: RoundingMode) -> Float {
    let (lg_a, sa) = a.clone().ln_abs_gamma();
    let (lg_b, sb) = b.clone().ln_abs_gamma();
    let ab = Float::with_val(p_oracle, a + b);
    let (lg_ab, sab) = ab.ln_abs_gamma();
    let ln_mag = lg_a + lg_b - lg_ab;
    let mag = Float::with_val(p_oracle, ln_mag.exp());
    let negatives = [sa, sb, sab]
        .iter()
        .filter(|o| **o == Ordering::Less)
        .count();
    let signed = if negatives % 2 == 1 { -mag } else { mag };
    Float::with_val_round(
        p,
        signed,
        mpfr_round_of(mode).expect("NE-only lane: NearestEven has an MPFR equivalent (pf-suo)"),
    )
    .0
}

/// `|x − y| ≤ 2 · 2^(−(p−2)) · |y|`: the loose 2-ULP relative
/// tolerance shared by both arms.
fn within_two_ulp(got: &Float, oracle: &Float, p: u32, p_oracle: u32) -> bool {
    let diff = Float::with_val(p_oracle, got - oracle).abs();
    let ulp_scale = 2.0_f64.powi(-(p as i32 - 2));
    let ulp = Float::with_val(p_oracle, oracle).abs() * Float::with_val(p_oracle, ulp_scale);
    diff <= ulp
}

#[test]
fn beta_matches_mpfr_lgamma_composition_loosely() {
    let mut state: u64 = u64::from_le_bytes(*b"pfloat6e");
    let cases = sweep_size().min(100);

    // p=113 only: the lgamma-composition oracle compounds rounding
    // at three lgamma calls plus exp, so it is a loose smoke gate
    // (finite positive result for valid integer inputs, agreement
    // within 2 ULP), not a bit-exact lane.
    let p: u32 = 113;
    let p_oracle = p + 64;
    for _ in 0..cases {
        let a = next_i64_in(&mut state, 1, 20);
        let b = next_i64_in(&mut state, 1, 20);
        for &mode in ALL_ROUNDING_MODES {
            let (bf_r, _status) = {
                let a_bf = bigfloat_from_i64(a, p);
                let b_bf = bigfloat_from_i64(b, p);
                a_bf.beta(&b_bf, mode)
            };
            let rug_r = {
                let a_rg = rug_from_i64(a, p_oracle);
                let b_rg = rug_from_i64(b, p_oracle);
                beta_via_lgamma(&a_rg, &b_rg, p, p_oracle, mode)
            };
            assert!(
                bf_r.is_finite(),
                "beta({a}, {b}) at p={p}: pfloat got {bf_r}"
            );
            assert!(bf_r.is_sign_positive(), "beta({a}, {b}): positive integers");
            let bf_as_rug = bigfloat_to_rug(&bf_r);
            assert!(
                within_two_ulp(&bf_as_rug, &rug_r, p, p_oracle),
                "beta({a}, {b}) at p={p}, mode={mode:?}: pfloat={bf_as_rug}, oracle={rug_r}"
            );
        }
    }
}

/// ADR-0030 case 2: `a` a negative non-integer, `b` a (positive or
/// negative) non-integer, `a + b` a non-integer — all off the Γ
/// poles, so `β` is finite and signed. `a = −(2i+1)/2`,
/// `b = ±(2j+1)/4`, so `a + b` has an odd numerator over 4 and is
/// never an integer (never a denominator pole, never case 5).
#[test]
fn beta_negative_non_integer_matches_signed_lgamma_oracle() {
    let mut state: u64 = u64::from_le_bytes(*b"pfbeta8a");
    let cases = sweep_size().min(80);
    let p: u32 = 113;
    let p_oracle = p + 64;
    let ne = RoundingMode::NearestEven;

    for _ in 0..cases {
        let i = next_i64_in(&mut state, 1, 8);
        let j = next_i64_in(&mut state, 0, 8);
        let b_negative = next_i64_in(&mut state, 0, 1) == 1;
        let a_num = -(2 * i + 1); // /2  → negative non-integer
        let b_num = if b_negative { -(2 * j + 1) } else { 2 * j + 1 }; // /4

        for &mode in ALL_ROUNDING_MODES {
            let (bf_r, status) = {
                let two = bigfloat_from_i64(2, p);
                let four = bigfloat_from_i64(4, p);
                let a_bf = bigfloat_from_i64(a_num, p).div(&two, ne).0;
                let b_bf = bigfloat_from_i64(b_num, p).div(&four, ne).0;
                a_bf.beta(&b_bf, mode)
            };
            let rug_r = {
                let a_rg = rug_from_i64(a_num, p_oracle) / rug_from_i64(2, p_oracle);
                let b_rg = rug_from_i64(b_num, p_oracle) / rug_from_i64(4, p_oracle);
                beta_via_lgamma(&a_rg, &b_rg, p, p_oracle, mode)
            };

            // Case 2 is always finite and never raises an exception.
            assert!(
                bf_r.is_finite() && !bf_r.is_nan(),
                "beta({a_num}/2, {b_num}/4): pfloat got {bf_r}"
            );
            assert!(
                !status.invalid() && !status.div_by_zero(),
                "beta({a_num}/2, {b_num}/4): unexpected status {status:?}"
            );
            // Independent sign cross-check against MPFR's lgamma
            // sign (this is what exercises gamma_sign_of).
            assert_eq!(
                bf_r.is_sign_negative(),
                rug_r.is_sign_negative(),
                "beta({a_num}/2, {b_num}/4) sign: pfloat={bf_r}, oracle={rug_r}"
            );
            let bf_as_rug = bigfloat_to_rug(&bf_r);
            assert!(
                within_two_ulp(&bf_as_rug, &rug_r, p, p_oracle),
                "beta({a_num}/2, {b_num}/4) at p={p}, mode={mode:?}: pfloat={bf_as_rug}, oracle={rug_r}"
            );
        }
    }
}
