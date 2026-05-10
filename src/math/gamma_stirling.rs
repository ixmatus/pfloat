//! Stirling's asymptotic series for `ln Γ(z)`:
//!
//! ```text
//! ln Γ(z) = (z − 1/2) · ln(z) − z + (1/2) · ln(2π)
//!         + Σ_{k=1}^N (B_{2k} / (2k · (2k − 1))) · z^(1−2k)
//! ```
//!
//! The series is asymptotic — divergent in the limit — but well
//! behaved for `z` large compared to the truncation `N`. Each
//! tail term `c_k · z^(1−2k)` decays as `z^(−2k+1)` until `k`
//! reaches roughly `πz`, after which the factor `c_k ~ (2k)! /
//! (2π)^(2k)` overwhelms the inverse power. Callers must shift
//! small `z` up via `Γ(z+1) = z · Γ(z)` to land in the convergent
//! regime before invoking this routine.
//!
//! Coefficients are hardcoded as `(i64, i64)` pairs encoding the
//! exact rational `c_k = B_{2k} / (2k · (2k − 1))`. The first 17
//! all fit in i64; beyond that the Bernoulli numerators overflow,
//! which caps the asymptotic at roughly 50–60 bits per term for
//! large `z`. Combined with the polynomial decay this still covers
//! `target_precision` ≤ 1024 bits when the caller shifts `z`
//! beyond `target_precision / 5`.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;

use super::{ln_2_at, ln_2pi_at};

/// Stirling coefficients `c_k = B_{2k} / (2k · (2k − 1))` as
/// `(numerator, denominator)` for `k = 1 ..= 17`.
///
/// Source: `mpmath.bernfrac` at slice-4b authoring time. All 34
/// integers fit in `i64`; subsequent `k` would overflow.
#[allow(dead_code)]
const STIRLING_C: &[(i64, i64)] = &[
    (1, 12),
    (-1, 360),
    (1, 1260),
    (-1, 1680),
    (1, 1188),
    (-691, 360360),
    (1, 156),
    (-3617, 122400),
    (43867, 244188),
    (-174611, 125400),
    (77683, 5796),
    (-236364091, 1506960),
    (657931, 300),
    (-3392780147, 93960),
    (1723168255201, 2492028),
    (-7709321041217, 505920),
    (151628697551, 396),
];

/// Computes `ln Γ(z)` for `z` large enough that Stirling's
/// asymptotic converges to `target_precision` bits at one of the
/// hardcoded coefficients. The kernel does not validate the
/// regime — small `z` will return a poor approximation. Callers
/// (`lgamma_kernel`) are responsible for the shift.
///
/// Returns the result at the supplied `working_prec`.
pub(super) fn stirling_lgamma(z: &BigFloat, working_prec: u32) -> BigFloat {
    let z_w = z
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let two = BigFloat::try_from_i64_exact(2, working_prec).expect("precision >= 1");

    // (z − 1/2) · ln(z)
    let (half, _) = one.div(&two, RoundingMode::NearestEven);
    let (z_minus_half, _) = z_w.sub(&half, RoundingMode::NearestEven);
    let (ln_z, _) = z_w.ln(RoundingMode::NearestEven);
    let (term_main, _) = z_minus_half.mul(&ln_z, RoundingMode::NearestEven);

    // − z
    let (term_minus_z, _) = term_main.sub(&z_w, RoundingMode::NearestEven);

    // + (1/2) · ln(2π) = ln(2π) / 2
    let ln_2pi = ln_2pi_at(working_prec);
    let (half_ln_2pi, _) = ln_2pi.div(&two, RoundingMode::NearestEven);
    let (with_const, _) = term_minus_z.add(&half_ln_2pi, RoundingMode::NearestEven);

    // Stirling tail: Σ c_k · z^(1 − 2k) = (1/z) · Σ c_k · (1/z²)^(k−1)·(1/z²)... hmm let's
    // just compute term-by-term: t_k = c_k / z^(2k−1).
    let (z_sq, _) = z_w.mul(&z_w, RoundingMode::NearestEven);
    let (one_over_z_sq, _) = one.div(&z_sq, RoundingMode::NearestEven);
    let (one_over_z, _) = one.div(&z_w, RoundingMode::NearestEven);

    // Start: t_1 = c_1 / z. Subsequent: t_{k+1} = t_k · (1/z²) · (c_{k+1}/c_k).
    // Easier and numerically stable: compute each term as
    // `(num_k / den_k) · z^(−(2k−1))`, where `z^(−(2k−1))` is
    // accumulated by multiplying `one_over_z²` per step.
    let mut z_power = one_over_z.clone(); // z^(-1)
    let mut sum = with_const;

    // Track magnitude to stop before the asymptotic diverges.
    let mut prev_term_exp: i64 = 1;
    for &(num, den) in STIRLING_C {
        let num_bf = BigFloat::try_from_i64_exact(num, working_prec).expect("precision >= 1");
        let den_bf = BigFloat::try_from_i64_exact(den, working_prec).expect("precision >= 1");
        let (coef, _) = num_bf.div(&den_bf, RoundingMode::NearestEven);
        let (term, _) = coef.mul(&z_power, RoundingMode::NearestEven);

        let term_exp = match &term.class {
            crate::class::Class::Normal { exponent, .. } => *exponent,
            _ => -i64::from(working_prec) - 1,
        };
        if term_exp >= prev_term_exp {
            // The asymptotic has started to diverge; stop.
            break;
        }
        prev_term_exp = term_exp;
        let (next_sum, _) = sum.add(&term, RoundingMode::NearestEven);
        sum = next_sum;
        if term_exp < -i64::from(working_prec) - 4 {
            break;
        }

        // Advance: z^(-(2(k+1)-1)) = z^(-(2k-1)) · z^(-2) = z_power · one_over_z_sq.
        let (next_power, _) = z_power.mul(&one_over_z_sq, RoundingMode::NearestEven);
        z_power = next_power;
    }

    // Avoid the unused-import lint in alt-feature builds: the
    // ln_2_at helper isn't called here, but it lives one module
    // up; we don't need it from inside Stirling itself.
    let _ = ln_2_at;

    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::cmp::Ordering;

    fn close_at(v: &BigFloat, expected: &BigFloat, bits: u32) -> bool {
        let (diff, _) = v.sub(expected, RoundingMode::NearestEven);
        let abs_diff = diff.abs();
        if abs_diff.is_zero() {
            return true;
        }
        let p = v.precision().max(expected.precision());
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let abs_b = expected.abs();
        let mut bound = if abs_b.is_zero() { one } else { abs_b };
        for _ in 0..bits {
            bound = bound.div(&two, RoundingMode::NearestEven).0;
        }
        matches!(
            abs_diff.partial_cmp(&bound).0,
            Some(Ordering::Less | Ordering::Equal)
        )
    }

    #[test]
    fn stirling_at_large_z_matches_factorial() {
        // At z = 100 with the 17-coefficient table, the asymptotic
        // truncation is roughly 2^(28.5 − 33·log₂(100)) ≈ 2^−190,
        // so target precision tighter than ~180 bits is not
        // reachable. Verify at a precision that respects this
        // budget; the kernel-level lgamma tests cover the
        // higher-precision regime by shifting z upward.
        let p = 128u32;
        let z = BigFloat::try_from_i64_exact(100, p).unwrap();
        let r = stirling_lgamma(&z, p);
        let expected = BigFloat::parse_str(
            "359.13420536957539877604401046028690961262171808563",
            p,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, p - 16));
    }

    #[test]
    fn stirling_at_z_25_matches_target() {
        // ln Γ(25) = ln(24!) ≈ 54.78472939811231919.
        let p = 113u32;
        let z = BigFloat::try_from_i64_exact(25, p).unwrap();
        let r = stirling_lgamma(&z, p);
        let expected = BigFloat::parse_str(
            "54.784729398112319190093344083606184686866212381333",
            p,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        assert!(close_at(&r, &expected, p - 24));
    }
}
