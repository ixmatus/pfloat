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

use super::{ln_2_at, ln_2pi_at, pi_at};

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

/// Computes `ψ(z) = digamma(z)` for `z` large enough that the
/// digamma asymptotic
///
/// ```text
/// ψ(z) = ln(z) − 1/(2z) − Σ_{k=1}^N B_{2k} / (2k · z^(2k))
/// ```
///
/// converges to `working_prec` bits. The coefficients reuse the
/// hardcoded `STIRLING_C` table: `B_{2k}/(2k) = c_k · (2k − 1)`.
/// Same truncation cap (17 terms) as `stirling_lgamma`; callers
/// shift small `z` up before invoking.
pub(super) fn stirling_digamma(z: &BigFloat, working_prec: u32) -> BigFloat {
    let z_w = z
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let two = BigFloat::try_from_i64_exact(2, working_prec).expect("precision >= 1");

    // Leading: ln(z) − 1/(2z).
    let (ln_z, _) = z_w.ln(RoundingMode::NearestEven);
    let (two_z, _) = two.mul(&z_w, RoundingMode::NearestEven);
    let (half_over_z, _) = one.div(&two_z, RoundingMode::NearestEven);
    let (with_correction, _) = ln_z.sub(&half_over_z, RoundingMode::NearestEven);

    // Tail: −Σ c'_k / z^(2k) where c'_k = B_{2k}/(2k) = c_k·(2k−1).
    let (z_sq, _) = z_w.mul(&z_w, RoundingMode::NearestEven);
    let (one_over_z_sq, _) = one.div(&z_sq, RoundingMode::NearestEven);

    let mut z_power = one_over_z_sq.clone(); // z^(-2)
    let mut sum = with_correction;

    let mut prev_term_exp: i64 = 1;
    for (k_minus_1, &(num, den)) in STIRLING_C.iter().enumerate() {
        let k = (k_minus_1 + 1) as i64;
        let multiplier = 2 * k - 1; // (2k − 1)
        let num_scaled = num.saturating_mul(multiplier);
        let num_bf =
            BigFloat::try_from_i64_exact(num_scaled, working_prec).expect("precision >= 1");
        let den_bf = BigFloat::try_from_i64_exact(den, working_prec).expect("precision >= 1");
        let (coef, _) = num_bf.div(&den_bf, RoundingMode::NearestEven);
        let (term, _) = coef.mul(&z_power, RoundingMode::NearestEven);

        let term_exp = match &term.class {
            crate::class::Class::Normal { exponent, .. } => *exponent,
            _ => -i64::from(working_prec) - 1,
        };
        if term_exp >= prev_term_exp {
            break;
        }
        prev_term_exp = term_exp;
        let (next_sum, _) = sum.sub(&term, RoundingMode::NearestEven);
        sum = next_sum;
        if term_exp < -i64::from(working_prec) - 4 {
            break;
        }

        // Advance: z^(-2(k+1)) = z^(-2k) · z^(-2).
        let (next_power, _) = z_power.mul(&one_over_z_sq, RoundingMode::NearestEven);
        z_power = next_power;
    }

    sum
}

/// Computes `ln Γ(z)` via Spouge's approximation
/// (Spouge, J.L. "Computation of the Gamma, Digamma, and Trigamma
/// Functions" SIAM J. Numer. Anal. 31:1, 1994). For positive real `z`
/// and a chosen integer parameter `a > 0`:
///
/// ```text
/// Γ(z+1) = (z+a)^(z+1/2) · e^(−(z+a)) · S(z, a)
///   S(z, a) = √(2π) + Σ_{k=1}^{a−1} c_k / (z+k)
///   c_k  = (−1)^(k−1) / (k−1)! · (a−k)^(k−1/2) · e^(a−k)
/// ```
///
/// and `ln Γ(z) = ln Γ(z+1) − ln(z)`.
///
/// Unlike Stirling's asymptotic series (which requires the caller
/// to shift `z` up beyond `z_min ≈ target_precision/log_2(z_min)`
/// for the truncation error to fall below ULP), Spouge's
/// approximation works directly for any positive `z` with cost
/// linear in the parameter `a`. For binary precision `p`, the
/// truncation bound `|ε| ≤ a^(1/2 − a)` gives `2^(−p)` accuracy
/// when `a · log_2(a) ≥ p`. [`spouge_a_for`] selects a conservative
/// `a` with safety margin.
///
/// The lgamma kernel routes to this function for `target_precision`
/// past the 17-Bernoulli-pair Stirling table's reach (~895 bits).
/// Below that target, [`stirling_lgamma`] with upward-shift remains
/// the faster path. Phase 1f slice closes pf-l6s5 by dispatching at
/// the boundary.
///
/// References:
/// - Spouge, J.L. (1994), op. cit.
/// - Pugh, G.R. "An Analysis of the Lanczos Gamma Approximation"
///   `PhD` thesis, UBC (2004), §3 (error analysis for Spouge).
/// - Toth, V.T. "Programmable Calculators: The Gamma Function"
///   (2005), reference implementation pattern.
#[allow(dead_code)]
pub(super) fn spouge_lgamma(z: &BigFloat, working_prec: u32) -> BigFloat {
    let a = spouge_a_for(working_prec);
    let z_w = z
        .round_to_precision(working_prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0;

    let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let two = BigFloat::try_from_i64_exact(2, working_prec).expect("precision >= 1");
    let (half, _) = one.div(&two, RoundingMode::NearestEven);

    // √(2π) at working precision: sqrt(2 · π).
    let pi = pi_at(working_prec);
    let (two_pi, _) = two.mul(&pi, RoundingMode::NearestEven);
    let (sqrt_2pi, _) = two_pi.sqrt(RoundingMode::NearestEven);

    // S(z, a) = √(2π) + Σ_{k=1}^{a−1} c_k / (z+k). The c_k values
    // depend only on (working_prec, a) and are memoized per
    // working_prec; the per-call cost reduces to the partial-sum
    // accumulation plus the leading-factor evaluation below.
    let coefficients = spouge_coefficients(working_prec, a);
    let mut sum = sqrt_2pi;
    for (k_minus_1, c_k) in coefficients.iter().enumerate() {
        let k = (k_minus_1 + 1) as u32;
        let k_bf =
            BigFloat::try_from_i64_exact(i64::from(k), working_prec).expect("precision >= 1");
        let (z_plus_k, _) = z_w.add(&k_bf, RoundingMode::NearestEven);
        let (term, _) = c_k.div(&z_plus_k, RoundingMode::NearestEven);
        let (next_sum, _) = sum.add(&term, RoundingMode::NearestEven);
        sum = next_sum;
    }

    // Leading factor in log: (z + 1/2) · ln(z + a) − (z + a).
    let a_bf = BigFloat::try_from_i64_exact(i64::from(a), working_prec).expect("precision >= 1");
    let (z_plus_a, _) = z_w.add(&a_bf, RoundingMode::NearestEven);
    let (z_plus_half, _) = z_w.add(&half, RoundingMode::NearestEven);
    let (ln_z_plus_a, _) = z_plus_a.ln(RoundingMode::NearestEven);
    let (term1, _) = z_plus_half.mul(&ln_z_plus_a, RoundingMode::NearestEven);
    let (leading, _) = term1.sub(&z_plus_a, RoundingMode::NearestEven);

    // ln Γ(z+1) = leading + ln(S(z, a)).
    let (ln_sum, _) = sum.ln(RoundingMode::NearestEven);
    let (ln_gamma_z_plus_1, _) = leading.add(&ln_sum, RoundingMode::NearestEven);

    // ln Γ(z) = ln Γ(z+1) − ln(z).
    let (ln_z, _) = z_w.ln(RoundingMode::NearestEven);
    let (ln_gamma_z, _) = ln_gamma_z_plus_1.sub(&ln_z, RoundingMode::NearestEven);

    ln_gamma_z
}

/// Compute the Spouge coefficient vector `[c_1, c_2, …, c_{a-1}]`
/// at the supplied working precision, memoized per `working_prec`
/// (since `a = spouge_a_for(working_prec)` is a function of
/// `working_prec`, the cache key reduces to one dimension).
///
/// Each `c_k` is computed via the log form to avoid the
/// fractional-exponent pow that the textbook
/// `(a−k)^(k−1/2) · e^(a−k)` form would require: one ln call for
/// `ln(a−k)`, one exp call for `c_k = ±exp(ln|c_k|)`, plus a few
/// mul/sub. Total cost is O(a) ln + O(a) exp calls, dominated by
/// the ln/exp at `working_prec`. Without memoization the
/// `spouge_lgamma` runtime at p=1024 would dominate the
/// `differential_zeta` lane's wall-clock; the cache amortizes
/// across the lane's 15 dyadic inputs to a single computation per
/// `working_prec`.
#[allow(dead_code)]
fn spouge_coefficients(working_prec: u32, a: u32) -> alloc::vec::Vec<BigFloat> {
    spouge_cache::memoized(working_prec, || {
        spouge_coefficients_compute(working_prec, a)
    })
}

#[allow(dead_code)]
fn spouge_coefficients_compute(working_prec: u32, a: u32) -> alloc::vec::Vec<BigFloat> {
    let two = BigFloat::try_from_i64_exact(2, working_prec).expect("precision >= 1");
    let mut coefficients = alloc::vec::Vec::with_capacity((a - 1) as usize);
    let mut ln_factorial =
        BigFloat::try_new_zero(crate::sign::Sign::Positive, working_prec).expect("precision >= 1"); // ln(0!) = 0.
    for k in 1u32..a {
        let a_minus_k =
            BigFloat::try_from_i64_exact(i64::from(a - k), working_prec).expect("precision >= 1");
        let k_minus_half_int = BigFloat::try_from_i64_exact(i64::from(2 * k) - 1, working_prec)
            .expect("precision >= 1");
        let (k_minus_half, _) = k_minus_half_int.div(&two, RoundingMode::NearestEven);

        // ln|c_k| = (k − 1/2)·ln(a − k) + (a − k) − ln((k − 1)!).
        let (ln_a_minus_k, _) = a_minus_k.ln(RoundingMode::NearestEven);
        let (term_a, _) = k_minus_half.mul(&ln_a_minus_k, RoundingMode::NearestEven);
        let (term_b, _) = term_a.add(&a_minus_k, RoundingMode::NearestEven);
        let (ln_c_k_abs, _) = term_b.sub(&ln_factorial, RoundingMode::NearestEven);

        // |c_k| = exp(ln|c_k|), sign = (−1)^(k−1).
        let (c_k_abs, _) = ln_c_k_abs.exp(RoundingMode::NearestEven);
        let c_k = if k % 2 == 0 {
            c_k_abs.negated()
        } else {
            c_k_abs
        };
        coefficients.push(c_k);

        // Advance: ln(k!) = ln((k − 1)!) + ln(k).
        let k_bf =
            BigFloat::try_from_i64_exact(i64::from(k), working_prec).expect("precision >= 1");
        let (ln_k, _) = k_bf.ln(RoundingMode::NearestEven);
        let (next_ln_fact, _) = ln_factorial.add(&ln_k, RoundingMode::NearestEven);
        ln_factorial = next_ln_fact;
    }
    coefficients
}

/// Thread-local cache for Spouge coefficient vectors, keyed by
/// `working_prec`. Mirrors the `agm_constants::cache` pattern: a
/// small RefCell-protected list under `std`, transparent
/// passthrough under `no_std`. `CACHE_CAP` bounds memory for
/// callers that sweep unboundedly many precisions.
#[cfg(feature = "std")]
mod spouge_cache {
    use super::BigFloat;
    use std::cell::RefCell;

    /// Soft bound on retained working precisions per thread. lgamma
    /// touches a few distinct working precisions (target+64 across
    /// the few targets a sweep exercises), so this never bites in
    /// normal use; it caps memory under sweeps over many precisions.
    const CACHE_CAP: usize = 16;

    std::thread_local! {
        static CACHE: RefCell<Vec<(u32, Vec<BigFloat>)>> =
            const { RefCell::new(Vec::new()) };
    }

    pub(super) fn memoized(
        working_prec: u32,
        compute: impl FnOnce() -> Vec<BigFloat>,
    ) -> Vec<BigFloat> {
        if let Some(hit) = CACHE.with(|c| {
            c.borrow()
                .iter()
                .find(|(p, _)| *p == working_prec)
                .map(|(_, v)| v.clone())
        }) {
            return hit;
        }
        let value = compute();
        CACHE.with(|c| {
            let mut entries = c.borrow_mut();
            if entries.len() < CACHE_CAP {
                entries.push((working_prec, value.clone()));
            }
        });
        value
    }
}

#[cfg(not(feature = "std"))]
mod spouge_cache {
    use super::BigFloat;
    use alloc::vec::Vec;

    pub(super) fn memoized(
        _working_prec: u32,
        compute: impl FnOnce() -> Vec<BigFloat>,
    ) -> Vec<BigFloat> {
        compute()
    }
}

/// Pick the Spouge parameter `a` to deliver `working_prec` bits of
/// accuracy with safety margin.
///
/// Spouge's truncation bound `|ε| ≤ a^(1/2 − a)` requires
/// `(a − 1/2) · log_2(a) ≥ working_prec` for `2^(−working_prec)`
/// relative error. For moderately large `a` this is well
/// approximated by `a · log_2(a) ≥ working_prec`. The function
/// `a · log_2(a)` is monotone increasing; this helper picks `a`
/// with explicit margin to absorb cancellation in the partial sum
/// and the leading-factor logarithms.
///
/// Empirical formula: `a = max(20, ceil(working_prec / 5) + 20)`.
/// For `working_prec = 1024` this gives `a = 225`
/// (`a·log_2(a) ≈ 1759`, margin > 700 bits). For
/// `working_prec = 4096` this gives `a = 840` (`a·log_2(a) ≈ 8198`,
/// margin > 4000 bits). The margin is asymptotically wasteful but
/// cost is linear in `a`; this trades CPU for confidence in the
/// bit-exactness gate.
#[allow(dead_code)]
pub(super) fn spouge_a_for(working_prec: u32) -> u32 {
    (working_prec / 5).saturating_add(20).max(20)
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
