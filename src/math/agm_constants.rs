//! On-the-fly computation of the transcendental constants pfloat
//! consumes (`π`, `ln(2)`, derived `2/π`, `2/√π`, `ln(2π)`, `ln(10)`).
//!
//! Slice 7b lifts the hardcoded 1024-bit table ceiling that
//! slice-6h's MPFR differential lane documented as limitation #2.
//! For precisions at or below the table size, the
//! `super::pi_at`/`ln_2_at`/... dispatchers keep returning the
//! rounded table value. For higher precisions they call into this
//! module, which computes the constant at the requested precision
//! via AGM-based formulas.
//!
//! Algorithms:
//!
//! - **π** uses the Brent–Salamin iteration:
//!   ```text
//!   a_0 = 1, b_0 = 1/√2, t_0 = 1/4, p_0 = 1
//!   a_{n+1} = (a_n + b_n) / 2
//!   b_{n+1} = sqrt(a_n · b_n)
//!   t_{n+1} = t_n − p_n · (a_n − a_{n+1})²
//!   p_{n+1} = 2 · p_n
//!   ```
//!   Convergence is quadratic; `O(log p)` iterations suffice at
//!   working precision `p`. The result is
//!   `π ≈ (a_n + b_n)² / (4 · t_n)`.
//!
//! - **ln(2)** uses the atanh series identity
//!   `ln(2) = 2 · atanh(1/3) = 2 · (1/3 + 1/(3·27) + 1/(5·243) + …)`.
//!   The series converges roughly `log₂(9) ≈ 3.17` bits per term;
//!   `O(p)` terms suffice at working precision `p`.
//!
//! - **ln(10)** is bootstrapped from `ln(2)` via the identity
//!   `ln(10) = 3·ln(2) + 2·atanh(1/9)`, exploiting
//!   `ln(5/4) = 2·atanh(1/9)` and `ln(10) = ln(8) + ln(5/4)`. The
//!   atanh(1/9) series converges `log₂(81) ≈ 6.34` bits per term.
//!
//! Derived constants (`2/π`, `2/√π`, `ln(2π)`) compose `π` and
//! `ln(2)` via the standard arithmetic kernels at working precision.
//!
//! No special-case handling is required: every algorithm in this
//! module is invoked with `prec >= 1` and produces a finite normal
//! result by construction.
//!
//! ADR-0017 records the design.

use crate::big::BigFloat;
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

/// Internal-precision guard above the caller's requested precision.
/// The Brent–Salamin loop and the atanh series each accumulate at
/// most a few ULPs per iteration; 64 bits absorbs the compounded
/// error at the precisions pfloat supports.
const GUARD_BITS: u32 = 64;

/// Identifies which constant a cache entry holds, so distinct
/// constants requested at the same precision do not collide on the
/// `(kind, prec)` key.
#[allow(dead_code)] // some variants are unused under narrow feature combos
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Pi,
    Ln2,
    Ln10,
    TwoOverPi,
    TwoOverSqrtPi,
    Ln2Pi,
    EulerGamma,
}

/// Thread-local memoization of the AGM-computed constants.
///
/// Every transcendental kernel that needs `π`/`ln(2)`/… recomputes it
/// from scratch per call; the Brent–Salamin and atanh series are
/// `O(p)`–`O(p log p)` and dominate a high-precision differential
/// sweep (the `~10⁴ iterations` p=1024 lane was hour-scale). The
/// kernels request a small, stable set of working precisions, so a
/// per-thread `(kind, prec) → value` table collapses that cost to one
/// computation per distinct precision.
///
/// Keying by the *exact* requested precision is what makes this a
/// pure optimization: each `*_compute` returns the value already
/// rounded to `prec`, so a hit hands back the bit-identical `BigFloat`
/// the uncached path would have produced. There is no double rounding
/// because nothing is re-rounded on the hit path.
///
/// Under `no_std` there is no thread-local storage, so [`memoized`]
/// degrades to a transparent passthrough: identical results, only the
/// recompute cost returns.
#[cfg(feature = "std")]
mod cache {
    use super::{BigFloat, Kind};
    use std::cell::RefCell;

    /// Soft bound on retained `(kind, prec)` entries per thread. The
    /// kernels touch only a handful of distinct precisions, so this
    /// never bites in normal use; it caps memory under a caller that
    /// sweeps unboundedly many distinct precisions, degrading
    /// gracefully to recompute-on-miss rather than growing without
    /// limit (the relevant unbounded-input concern here).
    const CACHE_CAP: usize = 64;

    std::thread_local! {
        static CACHE: RefCell<Vec<(Kind, u32, BigFloat)>> =
            const { RefCell::new(Vec::new()) };
    }

    /// Returns the cached constant for `(kind, prec)`, or computes it
    /// via `compute` and inserts it on a miss. The lookup borrow is
    /// released before `compute` runs, so a recursive constant
    /// (`ln(10)` calling `ln(2)`, `2/π` calling `π`) cannot trip a
    /// `RefCell` double borrow.
    pub(super) fn memoized(kind: Kind, prec: u32, compute: impl FnOnce() -> BigFloat) -> BigFloat {
        if let Some(hit) = CACHE.with(|c| {
            c.borrow()
                .iter()
                .find(|(k, p, _)| *k == kind && *p == prec)
                .map(|(_, _, v)| v.clone())
        }) {
            return hit;
        }
        let value = compute();
        CACHE.with(|c| {
            let mut entries = c.borrow_mut();
            if entries.len() < CACHE_CAP {
                entries.push((kind, prec, value.clone()));
            }
        });
        value
    }

    #[cfg(test)]
    pub(super) fn entry_count() -> usize {
        CACHE.with(|c| c.borrow().len())
    }

    #[cfg(test)]
    pub(super) fn reset() {
        CACHE.with(|c| c.borrow_mut().clear());
    }
}

#[cfg(not(feature = "std"))]
mod cache {
    use super::{BigFloat, Kind};

    /// No thread-local storage without `std`: transparent
    /// passthrough. Correctness is identical to the cached path; only
    /// the per-call recompute cost differs.
    #[inline]
    pub(super) fn memoized(
        _kind: Kind,
        _prec: u32,
        compute: impl FnOnce() -> BigFloat,
    ) -> BigFloat {
        compute()
    }
}

/// Compute `π` at the requested precision via Brent–Salamin.
#[allow(dead_code)]
pub(super) fn pi_via_agm(prec: u32) -> BigFloat {
    cache::memoized(Kind::Pi, prec, || pi_compute(prec))
}

#[allow(dead_code)]
fn pi_compute(prec: u32) -> BigFloat {
    let working = prec.saturating_add(GUARD_BITS);

    let one = BigFloat::try_from_i64_exact(1, working).expect("precision >= 1");
    let two = BigFloat::try_from_i64_exact(2, working).expect("precision >= 1");
    let four = BigFloat::try_from_i64_exact(4, working).expect("precision >= 1");

    // a_0 = 1; b_0 = 1/√2; t_0 = 1/4; p_0 = 1.
    let mut a = one.clone();
    let (sqrt_two, _) = two.sqrt(RoundingMode::NearestEven);
    let (mut b, _) = one.div(&sqrt_two, RoundingMode::NearestEven);
    let (mut t, _) = one.div(&four, RoundingMode::NearestEven);
    let mut p = one.clone();

    // Quadratic convergence: bit agreement doubles each step.
    let max_iter = 64u32;
    let convergence_floor = -i64::from(working) - 4;

    for _ in 0..max_iter {
        let (diff, _) = a.sub(&b, RoundingMode::NearestEven);
        let abs_diff = diff.abs();
        let converged = match &abs_diff.class {
            Class::Zero { .. } => true,
            Class::Normal { exponent, .. } => *exponent < convergence_floor,
            _ => false,
        };
        if converged {
            break;
        }

        let (sum, _) = a.add(&b, RoundingMode::NearestEven);
        let (a_next, _) = sum.div(&two, RoundingMode::NearestEven);
        let (prod, _) = a.mul(&b, RoundingMode::NearestEven);
        let (b_next, _) = prod.sqrt(RoundingMode::NearestEven);
        let (c, _) = a.sub(&a_next, RoundingMode::NearestEven);
        let (c_sq, _) = c.mul(&c, RoundingMode::NearestEven);
        let (p_c_sq, _) = p.mul(&c_sq, RoundingMode::NearestEven);
        let (t_next, _) = t.sub(&p_c_sq, RoundingMode::NearestEven);
        let (p_next, _) = p.mul(&two, RoundingMode::NearestEven);

        a = a_next;
        b = b_next;
        t = t_next;
        p = p_next;
    }

    // π ≈ (a + b)² / (4 t).
    let (sum, _) = a.add(&b, RoundingMode::NearestEven);
    let (sum_sq, _) = sum.mul(&sum, RoundingMode::NearestEven);
    let (four_t, _) = four.mul(&t, RoundingMode::NearestEven);
    let (pi, _) = sum_sq.div(&four_t, RoundingMode::NearestEven);

    pi.round_to_precision(prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}

/// Compute `ln(2)` at the requested precision via the atanh series
/// `ln(2) = 2 · atanh(1/3)`.
pub(super) fn ln_2_via_atanh(prec: u32) -> BigFloat {
    cache::memoized(Kind::Ln2, prec, || ln_2_compute(prec))
}

fn ln_2_compute(prec: u32) -> BigFloat {
    let working = prec.saturating_add(GUARD_BITS);
    let atanh_third = atanh_one_over(3, working);
    let two = BigFloat::try_from_i64_exact(2, working).expect("precision >= 1");
    let (ln_2, _) = two.mul(&atanh_third, RoundingMode::NearestEven);
    ln_2.round_to_precision(prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}

/// Compute `ln(10)` at the requested precision via
/// `ln(10) = 3·ln(2) + 2·atanh(1/9)`.
pub(super) fn ln_10_via_atanh(prec: u32) -> BigFloat {
    cache::memoized(Kind::Ln10, prec, || ln_10_compute(prec))
}

fn ln_10_compute(prec: u32) -> BigFloat {
    let working = prec.saturating_add(GUARD_BITS);
    let ln_2 = ln_2_via_atanh(working);
    let three = BigFloat::try_from_i64_exact(3, working).expect("precision >= 1");
    let (three_ln_2, _) = three.mul(&ln_2, RoundingMode::NearestEven);
    let atanh_ninth = atanh_one_over(9, working);
    let two = BigFloat::try_from_i64_exact(2, working).expect("precision >= 1");
    let (two_atanh_ninth, _) = two.mul(&atanh_ninth, RoundingMode::NearestEven);
    let (result, _) = three_ln_2.add(&two_atanh_ninth, RoundingMode::NearestEven);
    result
        .round_to_precision(prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}

/// `atanh(1/n) = (1/n) + (1/n)³/3 + (1/n)⁵/5 + …` evaluated at
/// working precision for small positive integer `n`. Used by
/// [`ln_2_via_atanh`] (n = 3) and [`ln_10_via_atanh`] (n = 9).
fn atanh_one_over(n: i64, working_prec: u32) -> BigFloat {
    debug_assert!(n >= 2, "atanh(1/n) only valid for |n| >= 2");
    let one = BigFloat::try_from_i64_exact(1, working_prec).expect("precision >= 1");
    let n_big = BigFloat::try_from_i64_exact(n, working_prec).expect("precision >= 1");
    let (x, _) = one.div(&n_big, RoundingMode::NearestEven);
    let (x_squared, _) = x.mul(&x, RoundingMode::NearestEven);

    let mut sum = x.clone();
    let mut term = x;

    // Compute the maximum iteration count from the convergence
    // requirement: at iteration k, |term| ≈ x^(2k+1). For x = 1/n,
    // the magnitude is `n^(-(2k+1))`, and we want it below
    // `2^(-working_prec - 8)`. Solving for k gives
    // `k ≈ working_prec / (2 · log2(n)) + small_slack`. Cap by an
    // upper bound to keep this loop bounded for pathological inputs.
    let log2_n_times_8 = (63 - n.leading_zeros()) * 8; // floor(log2(n)) * 8
                                                       // saturating_mul: the `* 8` would otherwise overflow u32 at working
                                                       // precisions past ~5.4e8 bits (the sibling iteration-count helpers
                                                       // already saturate). Review 2026-05-29.
    let max_iter =
        (working_prec.saturating_add(64).saturating_mul(8) / log2_n_times_8.max(1)).max(64);
    for k in 1u32..=max_iter {
        let (next_term, _) = term.mul(&x_squared, RoundingMode::NearestEven);
        term = next_term;
        let divisor = BigFloat::try_from_i64_exact(i64::from(2 * k + 1), working_prec)
            .expect("precision >= 1");
        let (term_over_d, _) = term.div(&divisor, RoundingMode::NearestEven);
        let (next_sum, _) = sum.add(&term_over_d, RoundingMode::NearestEven);
        sum = next_sum;
    }

    sum
}

/// Compute `2/π` at the requested precision. Convenience wrapper:
/// `pi_via_agm(prec)` then `2/pi`. Returns a fresh `BigFloat` at
/// `prec` rounded under `NearestEven`.
#[allow(dead_code)]
pub(super) fn two_over_pi_via_agm(prec: u32) -> BigFloat {
    cache::memoized(Kind::TwoOverPi, prec, || two_over_pi_compute(prec))
}

#[allow(dead_code)]
fn two_over_pi_compute(prec: u32) -> BigFloat {
    let working = prec.saturating_add(GUARD_BITS);
    let pi = pi_via_agm(working);
    let two = BigFloat::try_from_i64_exact(2, working).expect("precision >= 1");
    let (result, _) = two.div(&pi, RoundingMode::NearestEven);
    result
        .round_to_precision(prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}

/// Compute `2/√π` at the requested precision via
/// `2 / sqrt(π)` at working precision.
// Consumed only by `two_over_sqrt_pi_at` (mod.rs), which is `specials`-
// gated; cfg-gated here so a trig-only / `big,agm` build is dead-code-free
// under clippy `-D warnings` (pf-gwum).
#[cfg(feature = "specials")]
pub(super) fn two_over_sqrt_pi_via_agm(prec: u32) -> BigFloat {
    cache::memoized(Kind::TwoOverSqrtPi, prec, || two_over_sqrt_pi_compute(prec))
}

#[cfg(feature = "specials")]
fn two_over_sqrt_pi_compute(prec: u32) -> BigFloat {
    let working = prec.saturating_add(GUARD_BITS);
    let pi = pi_via_agm(working);
    let (sqrt_pi, _) = pi.sqrt(RoundingMode::NearestEven);
    let two = BigFloat::try_from_i64_exact(2, working).expect("precision >= 1");
    let (result, _) = two.div(&sqrt_pi, RoundingMode::NearestEven);
    result
        .round_to_precision(prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}

/// Compute `ln(2π)` at the requested precision via
/// `ln(2) + ln(π)`. `ln(π)` uses the atanh-series identity
/// `ln(π) = 2 · atanh((π − 1) / (π + 1))` evaluated at working
/// precision; the convergence factor `((π−1)/(π+1))² ≈ 0.267`
/// gives roughly `log₂(1/0.267) ≈ 1.9` bits per term.
#[cfg(feature = "specials")]
pub(super) fn ln_2pi_via_agm(prec: u32) -> BigFloat {
    cache::memoized(Kind::Ln2Pi, prec, || ln_2pi_compute(prec))
}

#[cfg(feature = "specials")]
fn ln_2pi_compute(prec: u32) -> BigFloat {
    let working = prec.saturating_add(GUARD_BITS);
    let ln_2 = ln_2_via_atanh(working);
    let ln_pi = {
        let pi = pi_via_agm(working);
        let one = BigFloat::try_from_i64_exact(1, working).expect("precision >= 1");
        let (pi_minus_1, _) = pi.sub(&one, RoundingMode::NearestEven);
        let (pi_plus_1, _) = pi.add(&one, RoundingMode::NearestEven);
        let (u, _) = pi_minus_1.div(&pi_plus_1, RoundingMode::NearestEven);
        let (u_squared, _) = u.mul(&u, RoundingMode::NearestEven);

        let mut sum = u.clone();
        let mut term = u;
        let max_iter = working.saturating_mul(4).max(2048);
        for k in 1u32..=max_iter {
            let (next_term, _) = term.mul(&u_squared, RoundingMode::NearestEven);
            term = next_term;
            let divisor = BigFloat::try_from_i64_exact(i64::from(2 * k + 1), working)
                .expect("precision >= 1");
            let (term_over_d, _) = term.div(&divisor, RoundingMode::NearestEven);
            let (next_sum, _) = sum.add(&term_over_d, RoundingMode::NearestEven);
            sum = next_sum;
            if let Class::Normal { exponent, .. } = &term.class {
                if *exponent < -i64::from(working) - 4 {
                    break;
                }
            } else {
                break;
            }
        }
        let two = BigFloat::try_from_i64_exact(2, working).expect("precision >= 1");
        let (ln_pi, _) = two.mul(&sum, RoundingMode::NearestEven);
        ln_pi
    };
    let (result, _) = ln_2.add(&ln_pi, RoundingMode::NearestEven);
    result
        .round_to_precision(prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}

/// Compute the Euler–Mascheroni constant `γ` at the requested
/// precision via the Brent–McMillan algorithm B1.
///
/// Identity (Brent–McMillan 1980, algorithm B1): with
/// `t_k = (nᵏ / k!)²`, `H_k = Σ_{j=1}^{k} 1/j` (and `H_0 = 0`),
///
/// ```text
/// I(n) = Σ_{k≥0} t_k
/// S(n) = Σ_{k≥0} t_k · H_k
/// γ    = S(n)/I(n) − ln(n) + O(π·e^{−4n})
/// ```
///
/// The truncation error decays like `e^{−4n}`, so `n` is chosen so
/// `4n ≥ (working + slack)·ln 2`. Using the integer ratio
/// `7/40 > (ln 2)/4` keeps the choice no_std-clean and conservative.
/// The inner sums peak near `k = n` and then decay super-geometrically;
/// the loop runs until a term is negligible relative to the running
/// `I(n)` (and past the peak), with a hard cap so the loop stays
/// bounded for pathological inputs. This is a derivation from the
/// published identity, not a port of any implementation.
pub(super) fn euler_gamma_via_bm(prec: u32) -> BigFloat {
    cache::memoized(Kind::EulerGamma, prec, || euler_gamma_compute(prec))
}

fn euler_gamma_compute(prec: u32) -> BigFloat {
    let working = prec.saturating_add(GUARD_BITS);

    // n with 4n ≥ working·ln2: 7/40 = 0.175 > (ln2)/4 ≈ 0.17329.
    let n: i64 = i64::from((working.saturating_mul(7) / 40).saturating_add(3));

    let n_big = BigFloat::try_from_i64_exact(n, working).expect("precision >= 1");
    let (n_sq, _) = n_big.mul(&n_big, RoundingMode::NearestEven);

    // k = 0: t_0 = 1, H_0 = 0 ⇒ I starts at 1, S starts at 0.
    let mut term = BigFloat::try_from_i64_exact(1, working).expect("precision >= 1");
    let mut harmonic = BigFloat::try_new_zero(Sign::Positive, working).expect("precision >= 1");
    let mut i_sum = BigFloat::try_from_i64_exact(1, working).expect("precision >= 1");
    let mut s_sum = BigFloat::try_new_zero(Sign::Positive, working).expect("precision >= 1");

    // Negligible once a term sits `working + 8` bits below the running
    // I(n); the `k > n` guard keeps the rising prefix from stopping
    // early. The cap bounds the loop irrespective of convergence.
    let max_iter: i64 = n.saturating_mul(8).saturating_add(64);
    for k in 1..=max_iter {
        let k_big = BigFloat::try_from_i64_exact(k, working).expect("precision >= 1");
        let (k_sq, _) = k_big.mul(&k_big, RoundingMode::NearestEven);

        // t_k = t_{k-1} · n² / k².
        let (t_n2, _) = term.mul(&n_sq, RoundingMode::NearestEven);
        let (t_next, _) = t_n2.div(&k_sq, RoundingMode::NearestEven);
        term = t_next;

        // H_k = H_{k-1} + 1/k.
        let one = BigFloat::try_from_i64_exact(1, working).expect("precision >= 1");
        let (inv_k, _) = one.div(&k_big, RoundingMode::NearestEven);
        let (h_next, _) = harmonic.add(&inv_k, RoundingMode::NearestEven);
        harmonic = h_next;

        let (i_next, _) = i_sum.add(&term, RoundingMode::NearestEven);
        i_sum = i_next;
        let (t_h, _) = term.mul(&harmonic, RoundingMode::NearestEven);
        let (s_next, _) = s_sum.add(&t_h, RoundingMode::NearestEven);
        s_sum = s_next;

        if k > n {
            let negligible = match (&term.class, &i_sum.class) {
                (Class::Zero { .. }, _) => true,
                (
                    Class::Normal {
                        exponent: t_exp, ..
                    },
                    Class::Normal {
                        exponent: i_exp, ..
                    },
                ) => *t_exp < *i_exp - i64::from(working) - 8,
                _ => false,
            };
            if negligible {
                break;
            }
        }
    }

    let (ratio, _) = s_sum.div(&i_sum, RoundingMode::NearestEven);
    let (ln_n, _) = n_big.ln(RoundingMode::NearestEven);
    let (gamma, _) = ratio.sub(&ln_n, RoundingMode::NearestEven);
    gamma
        .round_to_precision(prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pins the regenerated `LN2_LIMBS_1024` table (slice 7b2) to two
    /// independent primary derivations of the correctly-rounded
    /// 1024-bit `ln(2)`: the authoritative 1100-digit decimal
    /// `LN2_REFERENCE` parsed by the bit-exact decimal parser, and
    /// the in-repo AGM atanh series computed at 2048 bits and rounded
    /// down to 1024. All three must agree bit-for-bit. Any future
    /// edit to the table, the parser, or the AGM kernel that disturbs
    /// the 1024-bit value fails here.
    #[test]
    fn ln_2_table_is_correctly_rounded_at_p1024() {
        use core::cmp::Ordering;
        let table = super::super::ln_2_via_table(1024);
        let reference = BigFloat::parse_str(LN2_REFERENCE, 1024, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let agm_hi = ln_2_via_atanh(2048)
            .round_to_precision(1024, RoundingMode::NearestEven)
            .unwrap()
            .0;
        assert_eq!(
            table.partial_cmp(&reference).0,
            Some(Ordering::Equal),
            "regenerated table must equal the authoritative decimal at p=1024:\n  table={table}\n  ref  ={reference}"
        );
        assert_eq!(
            table.partial_cmp(&agm_hi).0,
            Some(Ordering::Equal),
            "regenerated table must equal the independent AGM derivation at p=1024:\n  table={table}\n  agm  ={agm_hi}"
        );
    }

    fn close_within(a: &BigFloat, b: &BigFloat, bits: u32) -> bool {
        use core::cmp::Ordering;
        let (diff, _) = a.sub(b, RoundingMode::NearestEven);
        let abs_diff = diff.abs();
        if abs_diff.is_zero() {
            return true;
        }
        let p = a.precision().max(b.precision());
        let one = BigFloat::try_from_i64_exact(1, p).unwrap();
        let two = BigFloat::try_from_i64_exact(2, p).unwrap();
        let abs_b = b.abs();
        let mut bound = if abs_b.is_zero() { one } else { abs_b };
        for _ in 0..bits {
            bound = bound.div(&two, RoundingMode::NearestEven).0;
        }
        matches!(
            abs_diff.partial_cmp(&bound).0,
            Some(Ordering::Less | Ordering::Equal)
        )
    }

    /// Authoritative high-precision ln(2) decimal value (truncated
    /// at 1100 decimal digits) for cross-checking the AGM-based
    /// implementation. Source: mathematical-constant reference
    /// value (treated as a fact); cross-checked across multiple
    /// authoritative arbitrary-precision libraries.
    const LN2_REFERENCE: &str =
        "0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633\
         269964186875420014810205706857336855202357581305570326707516350759619307275708283714351903\
         070386238916734711233501153644979552391204751726815749320651555247341395258829504530070953\
         263666426541042391578149520437404303855008019441706416715186447128399681717845469570262716\
         310645461502572074024409266063812120472250259404571525907706178843853479326160793321461063\
         229630667004174887399405193455046752905726793181080078681655055956436321214105735100354087\
         787919495015321866823999361974478928553770776263995860723087366603058241647007054710076761\
         566064155715620770193768712050913923489876727076708502185918872946540390210712092834054812\
         920028324859014194977100138908915636011340902091881879477441960541933027880181437834884330\
         420106892098524138420296010060672085522859500336916517706660608854020525717637145872236018\
         9396977000179565457814057";

    /// Authoritative high-precision Euler–Mascheroni γ decimal
    /// (1102 fractional digits) for pinning `EULER_GAMMA_LIMBS_1024`
    /// and cross-checking the Brent–McMillan implementation. Source:
    /// OEIS A001620, generated by its published recipe
    /// `sympy.S.EulerGamma.n(d)` (mpmath backend); treated as a
    /// mathematical fact, the `LN2_REFERENCE` pattern.
    const EULER_GAMMA_REFERENCE: &str =
        "0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470\
         632917467495146314472498070824809605040144865428362241739976449235362535003337429373377376\
         739427925952582470949160087352039481656708532331517766115286211995015079847937450857057400\
         299213547861466940296043254215190587755352673313992540129674205137541395491116851028079842\
         348775872050384310939973613725530608893312676001724795378367592713515772261027349291394079\
         843010341777177808815495706610750101619166334015227893586796549725203621287922655595366962\
         817638879272680132431010476505963703947394957638906572967929601009015125195950922243501409\
         349871228247949747195646976318506676129063811051824197444867836380861749455169892792301877\
         391072945781554316005002182844096053772434203285478367015177394398700302370339518328690001\
         558193988042707411542227819716523011073565833967348717650491941812300040654693142999297779\
         569303100503086303418569803231083691640025892970890985486825777364288253954925873629596133\
         298574739302373438847070370284412920166417850248733379080562754998434590761643167103146710\
         722370021810745044418664";

    /// Pins `EULER_GAMMA_LIMBS_1024` (slice 6m0) to two independent
    /// primary derivations of the correctly-rounded 1024-bit γ: the
    /// authoritative OEIS A001620 decimal parsed by the bit-exact
    /// parser, and the in-repo Brent–McMillan computation at 2048
    /// bits rounded down to 1024. All three must agree bit-for-bit.
    /// Because Brent–McMillan is an independent code path, a
    /// transcription error in the reference cannot pass here.
    #[test]
    fn euler_gamma_table_is_correctly_rounded_at_p1024() {
        use core::cmp::Ordering;
        let table = super::super::euler_gamma_via_table(1024);
        let reference = BigFloat::parse_str(EULER_GAMMA_REFERENCE, 1024, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let bm_hi = euler_gamma_via_bm(2048)
            .round_to_precision(1024, RoundingMode::NearestEven)
            .unwrap()
            .0;
        assert_eq!(
            table.partial_cmp(&reference).0,
            Some(Ordering::Equal),
            "γ table must equal the OEIS A001620 reference at p=1024:\n  table={table}\n  ref  ={reference}"
        );
        assert_eq!(
            table.partial_cmp(&bm_hi).0,
            Some(Ordering::Equal),
            "γ table must equal the independent Brent–McMillan derivation at p=1024:\n  table={table}\n  bm   ={bm_hi}"
        );
    }

    /// Authoritative high-precision `2/√π` decimal (1115 fractional
    /// digits) for pinning `TWO_OVER_SQRT_PI_LIMBS_1024` and
    /// cross-checking `two_over_sqrt_pi_via_agm`. Generated by
    /// `mpmath` (Python) at `mp.dps = 1130` via
    /// `mpmath.mpf(2) / mpmath.sqrt(mpmath.pi)`; treated as a
    /// mathematical fact, the `LN2_REFERENCE` pattern. The exact
    /// reproduction recipe is recorded in the slice-6m-audit ADR.
    #[cfg(feature = "specials")]
    const TWO_OVER_SQRT_PI_REFERENCE: &str =
        "1.1283791670955125738961589031215451716881012586579977136881714434212849368829868289734873\
         204042147268860566958127234147033798629896523257327309790400355379865856752741191968795207\
         049287004359451424231604915456404411090170543464332444169266162227990255269089720461364753\
         818374903174932317026021327967155439987546683207155977523334881524660787604327012032872433\
         924701009166250638937589133125766516310432488690977314063797548617635563658967789502170018\
         369170684432635651786705036660240492451244474498945400677948625285993188527008566089807266\
         316078753919712163186756584411147658475764631584662115239295549365061803431236161190444592\
         352649307180801706885897250057894784328362385486195484511397575915580997496382738744793841\
         457212668495359398972191775260872674529117575030861618688394769665769827583507237913270184\
         826978506176607899308116821145082965496469503494840187939766833554297711783356674789971831\
         633002753719773724087928258145738579276147546534622368573576042310497324379438177936150629\
         902402949105435182596305442441264686935148305205123699645557899739050603338993937132772461\
         8440884972286626212220347940863310887";

    /// Authoritative high-precision `ln(2π)` decimal (1115 fractional
    /// digits) for pinning `LN_2PI_LIMBS_1024` and cross-checking
    /// `ln_2pi_via_agm`. Generated by `mpmath` (Python) at
    /// `mp.dps = 1130` via `mpmath.log(2 * mpmath.pi)`; treated as a
    /// mathematical fact, the `LN2_REFERENCE` pattern. The exact
    /// reproduction recipe is recorded in the slice-6m-audit ADR.
    #[cfg(feature = "specials")]
    const LN_2PI_REFERENCE: &str =
        "1.8378770664093454835606594728112352797227949472755668256343030809655313918545207953894865\
         972719083952440112932492686748927337257636815871443117518304453627872071214850947173380927\
         918119827616112603264697461892547492510365033899089548201917187027839632231962611480106953\
         907721299179844624279113855486999422005670391966389850627885412925913729488231249524260974\
         736305689987586887646607970258953093145638634759757061713788462725643079461672052950585309\
         829800787111999992074126943705144047152430700687247592054316975009722719076849626583582485\
         399922753679280302789575459100202066417683936712388159514332525411750507649724518605059042\
         160990362403936104519600917610771497670658882278136156555534754445076266765187901482804052\
         386787426337408944137118915686982655208159082601536796094035051774961877174911446465066877\
         848938559655749937054225161751623317487505801769689661835077881525919088198969357960783242\
         618144657028735729075124759420708690852634755752923440722283452753593767913238054014882609\
         582282799976925761217812723574091548090088859200013721780671774949241617759590438569372865\
         7385345545108582901661561895442972855";

    /// Pins `TWO_OVER_SQRT_PI_LIMBS_1024` (slice 4b, audited slice
    /// 6m-audit) to two independent primary derivations of the
    /// correctly-rounded 1024-bit `2/√π`: the authoritative `mpmath`
    /// decimal parsed by the bit-exact parser, and the in-repo
    /// `2 / sqrt(π)` AGM derivation computed at 2048 bits rounded down
    /// to 1024. All three must agree bit-for-bit. Because the AGM path
    /// is independent code, a transcription error in the reference
    /// cannot pass here. Closes the LN2-defect backlog item: the
    /// `prec <= 1024` fast-path table is now guarded against the
    /// "correct through ~450 bits then diverges" failure mode.
    #[cfg(feature = "specials")]
    #[test]
    fn two_over_sqrt_pi_table_is_correctly_rounded_at_p1024() {
        use core::cmp::Ordering;
        let table = super::super::two_over_sqrt_pi_via_table(1024);
        let reference =
            BigFloat::parse_str(TWO_OVER_SQRT_PI_REFERENCE, 1024, RoundingMode::NearestEven)
                .unwrap()
                .0;
        let agm_hi = two_over_sqrt_pi_via_agm(2048)
            .round_to_precision(1024, RoundingMode::NearestEven)
            .unwrap()
            .0;
        assert_eq!(
            table.partial_cmp(&reference).0,
            Some(Ordering::Equal),
            "2/√π table must equal the authoritative decimal at p=1024:\n  table={table}\n  ref  ={reference}"
        );
        assert_eq!(
            table.partial_cmp(&agm_hi).0,
            Some(Ordering::Equal),
            "2/√π table must equal the independent AGM derivation at p=1024:\n  table={table}\n  agm  ={agm_hi}"
        );
    }

    /// Pins `LN_2PI_LIMBS_1024` (slice 4b, audited slice 6m-audit) to
    /// two independent primary derivations of the correctly-rounded
    /// 1024-bit `ln(2π)`: the authoritative `mpmath` decimal parsed by
    /// the bit-exact parser, and the in-repo `ln(2) + ln(π)` atanh-AGM
    /// derivation computed at 2048 bits rounded down to 1024. All three
    /// must agree bit-for-bit. Because the AGM path is independent
    /// code, a transcription error in the reference cannot pass here.
    /// Closes the LN2-defect backlog item: the `prec <= 1024`
    /// fast-path table is now guarded against the "correct through
    /// ~450 bits then diverges" failure mode.
    #[cfg(feature = "specials")]
    #[test]
    fn ln_2pi_table_is_correctly_rounded_at_p1024() {
        use core::cmp::Ordering;
        let table = super::super::ln_2pi_via_table(1024);
        let reference = BigFloat::parse_str(LN_2PI_REFERENCE, 1024, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let agm_hi = ln_2pi_via_agm(2048)
            .round_to_precision(1024, RoundingMode::NearestEven)
            .unwrap()
            .0;
        assert_eq!(
            table.partial_cmp(&reference).0,
            Some(Ordering::Equal),
            "ln(2π) table must equal the authoritative decimal at p=1024:\n  table={table}\n  ref  ={reference}"
        );
        assert_eq!(
            table.partial_cmp(&agm_hi).0,
            Some(Ordering::Equal),
            "ln(2π) table must equal the independent AGM derivation at p=1024:\n  table={table}\n  agm  ={agm_hi}"
        );
    }

    /// Computing γ at 2048 bits and rounding to 512 must match a
    /// direct 512-bit computation: isolates any bug to the
    /// round-to-precision step or earlier (the
    /// `ln_2_via_atanh_self_consistent` pattern).
    #[test]
    fn euler_gamma_via_bm_self_consistent() {
        use core::cmp::Ordering;
        let high = euler_gamma_via_bm(2048)
            .round_to_precision(512, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let direct = euler_gamma_via_bm(512);
        assert_eq!(
            high.partial_cmp(&direct).0,
            Some(Ordering::Equal),
            "γ self-consistency at p=512:\n  rounded from p=2048: {high}\n  direct at p=512:     {direct}"
        );
    }

    #[test]
    fn ln_2_via_atanh_matches_reference_at_p1024() {
        let computed = ln_2_via_atanh(1024);
        let reference = BigFloat::parse_str(LN2_REFERENCE, 1024, RoundingMode::NearestEven)
            .unwrap()
            .0;
        // Tolerate up to 2 ULPs at p=1024.
        let (diff, _) = computed.sub(&reference, RoundingMode::NearestEven);
        let abs = diff.abs();
        let bound = {
            let p = 1024u32;
            let two = BigFloat::try_from_i64_exact(2, p).unwrap();
            let mut b = reference.abs();
            for _ in 0..1022 {
                b = b.div(&two, RoundingMode::NearestEven).0;
            }
            b
        };
        use core::cmp::Ordering;
        assert!(
            matches!(
                abs.partial_cmp(&bound).0,
                Some(Ordering::Less | Ordering::Equal)
            ),
            "atanh ln(2) at p=1024 vs reference:\n  abs  = {abs}\n  bound= {bound}"
        );
    }

    #[test]
    fn ln_2_via_atanh_self_consistent() {
        // Compute at p=2048, round to p=512, compare to direct
        // p=512 computation. Any mismatch isolates the bug to the
        // round-to-precision step (or earlier).
        let high = ln_2_via_atanh(2048);
        let rounded = high
            .round_to_precision(512, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let direct = ln_2_via_atanh(512);
        use core::cmp::Ordering;
        assert_eq!(
            rounded.partial_cmp(&direct).0,
            Some(Ordering::Equal),
            "self-consistency at p=512:\n  rounded from p=2048: {rounded}\n  direct at p=512:     {direct}"
        );
    }

    #[test]
    fn ln_2_via_atanh_at_p128() {
        // Quick mid-range probe: at p=128 the result should agree
        // with the parsed 50-decimal-digit reference value to within
        // 1 ULP.
        let computed = ln_2_via_atanh(128);
        let expected = BigFloat::parse_str(
            "0.69314718055994530941723212145817656807550013436025",
            128,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        let (diff, _) = computed.sub(&expected, RoundingMode::NearestEven);
        let abs = diff.abs();
        let bound = {
            let p = 128u32;
            let two = BigFloat::try_from_i64_exact(2, p).unwrap();
            let mut b = expected.abs();
            for _ in 0..120 {
                b = b.div(&two, RoundingMode::NearestEven).0;
            }
            b
        };
        use core::cmp::Ordering;
        assert!(
            matches!(
                abs.partial_cmp(&bound).0,
                Some(Ordering::Less | Ordering::Equal)
            ),
            "ln(2) at p=128: computed={computed}, expected={expected}, abs={abs}, bound={bound}"
        );
    }

    #[test]
    fn ln_2_table_matches_agm_at_cap_precision() {
        // Seamless-boundary guard: at `prec = LN2_TABLE_PRECISION_CAP`
        // (1024 post slice 7b2) the rounded table value and the AGM
        // atanh series must be bit-identical, so `ln_2_at` has no
        // discontinuity at the table/AGM dispatch boundary. Pre-7b2
        // this held only because the cap sat below the table's faulty
        // ~450-bit range; post-7b2 the table is correctly rounded to
        // the full 1024 bits, so the equality holds at the boundary
        // itself. Any future change to the table, the cap, or the AGM
        // kernel that reintroduces a boundary discontinuity fails here.
        let cap = super::super::LN2_TABLE_PRECISION_CAP;
        let table = super::super::ln_2_via_table(cap);
        let agm = ln_2_via_atanh(cap);
        use core::cmp::Ordering;
        assert_eq!(
            table.partial_cmp(&agm).0,
            Some(Ordering::Equal),
            "LN2 table vs AGM at p={cap}:\n  table={table}\n  agm  ={agm}"
        );
    }

    #[test]
    fn ln_2_via_atanh_at_low_precision() {
        // ln(2) ≈ 0.6931471805599453 — sanity check that the
        // series computes ln(2) at all at low precision before
        // we chase down the 1024-bit divergence.
        let computed = ln_2_via_atanh(53);
        let expected = BigFloat::parse_str("0.6931471805599453", 53, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let (diff, _) = computed.sub(&expected, RoundingMode::NearestEven);
        let abs = diff.abs();
        let one_ulp = BigFloat::parse_str("0.0000000000000001", 53, RoundingMode::NearestEven)
            .unwrap()
            .0;
        use core::cmp::Ordering;
        assert!(
            matches!(
                abs.partial_cmp(&one_ulp).0,
                Some(Ordering::Less | Ordering::Equal)
            ),
            "ln(2) via atanh at p=53: got {computed}, expected {expected}, diff {abs}"
        );
    }

    #[cfg(feature = "trig")]
    #[test]
    fn pi_agm_matches_hardcoded_at_1024_bits() {
        let computed = pi_via_agm(1024);
        let hardcoded = super::super::pi_via_table(1024);
        assert!(
            close_within(&computed, &hardcoded, 1000),
            "Brent-Salamin π at 1024 bits diverges from hardcoded table"
        );
    }

    #[test]
    fn ln_10_via_atanh_at_p128() {
        // ln(10) ≈ 2.302585092994045684017991454684364207601...
        let computed = ln_10_via_atanh(128);
        let expected = BigFloat::parse_str(
            "2.302585092994045684017991454684364207601",
            128,
            RoundingMode::NearestEven,
        )
        .unwrap()
        .0;
        let (diff, _) = computed.sub(&expected, RoundingMode::NearestEven);
        let abs = diff.abs();
        let bound = {
            let p = 128u32;
            let two = BigFloat::try_from_i64_exact(2, p).unwrap();
            let mut b = expected.abs();
            for _ in 0..120 {
                b = b.div(&two, RoundingMode::NearestEven).0;
            }
            b
        };
        use core::cmp::Ordering;
        assert!(
            matches!(
                abs.partial_cmp(&bound).0,
                Some(Ordering::Less | Ordering::Equal)
            ),
            "ln(10) at p=128: abs={abs}, bound={bound}"
        );
    }

    /// A second request at the same precision must be a cache hit:
    /// the entry count does not grow and the returned value is
    /// bit-identical to the first. `π` has no recursive constant, so
    /// one call yields exactly one entry.
    #[cfg(feature = "std")]
    #[test]
    fn memoization_hit_is_bit_identical_and_adds_no_entry() {
        cache::reset();
        assert_eq!(cache::entry_count(), 0);

        let first = pi_via_agm(256);
        assert_eq!(cache::entry_count(), 1, "miss should insert one entry");

        let second = pi_via_agm(256);
        assert_eq!(
            cache::entry_count(),
            1,
            "hit must not insert a duplicate entry"
        );

        use core::cmp::Ordering;
        assert_eq!(
            first.partial_cmp(&second).0,
            Some(Ordering::Equal),
            "cached π must be bit-identical to the freshly computed one"
        );
    }

    /// Recursive constants populate every level: `ln(2π)` at p=128
    /// computes `π` and `ln(2)` at the 192-bit working precision plus
    /// the p=128 result itself — three distinct `(kind, prec)`
    /// entries. A repeat call adds none and the lookup borrow
    /// releasing before recompute is what keeps this from a `RefCell`
    /// double borrow.
    // `ln_2pi_via_agm` is `specials`-gated (pf-gwum), so this recursive
    // memoization check needs both `std` (the cache thread-local) and
    // `specials`.
    #[cfg(all(feature = "std", feature = "specials"))]
    #[test]
    fn memoization_caches_recursive_constants() {
        cache::reset();
        let first = ln_2pi_via_agm(128);
        assert_eq!(
            cache::entry_count(),
            3,
            "ln(2π)@128 should cache Ln2Pi@128, Ln2@192, Pi@192"
        );

        let second = ln_2pi_via_agm(128);
        assert_eq!(
            cache::entry_count(),
            3,
            "repeat call must hit every level, adding no entries"
        );
        use core::cmp::Ordering;
        assert_eq!(
            first.partial_cmp(&second).0,
            Some(Ordering::Equal),
            "cached ln(2π) must be bit-identical to the first result"
        );
    }
}
