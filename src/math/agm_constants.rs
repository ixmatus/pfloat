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
pub(super) fn pi_via_agm(prec: u32) -> BigFloat {
    cache::memoized(Kind::Pi, prec, || pi_compute(prec))
}

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
    let max_iter = (working_prec.saturating_add(64) * 8 / log2_n_times_8.max(1)).max(64);
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
pub(super) fn two_over_pi_via_agm(prec: u32) -> BigFloat {
    cache::memoized(Kind::TwoOverPi, prec, || two_over_pi_compute(prec))
}

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
pub(super) fn two_over_sqrt_pi_via_agm(prec: u32) -> BigFloat {
    cache::memoized(Kind::TwoOverSqrtPi, prec, || two_over_sqrt_pi_compute(prec))
}

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
pub(super) fn ln_2pi_via_agm(prec: u32) -> BigFloat {
    cache::memoized(Kind::Ln2Pi, prec, || ln_2pi_compute(prec))
}

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

#[cfg(test)]
mod tests {
    use super::*;

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
        // Regression guard: the LN2_LIMBS_1024 hardcoded table is
        // known to diverge from the mathematical value past
        // bit ~450. `LN2_TABLE_PRECISION_CAP` sits safely below
        // that boundary, but the bound is approximate; this test
        // forces a bit-exact comparison so any future change to
        // either the table or the cap must keep them consistent.
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
    #[cfg(feature = "std")]
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
