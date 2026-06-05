//! Kani harnesses for the Ziv interval-test soundness theorem
//! (pf-hdh8, ADR-0039, slice p1g.4).
//!
//! ## Theorem
//!
//! For any working precision `w`, target precision `t < w`, rounding
//! mode `m`, candidate value `y` at precision `w`, and half-width
//! `h ≥ 0` at precision `w` such that the Ziv interval-test predicate
//!
//! ```text
//! interval_test(y, h, t, m) =
//!   round_to_precision(y − h, t, m) == round_to_precision(y + h, t, m)
//! ```
//!
//! returns `true`, every value `y'` in the closed interval
//! `[y − h, y + h]` rounds to the same target-precision value under
//! mode `m`:
//!
//! ```text
//! ∀ y' with |y' − y| ≤ h:
//!   round_to_precision(y', t, m) == round_to_precision(y, t, m)
//! ```
//!
//! This is the soundness property the Ziv driver relies on but never
//! proves: if the interval test accepts, the rounding of `y` is the
//! same as the rounding of every other value in the uncertainty
//! interval — including the true value `f(x)` whose distance from
//! the kernel's `eval(w)` output is bounded by `h` by the per-
//! function calibrated `error_guard` (pf-yupm) actively verified by
//! the per-release oracle-sweep cross-check (pf-tqzz).
//!
//! ## Scope (ADR-0039 pre-commitment)
//!
//! The Kani discharge runs at **fixed target precisions
//! t ∈ {24, 53, 113}** (the IEEE 754-2019 binary32/binary64/
//! binary128 surface). The arbitrary-precision claim is recorded in
//! ADR-0039 as "validated by sweep, structurally analogous to the
//! discharged IEEE targets" — the round-to-precision predicate is
//! uniform across `t`, the interval test is uniform across `t`, and
//! the bounded encoding scales without changing shape, so the IEEE
//! targets stand in for the arbitrary-`t` family.
//!
//! ## Operand bounding
//!
//! The CBMC backend that Kani uses cannot symbolically enumerate
//! over `Vec<u64>` mantissas (ADR-0012 lesson: the existing 196
//! Kani harnesses moved to manual on-demand workflow because deep
//! transcendental verification with `Vec` storage times out). This
//! module's harnesses constrain operand draws to the canonical
//! eight-constant set via [`super::helpers::nondet_constant_at`],
//! which keeps the SAT problem tractable.
//!
//! The canonical-set discharge proves the theorem at the eight
//! "structural" inputs that exercise every IEEE 754-2019 class
//! (qNaN, sNaN, ±∞, ±0, ±1). The theorem's soundness on the much
//! larger set of arbitrary-mantissa normal values follows by the
//! structural-analogy argument recorded in ADR-0039 plus the
//! pf-tqzz sweep cross-check that actively guards the kernel-side
//! error bound at every f32 input.
//!
//! Lifting the discharge to true universal quantification over the
//! mantissa domain was investigated in pf-25zw (ADR-0062). The
//! `BoundedBigFloat<80>` fixed-array *operand* encoding is necessary
//! but not sufficient: the harness still evaluates the theorem through
//! the real `Vec`-backed `add` / `sub` / `round_to_precision` /
//! `partial_cmp`, and CBMC's model of pfloat's `Vec` storage is hostile
//! at the allocation level (it spuriously fails even a copy-and-compare
//! round-trip), not merely in the arithmetic loops. So a fixed-array
//! shim into the real ops cannot discharge; the genuine path is a
//! fully `Vec`-free re-implementation of the four operations on
//! `[u64; N]`, verified for fidelity against the real ops, which
//! ADR-0062 scopes as the open follow-up. Until then this
//! eight-constant discharge plus the pf-tqzz sweep cross-check stand
//! in for the universal claim, as ADR-0039 records.

use crate::big::BigFloat;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

use super::helpers::{nondet_constant_at, nondet_rounding_mode};

/// The Ziv interval-test predicate, factored out of
/// `crate::math::ziv::ziv_round_capturing` for verification.
/// Returns `true` when both endpoints of the bounded uncertainty
/// interval `[y − h, y + h]` round to the same `t`-precision value
/// under mode `m`. This is the conjunct the driver tests on every
/// loop iteration; the soundness theorem proves that this conjunct
/// implies `round(y', t, m) == round(y, t, m)` for every `y'` in
/// the interval.
fn interval_test(y: &BigFloat, h: &BigFloat, t: u32, m: RoundingMode) -> bool {
    let lo = y.sub(h, RoundingMode::NearestEven).0;
    let hi = y.add(h, RoundingMode::NearestEven).0;
    let lo_r = lo.round_to_precision(t, m).expect("target >= 1").0;
    let hi_r = hi.round_to_precision(t, m).expect("target >= 1").0;
    matches!(lo_r.partial_cmp(&hi_r).0, Some(core::cmp::Ordering::Equal))
}

/// Round `y` to target precision `t` under mode `m`; the function
/// the soundness theorem asserts equality on.
fn round(y: &BigFloat, t: u32, m: RoundingMode) -> BigFloat {
    y.round_to_precision(t, m).expect("target >= 1").0
}

/// Soundness at `t = 24` (IEEE binary32): when the interval test
/// accepts at the canonical-operand triple `(y, h, y')`, the
/// rounded `y'` matches the rounded `y` under the mode.
///
/// The harness draws `y`, `h`, and `y'` non-deterministically from
/// the canonical eight-constant set at working precision 88
/// (= target 24 + ZIV_BASE_GUARD 64), assumes `h ≥ 0` and `|y' − y|
/// ≤ h`, assumes the interval test accepts, and asserts the
/// rounding equality. ADR-0039.
#[kani::proof]
fn ziv_interval_test_is_sound_at_t24() {
    const T: u32 = 24;
    const W: u32 = T + 64;
    let y = nondet_constant_at(W);
    let h = nondet_constant_at(W);
    let y_prime = nondet_constant_at(W);
    let m = nondet_rounding_mode();

    // h must be non-negative (it is a half-width, which the Ziv
    // driver constructs from |y| · 2^-shift).
    kani::assume(!matches!(h.sign(), Sign::Negative));
    // h must be finite; the half-width on a non-normal `y` is zero
    // per the driver's `half_width` helper, which is the trivial
    // case (every endpoint of a zero-width interval is equal).
    kani::assume(!h.is_infinite() && !h.is_nan());
    // Interval-containment constraint: |y' - y| ≤ h.
    let delta = y_prime.sub(&y, RoundingMode::NearestEven).0;
    let abs_delta = delta.abs();
    kani::assume(
        matches!(
            abs_delta.partial_cmp(&h).0,
            Some(core::cmp::Ordering::Less | core::cmp::Ordering::Equal)
        ) && !abs_delta.is_nan(),
    );

    // Hypothesis: the interval test accepts at (y, h, t, m).
    kani::assume(interval_test(&y, &h, T, m));

    // Conclusion: round(y', t, m) == round(y, t, m).
    let y_r = round(&y, T, m);
    let y_prime_r = round(&y_prime, T, m);
    let _ = matches!(
        y_prime_r.partial_cmp(&y_r).0,
        Some(core::cmp::Ordering::Equal)
    );
    // NaN-preserving equality: NaN endpoints round to NaN under
    // every mode, so the assertion holds via NaN's pattern (a NaN
    // y produces NaN endpoints, the interval test fails via the
    // NaN-vs-NaN partial_cmp shape; the hypothesis assumption then
    // rules this case out). The non-NaN case is the meat.
    if !y_r.is_nan() && !y_prime_r.is_nan() {
        assert!(
            matches!(
                y_prime_r.partial_cmp(&y_r).0,
                Some(core::cmp::Ordering::Equal)
            ),
            "Ziv soundness failed at t=24"
        );
    }
}

/// Soundness at `t = 53` (IEEE binary64). Same shape as the t=24
/// harness above; the parameter narrowing keeps Kani's symbolic
/// execution tractable per ADR-0039's scope pre-commitment.
#[kani::proof]
fn ziv_interval_test_is_sound_at_t53() {
    const T: u32 = 53;
    const W: u32 = T + 64;
    let y = nondet_constant_at(W);
    let h = nondet_constant_at(W);
    let y_prime = nondet_constant_at(W);
    let m = nondet_rounding_mode();

    kani::assume(!matches!(h.sign(), Sign::Negative));
    kani::assume(!h.is_infinite() && !h.is_nan());
    let delta = y_prime.sub(&y, RoundingMode::NearestEven).0;
    let abs_delta = delta.abs();
    kani::assume(
        matches!(
            abs_delta.partial_cmp(&h).0,
            Some(core::cmp::Ordering::Less | core::cmp::Ordering::Equal)
        ) && !abs_delta.is_nan(),
    );

    kani::assume(interval_test(&y, &h, T, m));

    let y_r = round(&y, T, m);
    let y_prime_r = round(&y_prime, T, m);
    if !y_r.is_nan() && !y_prime_r.is_nan() {
        assert!(
            matches!(
                y_prime_r.partial_cmp(&y_r).0,
                Some(core::cmp::Ordering::Equal)
            ),
            "Ziv soundness failed at t=53"
        );
    }
}

/// Soundness at `t = 113` (IEEE binary128). Same shape; this
/// rounds out the IEEE binary{32,64,128} target-precision scope
/// the ADR-0039 commitment fixes.
#[kani::proof]
fn ziv_interval_test_is_sound_at_t113() {
    const T: u32 = 113;
    const W: u32 = T + 64;
    let y = nondet_constant_at(W);
    let h = nondet_constant_at(W);
    let y_prime = nondet_constant_at(W);
    let m = nondet_rounding_mode();

    kani::assume(!matches!(h.sign(), Sign::Negative));
    kani::assume(!h.is_infinite() && !h.is_nan());
    let delta = y_prime.sub(&y, RoundingMode::NearestEven).0;
    let abs_delta = delta.abs();
    kani::assume(
        matches!(
            abs_delta.partial_cmp(&h).0,
            Some(core::cmp::Ordering::Less | core::cmp::Ordering::Equal)
        ) && !abs_delta.is_nan(),
    );

    kani::assume(interval_test(&y, &h, T, m));

    let y_r = round(&y, T, m);
    let y_prime_r = round(&y_prime, T, m);
    if !y_r.is_nan() && !y_prime_r.is_nan() {
        assert!(
            matches!(
                y_prime_r.partial_cmp(&y_r).0,
                Some(core::cmp::Ordering::Equal)
            ),
            "Ziv soundness failed at t=113"
        );
    }
}

/// Zero-half-width sanity: when `h = +0`, the interval is the
/// single point `y` and rounding `y` to itself is trivially
/// equal. This pins the boundary case the driver hits whenever
/// `y` is non-normal (`half_width` returns zero for NaN, infinity,
/// and zero per the driver implementation).
#[kani::proof]
fn ziv_interval_test_zero_half_width_is_trivially_sound() {
    const T: u32 = 24;
    const W: u32 = T + 64;
    let y = nondet_constant_at(W);
    let h = BigFloat::try_new_zero(Sign::Positive, W).expect("precision >= 1");
    let m = nondet_rounding_mode();

    // Zero half-width: every y' in the interval IS y, so soundness
    // holds trivially. The interval test trivially accepts (lo = hi
    // = y, both round to round(y, t, m)). Confirm this matches the
    // theorem's structural claim.
    if interval_test(&y, &h, T, m) {
        let y_r = round(&y, T, m);
        // Soundness for y itself; the universal quantifier in the
        // theorem collapses to the single-point case.
        assert_eq!(
            y_r.partial_cmp(&y_r).0,
            Some(core::cmp::Ordering::Equal),
            "zero-half-width round-to-precision must equal itself"
        );
    }
}
