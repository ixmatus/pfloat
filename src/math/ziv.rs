//! Ziv correct-rounding driver.
//!
//! The interval test, factored out of [`crate::math::pow`] at slice
//! p1.2 so the rest of the elementary-transcendental surface can move
//! off the slice-3a fixed-64-bit-guard convention without each kernel
//! reimplementing the loop. The pow kernel was the first caller
//! (ADR-0022, slice 7c); subsequent slices wire `exp`, `ln`, `tanh`,
//! and `lgamma` through [`ziv_round`] and inherit the same correctness
//! argument.
//!
//! The driver is the realization of DESIGN.md §"Ziv's strategy": the
//! caller supplies a closure that evaluates the function at variable
//! working precision; the driver rounds the value to the target under
//! the caller's mode, then asks whether both ends of the bounded
//! evaluation-error interval round to the same value. If they do, the
//! true value rounds there too and the result is correctly rounded;
//! otherwise a rounding boundary lies inside the uncertainty, the
//! guard doubles, and the loop retries, bounded by
//! [`ZIV_MAX_ITERS`]. The cap is the honest measure-zero caveat MPFR
//! also documents.
//!
//! The interval test is the sound termination criterion. Comparing
//! two adjacent guards' rounded values would false-converge on a
//! hard-to-round input (both insufficient guards agree on the same
//! wrong value); the interval test does not.

use core::cmp::Ordering;

use crate::big::BigFloat;
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;
use crate::status::{auto_raise, Status};

/// First Ziv guard: the initial evaluation uses
/// `target + ZIV_BASE_GUARD` extra bits. `pub(super)` so callers (like
/// `sin_kernel` / `cos_kernel`) can derive their range-cap pre-checks
/// from the same constant the driver starts from (pf-1axr).
pub(super) const ZIV_BASE_GUARD: u32 = 64;

/// Maximum extra guard bits above the target precision. The doubling
/// schedule (64, 128, 256, 512, 1024) reaches this at the last
/// iteration.
const ZIV_GUARD_CAP: u32 = 1024;

/// Maximum guard-doubling iterations. On the measure-zero exact-tie
/// inputs that exhaust this many iterations the result may be 1 ULP
/// off in directed modes — the honest caveat MPFR also documents
/// (DESIGN.md §"Ziv's strategy", lines 287-299).
pub(crate) const ZIV_MAX_ITERS: u32 = 5;

// The per-kernel `error_guard` argument supplied at every `ziv_round`
// call site is the slack, in bits below the working precision, the
// caller charges to `eval`'s accumulated `NearestEven` rounding
// error. The half-width `|y|·2^-(working - error_guard)` is then a
// sound upper bound on `|eval(working) − f(x)|` for that kernel.
// The pre-Phase-1g driver carried a single global `ZIV_ERROR_GUARD =
// 24` as documentation-tier assumption (DESIGN.md "Caveats and open
// questions" §1); Phase 1g moves the bound to per-kernel calibrated
// values in `crate::math::ziv_calibration` and forces every caller
// to opt in by name (pf-yupm, ADR-0039). Active sweep-time
// verification of each kernel's bound against the rigorous Arb
// midpoint is pf-tqzz (slice p1g.3).

/// `|y| · 2^-shift`, formed by decrementing the binary exponent (the
/// exact power-of-two scaling used elsewhere in the crate, e.g.
/// `math::mod::pi_over_2_at`). Non-normal `y` has no boundary
/// uncertainty, so a zero half-width is returned.
fn half_width(y: &BigFloat, shift: i64) -> BigFloat {
    match &y.class {
        Class::Normal {
            exponent, mantissa, ..
        } => BigFloat {
            class: Class::Normal {
                sign: Sign::Positive,
                // Saturating (pf-a77o, ADR-0107): a value near the
                // bottom rim made the raw subtraction overflow (a
                // debug panic, wrapped garbage in release). The clamp
                // at i64::MIN OVERSTATES the half-width — the sound
                // direction: an overstated width only refuses
                // certification.
                exponent: exponent.saturating_sub(shift),
                mantissa: mantissa.clone(),
            },
            precision: y.precision,
        },
        _ => BigFloat::try_new_zero(Sign::Positive, y.precision).expect("precision >= 1"),
    }
}

/// Correctly round `eval`'s value to `target` precision under `mode`
/// by the Ziv interval test (DESIGN.md §"Ziv's strategy").
///
/// `eval(working)` returns the kernel's value computed at the working
/// precision with `NearestEven` internal rounding; its error against
/// the true value is bounded by the half-width
/// `|y|·2^-(working − error_guard)`. If both ends of that uncertainty
/// interval round to the same `target`-precision value under `mode`,
/// every point in the interval — including the true value — rounds
/// there too, so that value is correctly rounded. Otherwise a
/// rounding boundary lies within the uncertainty: the guard doubles
/// (capped at [`ZIV_GUARD_CAP`]) and the loop retries, bounded by
/// [`ZIV_MAX_ITERS`]. Comparing two *adjacent* guards would falsely
/// converge on a hard-to-round input (both insufficient guards agree
/// on the wrong value); the interval test does not.
///
/// `error_guard` is the kernel's per-function calibrated bound from
/// [`crate::math::ziv_calibration`] (pf-yupm, ADR-0039). Every call
/// site supplies an explicitly-named constant; the driver carries no
/// implicit default. Active sweep-time verification of each kernel's
/// bound is pf-tqzz (slice p1g.3) via [`ziv_round_capturing`].
pub(crate) fn ziv_round(
    eval: impl Fn(u32) -> BigFloat,
    target: u32,
    mode: RoundingMode,
    error_guard: u32,
) -> (BigFloat, Status) {
    let (cand, status, _converged_working, _eval_intermediate) =
        ziv_round_capturing(eval, target, mode, error_guard);
    (cand, status)
}

/// [`ziv_round`] with an input-derived certification depth, evaluated
/// lazily on cap exhaustion (pf-jl35, ADR-0103).
///
/// The fixed `ZIV_GUARD_CAP` bounds certification at `target + 1024`
/// bits, but a kernel whose input encodes proximity deeper than that
/// (ζ near its pole, trig near a multiple of π/2, the tiny-x series
/// tails) knows — from the input alone — how far the truth can sit
/// from a rounding boundary. The half-width model
/// `|y|·2^-(working − error_guard)` is a per-kernel claim about `eval`
/// valid at *any* working precision, so one further iteration at
/// `target + depth + 64` certifies legitimately where the legacy
/// schedule exhausts. `depth_hint` is called only on exhaustion (zero
/// cost on the common path) and must be input-proportional (the
/// DoS-budget posture); a hint that still leaves the boundary
/// unresolved falls back exactly as before — the measure-zero caveat,
/// one layer deeper.
pub(crate) fn ziv_round_with_depth(
    eval: impl Fn(u32) -> BigFloat,
    target: u32,
    mode: RoundingMode,
    error_guard: u32,
    depth_hint: impl Fn() -> u32,
) -> (BigFloat, Status) {
    let (cand, status, _w, _y) =
        ziv_round_capturing_with_depth(eval, target, mode, error_guard, depth_hint);
    (cand, status)
}

/// Trace returned by [`ziv_round_capturing`]: the rounded candidate
/// at the caller's target precision, the IEEE status, the working
/// precision at which the Ziv interval test converged (or capped),
/// and the `eval(working)` intermediate at the converging iteration.
///
/// The trailing two fields are the quantities pf-tqzz (slice p1g.3)
/// asserts against the rigorous Arb midpoint on every f32 input.
/// Production callers use [`ziv_round`], which destructures with `_`
/// for those fields; the compiler discards the trailing `BigFloat`
/// allocation in the success path. The fallback (cap-exhaustion)
/// path holds one extra `BigFloat` across iterations bounded by
/// [`ZIV_MAX_ITERS`].
pub type ZivTrace = (BigFloat, Status, u32, BigFloat);

// Thread-local capture of the last `ziv_round_capturing` trace.
// Enabled only under the `ziv-instrumented` feature (and unit-test
// builds); off in production so the capture costs nothing. The
// pf-tqzz cross-check harness (slice p1g.3) drains this via
// `take_last_trace` after every kernel call routed through the
// public API (`BigFloat::<fn>_round`, etc.), so the cross-check
// stays generic across the 47 v1.0 kernels without per-kernel
// `_round_capturing` wrapper boilerplate.
#[cfg(all(feature = "std", any(test, feature = "ziv-instrumented")))]
std::thread_local! {
    static LAST_TRACE: core::cell::RefCell<Option<ZivTrace>> =
        const { core::cell::RefCell::new(None) };
}

/// Drain the thread-local trace populated by the most-recent
/// `ziv_round_capturing` (and thus `ziv_round`, which wraps it).
/// Returns `None` when no Ziv-routed call has been made on this
/// thread since the last drain. Production builds without the
/// `ziv-instrumented` feature return `None` unconditionally; the
/// thread-local does not exist there.
///
/// pf-tqzz, slice p1g.3, ADR-0039.
#[cfg(all(feature = "std", any(test, feature = "ziv-instrumented")))]
pub fn take_last_trace() -> Option<ZivTrace> {
    LAST_TRACE.with(|t| t.borrow_mut().take())
}

/// Same as [`ziv_round`] but additionally returns the working
/// precision at which the interval test converged and the
/// `eval(working)` intermediate at that iteration. This is the
/// shape the pf-tqzz cross-check (slice p1g.3) consumes to assert
/// `|eval(working) − rigorous_midpoint| ≤ 2^(error_guard − working)
/// · |rigorous_midpoint|` for every swept f32 input. Production
/// callers use the thin [`ziv_round`] wrapper above.
pub(crate) fn ziv_round_capturing(
    eval: impl Fn(u32) -> BigFloat,
    target: u32,
    mode: RoundingMode,
    error_guard: u32,
) -> ZivTrace {
    ziv_round_capturing_with_depth(eval, target, mode, error_guard, || 0)
}

/// The shared loop behind [`ziv_round`], [`ziv_round_capturing`], and
/// [`ziv_round_with_depth`]: the legacy doubling schedule, then — only
/// on exhaustion, and only when the lazily-computed hint exceeds the
/// legacy cap — one further iteration at `target + depth_hint() + 64`
/// under the identical interval test (ADR-0103).
pub(crate) fn ziv_round_capturing_with_depth(
    eval: impl Fn(u32) -> BigFloat,
    target: u32,
    mode: RoundingMode,
    error_guard: u32,
    depth_hint: impl Fn() -> u32,
) -> ZivTrace {
    let mut guard = ZIV_BASE_GUARD;
    let mut iter = 0u32;
    let mut deep_tried = false;
    let mut legacy_fallback: Option<ZivTrace> = None;
    // Cap (and any deep rung) exhausted breaks out with a best-effort
    // attempt on a pathologically hard input.
    let trace = loop {
        let working = target.saturating_add(guard);
        let y = eval(working);
        let (cand, status) = y
            .round_to_precision(target, mode)
            .expect("target precision >= 1");

        let shift = i64::from(working) - i64::from(error_guard);
        let d = half_width(&y, shift);
        let lo = y.sub(&d, RoundingMode::NearestEven).0;
        let hi = y.add(&d, RoundingMode::NearestEven).0;
        let lo_r = lo.round_to_precision(target, mode).expect("target >= 1").0;
        let hi_r = hi.round_to_precision(target, mode).expect("target >= 1").0;
        if matches!(lo_r.partial_cmp(&hi_r).0, Some(Ordering::Equal)) {
            // The whole uncertainty interval rounds to one value:
            // correct rounding is settled.
            auto_raise(status);
            let trace = (cand, status, working, y);
            #[cfg(all(feature = "std", any(test, feature = "ziv-instrumented")))]
            LAST_TRACE.with(|t| *t.borrow_mut() = Some(trace.clone()));
            return trace;
        }

        let attempt = (cand, status, working, y);
        iter += 1;
        if iter < ZIV_MAX_ITERS {
            legacy_fallback = Some(attempt);
            guard = guard.saturating_mul(2).min(ZIV_GUARD_CAP);
        } else if !deep_tried {
            // Legacy schedule exhausted: ask the kernel how deep the
            // input-encoded structure can park the truth next to a
            // boundary, and certify there if that is past the cap.
            deep_tried = true;
            // The ceiling bounds exponent-encoded depths (an i64 an
            // adversary writes for free): past it the deep rung is
            // refused and the fall-through keeps the documented
            // 1-ulp INEXACT caveat (pf-a77o, ADR-0107). Without it,
            // a saturated hint would request a u32::MAX-bit
            // evaluation — a de facto hang from a 16-byte input.
            // Input-precision-proportional hints below the ceiling
            // (16 Mbit ≈ a 2 MB operand) still certify.
            const ZIV_DEEP_GUARD_CEILING: u32 = 1 << 24;
            let deep_guard = depth_hint().saturating_add(64).min(ZIV_DEEP_GUARD_CEILING);
            if deep_guard <= ZIV_GUARD_CAP {
                break attempt;
            }
            legacy_fallback = Some(attempt);
            guard = deep_guard;
        } else {
            // The deep rung evaluated but did not certify. Its value
            // is the more accurate fallback — unless the deeper
            // working pushed the eval outside the kernel's internal
            // envelope and it degenerated to NaN, in which case the
            // legacy attempt (computed inside the tested envelope)
            // stands.
            if matches!(attempt.0.class, Class::Nan { .. }) {
                break legacy_fallback.expect("the legacy schedule ran to exhaustion");
            }
            break attempt;
        }
    };
    auto_raise(trace.1);
    #[cfg(all(feature = "std", any(test, feature = "ziv-instrumented")))]
    LAST_TRACE.with(|t| *t.borrow_mut() = Some(trace.clone()));
    trace
}

/// Binary exponent of a finite value, or `i64::MIN` for zero / special
/// values (so it sorts below any real operand scale).
// Unused under feature combos that enable a transcendental but no
// special-function consumer (e.g. `big,agm`), matching the calibration
// constants' tolerated-dead-code posture.
#[allow(dead_code)]
pub(crate) fn value_exponent(v: &BigFloat) -> i64 {
    match &v.class {
        Class::Normal { exponent, .. } => *exponent,
        _ => i64::MIN,
    }
}

/// Re-evaluate a cancellation-prone composition at a working precision
/// raised by the realised cancellation depth.
///
/// `eval(w)` returns `(value, operand_scale)`, where `operand_scale` is
/// the binary exponent of the largest magnitude that cancelled to form
/// `value`. When the value falls far below that scale — a near-zero
/// result produced by catastrophic cancellation, such as `lgamma` and
/// `digamma` near their negative-axis roots or `li` near its zero — the
/// value computed at `working_prec` is dominated by accumulated rounding
/// error, so the Ziv interval test's *relative* half-width
/// (`|y|·2^-(working-guard)`) understates the true *absolute* error and
/// would certify a wrong result. We raise the precision by the lost-bit
/// count and re-evaluate, iterating because a probe taken below the
/// cancellation depth under-reports it (the value collapses toward
/// zero). The returned value carries `working_prec` accurate bits
/// relative to its own magnitude, which is exactly the premise the Ziv
/// half-width assumes. Review 2026-05-29 (root cause 2).
///
/// The cancellation depth is bounded by the input's proximity to the
/// zero, itself bounded by the input precision; the iteration cap is a
/// backstop and the Ziv driver remains the outer correctness gate.
#[allow(dead_code)] // unused under transcendental-without-specials combos
pub(crate) fn cancellation_boosted(
    working_prec: u32,
    eval: impl Fn(u32) -> (BigFloat, i64),
) -> BigFloat {
    let mut w = working_prec;
    let mut last = None;
    for _ in 0..12 {
        let (value, operand_scale) = eval(w);
        let result_exp = match &value.class {
            Class::Normal { exponent, .. } => *exponent,
            _ => {
                // Collapsed to zero/special at this precision: the
                // cancellation exceeds w. Double and retry.
                last = Some(value);
                w = w.saturating_mul(2);
                continue;
            }
        };
        let cancel = operand_scale.saturating_sub(result_exp).max(0);
        let cancel = u32::try_from(cancel).unwrap_or(u32::MAX);
        let needed = working_prec.saturating_add(cancel).saturating_add(8);
        if w >= needed {
            return value;
        }
        last = Some(value);
        w = needed;
    }
    last.unwrap_or_else(|| eval(w).0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::ziv_calibration::DEFAULT_ERROR_GUARD;

    #[test]
    fn ziv_round_exact_value_is_bit_exact() {
        // An exact value (8) is identical at every working precision,
        // so the first two guards agree immediately and the result is
        // returned bit-exactly — the MPFR integer-fast-path parity
        // the pow integer branch relies on.
        for &mode in &[
            RoundingMode::NearestEven,
            RoundingMode::TowardPositive,
            RoundingMode::TowardNegative,
            RoundingMode::TowardZero,
            RoundingMode::NearestAway,
        ] {
            let (r, _) = ziv_round(
                |w| BigFloat::try_from_i64_exact(8, w).expect("w >= 1"),
                53,
                mode,
                DEFAULT_ERROR_GUARD,
            );
            let eight = BigFloat::try_from_i64_exact(8, 53).unwrap();
            assert_eq!(
                r.partial_cmp(&eight).0,
                Some(Ordering::Equal),
                "ziv_round exact 8 under {mode:?} = {r}"
            );
            assert_eq!(r.precision, 53);
        }
    }

    #[test]
    fn ziv_round_converges_to_correctly_rounded() {
        // A transcendental-shaped constant: as the guard grows the
        // target-rounding stabilizes, and ziv_round must land on the
        // value obtained by rounding the constant directly to the
        // target precision under the same mode (correct rounding).
        let digits = "1.41421356237309504880168872420969807856967187537694\
                       8073176679737990732478462107038850387534327641572735";
        for &mode in &[
            RoundingMode::NearestEven,
            RoundingMode::TowardPositive,
            RoundingMode::TowardNegative,
            RoundingMode::TowardZero,
            RoundingMode::NearestAway,
        ] {
            let target = 113u32;
            let (r, _) = ziv_round(
                |w| {
                    BigFloat::parse_str(digits, w, RoundingMode::NearestEven)
                        .expect("parse")
                        .0
                },
                target,
                mode,
                DEFAULT_ERROR_GUARD,
            );
            let direct = BigFloat::parse_str(digits, target, mode).expect("parse").0;
            assert_eq!(
                r.partial_cmp(&direct).0,
                Some(Ordering::Equal),
                "ziv_round vs direct round under {mode:?}: ziv={r}, direct={direct}"
            );
        }
    }

    #[test]
    fn ziv_round_is_idempotent_under_recompute() {
        // Stable input → stable output across an independent re-run:
        // the driver does not introduce nondeterminism.
        let eval = |w: u32| {
            BigFloat::parse_str(
                "2.7182818284590452353602874713527",
                w,
                RoundingMode::NearestEven,
            )
            .expect("parse")
            .0
        };
        let (a, _) = ziv_round(eval, 80, RoundingMode::NearestEven, DEFAULT_ERROR_GUARD);
        let (b, _) = ziv_round(eval, 80, RoundingMode::NearestEven, DEFAULT_ERROR_GUARD);
        assert_eq!(a.partial_cmp(&b).0, Some(Ordering::Equal));
    }
}
