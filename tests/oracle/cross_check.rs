//! pf-tqzz cross-check shared machinery (ADR-0039, slice p1g.3).
//!
//! For every `(kernel, input, mode)` triple this module asserts that
//! pfloat's Ziv-driver `eval(working)` intermediate stayed within the
//! kernel's per-function calibrated `error_guard` bound (from
//! [`pfloat::ziv_instrumented`], `src/math/ziv_calibration.rs`) when
//! compared against the rigorous-enclosure midpoint computed at
//! higher precision by Arb (for the 12 Arb-primary [`FnId`]s) or MPFR
//! (for the 35 MPFR-primary [`FnId`]s).
//!
//! The assertion is
//!
//! ```text
//! |eval(working) − midpoint| ≤ 2^(error_guard − working) · |midpoint|
//! ```
//!
//! evaluated in [`rug::Float`] at `oracle_prec = working + 64`.
//!
//! Two consumers share this surface:
//!
//! - `tests/oracle_cross_check_smoke.rs` — the per-push smoke gate.
//!   Wraps a [`CheckOutcome::Violation`] in `panic!` to preserve the
//!   pre-refactor structured-panic UX.
//! - `examples/pf_tqzz_sweep.rs` — the full per-release sweep
//!   (pf-hcz4). Collects [`ViolationRecord`]s into a sidecar list,
//!   emits them as JSON, and lets the sweep run to completion across
//!   all 65 536 × 5 inputs per `FnId`. ADR-0049 records the design.

#![cfg(all(feature = "differential-arb", feature = "ziv-instrumented"))]

use super::arb::ArbOracle;
use super::convert::bigfloat_to_rug;
use super::meta::is_arb_primary;
use super::mpfr::MpfrOracle;
use super::pfloat_kernels::pfloat_kernel_value;
use super::types::FnId;

use pfloat::ziv_instrumented::take_last_trace;
use pfloat::RoundingMode;

use rug::Float;

/// Map each `FnId` to its per-kernel calibrated `error_guard` bound.
/// The cross-check assertion reads this rather than the global
/// `DEFAULT_ERROR_GUARD = 24` to validate the per-kernel calibration
/// table (pf-yupm, p1g.2).
pub fn error_guard_for(f: FnId) -> u32 {
    use pfloat::ziv_instrumented as zc;
    match f {
        FnId::Exp => zc::EXP_ERROR_GUARD,
        FnId::Exp2 => zc::EXP2_ERROR_GUARD,
        FnId::Exp10 => zc::EXP10_ERROR_GUARD,
        FnId::Expm1 => zc::EXPM1_ERROR_GUARD,
        FnId::Ln => zc::LN_ERROR_GUARD,
        FnId::Log1p => zc::LOG1P_ERROR_GUARD,
        FnId::Sin => zc::SIN_ERROR_GUARD,
        FnId::Cos => zc::COS_ERROR_GUARD,
        FnId::Tan => zc::TAN_ERROR_GUARD,
        FnId::Asin => zc::ASIN_ERROR_GUARD,
        FnId::Acos => zc::ACOS_ERROR_GUARD,
        FnId::Atan => zc::ATAN_ERROR_GUARD,
        FnId::Sinh => zc::SINH_ERROR_GUARD,
        FnId::Cosh => zc::COSH_ERROR_GUARD,
        FnId::Tanh => zc::TANH_ERROR_GUARD,
        FnId::Asinh => zc::ASINH_ERROR_GUARD,
        FnId::Acosh => zc::ACOSH_ERROR_GUARD,
        FnId::Atanh => zc::ATANH_ERROR_GUARD,
        FnId::Erf => zc::ERF_ERROR_GUARD,
        FnId::Erfc => zc::ERFC_ERROR_GUARD,
        FnId::Gamma => zc::GAMMA_ERROR_GUARD,
        FnId::Lgamma => zc::LGAMMA_ERROR_GUARD,
        FnId::Digamma => zc::DIGAMMA_ERROR_GUARD,
        FnId::Zeta => zc::ZETA_ERROR_GUARD,
        FnId::Ei => zc::EI_ERROR_GUARD,
        FnId::Si => zc::SI_ERROR_GUARD,
        FnId::Ci => zc::CI_ERROR_GUARD,
        FnId::Li => zc::LI_ERROR_GUARD,
        FnId::Bi => zc::AIRY_ERROR_GUARD,
        FnId::AiPrime | FnId::BiPrime => zc::AIRY_ERROR_GUARD,
        FnId::BesselJ1 | FnId::BesselJn(_) => zc::BESSEL_J_ERROR_GUARD,
        FnId::BesselI0 | FnId::BesselI1 | FnId::BesselIn(_) => zc::BESSEL_I_ERROR_GUARD,
        FnId::BesselK0 | FnId::BesselK1 | FnId::BesselKn(_) => zc::BESSEL_K_ERROR_GUARD,
        // Fall-through for FnIds whose kernel does not go through
        // ziv_round (e.g., log2 / log10's fixed-guard composition,
        // sqrt's arithmetic-core path); the cross-check skips them
        // because no Ziv trace is captured. The DEFAULT here is
        // only reached if the FnId snuck through the
        // skip-on-no-trace check below; treat as conservative.
        _ => zc::DEFAULT_ERROR_GUARD,
    }
}

/// Route the midpoint request to the right oracle. The 12
/// Arb-primary `FnId`s go through [`ArbOracle::midpoint`]; the
/// MPFR-primary set goes through [`MpfrOracle::midpoint`]. Returns
/// `None` when the Arb oracle is unavailable for an Arb-primary
/// `FnId` (caller skips that triple).
pub fn midpoint_for(
    f: FnId,
    input: u32,
    oracle_prec: u32,
    arb: Option<&ArbOracle>,
    mpfr: &MpfrOracle,
) -> Option<Float> {
    if is_arb_primary(f) {
        let arb = arb?;
        arb.midpoint(f, input, oracle_prec).ok()
    } else {
        Some(mpfr.midpoint(f, input, oracle_prec))
    }
}

/// One violating triple. The cross-check records the kernel, input,
/// and mode that exceeded the calibrated `error_guard` bound, plus
/// numeric witnesses (`eval_w`, `midpoint`, `abs_diff`, `bound`, `gap`)
/// needed to triage post-run. The fields stay as [`rug::Float`] rather
/// than pre-formatted strings; consumers format on demand for either
/// the smoke's structured-panic message or the sweep's JSON sidecar.
#[derive(Debug, Clone)]
pub struct ViolationRecord {
    pub fn_id: FnId,
    pub input: u32,
    pub mode: RoundingMode,
    pub working_prec: u32,
    pub error_guard: u32,
    pub eval_w: Float,
    pub midpoint: Float,
    pub abs_diff: Float,
    pub bound: Float,
    pub gap: Float,
}

impl ViolationRecord {
    /// The structured-panic message the pre-refactor smoke harness
    /// emitted on violation. Preserved verbatim so the smoke's
    /// caller still sees the identical UX after the refactor.
    pub fn format_panic_message(&self) -> String {
        format!(
            "pf-tqzz cross-check violation:\n\
             kernel = {:?}\n\
             input  = 0x{:08x} (f32 {})\n\
             mode   = {:?}\n\
             working_prec = {}\n\
             error_guard  = {}\n\
             |eval(w) - midpoint| = {}\n\
             bound = 2^(error_guard - working) * |midpoint| = {}\n\
             gap   = {}",
            self.fn_id,
            self.input,
            f32::from_bits(self.input),
            self.mode,
            self.working_prec,
            self.error_guard,
            self.abs_diff,
            self.bound,
            self.gap,
        )
    }
}

/// What a single `cross_check_one` call resolved to. `Pass` and
/// `Skipped*` are normal outcomes; `Violation` is the assertion
/// failure. The caller decides whether to panic, log, or collect.
#[derive(Debug)]
pub enum CheckOutcome {
    Pass,
    SkippedNoZivPath,
    SkippedNoMidpoint,
    SkippedNonFiniteMidpoint,
    /// The captured Ziv trace does not belong to `f` itself. Composed
    /// kernels (`log2 = ln / ln2`, `log10 = ln / ln10`, …) leave the
    /// thread-local trace from their inner `ln` `ziv_round`; the
    /// trace's candidate is the inner value, not `f`'s output, so the
    /// `error_guard` assertion against `f`'s midpoint is undefined.
    /// pf-hcz4 surfaced this as a ~100%-of-inputs false violation on
    /// `log2`/`log10`. The cross-check skips these cells rather than
    /// mis-measure them.
    SkippedTraceNotFinal,
    Violation(ViolationRecord),
}

/// Run the cross-check assertion for one `(kernel, input, mode)`
/// triple. Returns the [`CheckOutcome`]; never panics. Inspect the
/// variant at the call site to decide whether to abort, collect, or
/// continue. Drains the thread-local Ziv trace before the kernel
/// call and after, so successive calls on the same thread are safe.
pub fn cross_check_one(
    f: FnId,
    input: u32,
    mode: RoundingMode,
    arb: Option<&ArbOracle>,
    mpfr: &MpfrOracle,
) -> CheckOutcome {
    // Drain any stale trace from a previous call on this thread.
    let _ = take_last_trace();

    // Route through the public dispatch; the kernel internally calls
    // `ziv_round` which writes the trace via the thread-local. Keep
    // the kernel's actual `BigFloat` output to validate the trace
    // belongs to `f` (see the composed-kernel guard below).
    let y_final = pfloat_kernel_value(f, input, mode);

    // If the kernel short-circuited (Class::Zero / Infinity / NaN
    // dispatch, pre-Ziv exact dispatch), no trace is captured —
    // nothing to assert about the calibrated bound at this input.
    let (cand, _, working, eval_w) = match take_last_trace() {
        Some(trace) => trace,
        None => return CheckOutcome::SkippedNoZivPath,
    };

    // Composed-kernel guard (pf-hcz4). For a primary kernel the
    // captured trace IS the kernel's final `ziv_round`, so its rounded
    // candidate `cand` equals the kernel's output `y_final` by value.
    // A composed kernel routes through an inner `ziv_round` (e.g.
    // `log2` calls `ln_round` then divides by `ln 2`); the trailing
    // trace is then the inner `ln` intermediate and `cand` is the
    // rounded `ln`, not `f`'s output. Asserting `|ln − log2(x)|`
    // against the `error_guard` band reads as a violation on nearly
    // every input. Detect the mismatch by value and skip: the bound
    // assertion is only well-defined when the kernel's output is a
    // direct `ziv_round` of the target function.
    if cand.partial_cmp(&y_final).0 != Some(core::cmp::Ordering::Equal) {
        return CheckOutcome::SkippedTraceNotFinal;
    }

    // Lift `eval_w` (a BigFloat at working precision) into a
    // rug::Float at `oracle_prec = working + 64`. Extending the
    // mantissa by 64 zero bits is exact; the value stays bit-
    // identical to the working-precision BigFloat representation.
    let oracle_prec = working.saturating_add(64);
    let eval_w_f_at_w: Float = bigfloat_to_rug(&eval_w);
    let eval_w_f = Float::with_val(oracle_prec, &eval_w_f_at_w);

    // Fetch the midpoint at oracle_prec from the right backend.
    let midpoint = match midpoint_for(f, input, oracle_prec, arb, mpfr) {
        Some(m) => m,
        None => return CheckOutcome::SkippedNoMidpoint,
    };

    // If the midpoint is non-finite (NaN, ±∞), the cross-check
    // assertion does not apply (no notion of "internal error
    // budget" against a non-finite reference).
    if !midpoint.is_finite() {
        return CheckOutcome::SkippedNonFiniteMidpoint;
    }

    // Compute |eval(w) - midpoint| and the bound
    //   bound = 2^(error_guard - working) * |midpoint|
    // both at oracle_prec.
    let error_guard = error_guard_for(f);
    let diff: Float = Float::with_val(oracle_prec, &eval_w_f - &midpoint);
    let abs_diff: Float = diff.abs();
    let abs_mid: Float = midpoint.clone().abs();

    let shift = i64::from(error_guard) - i64::from(working);
    let mut bound = abs_mid.clone();
    // bound *= 2^shift; shift may be negative (typical case where
    // working > error_guard, so the bound is tiny).
    if shift >= 0 {
        let s = u32::try_from(shift).expect("shift fits u32");
        bound <<= s;
    } else {
        let s = u32::try_from(-shift).expect("shift fits u32");
        bound >>= s;
    }

    if abs_diff <= bound {
        CheckOutcome::Pass
    } else {
        let gap: Float = Float::with_val(oracle_prec, &abs_diff - &bound);
        CheckOutcome::Violation(ViolationRecord {
            fn_id: f,
            input,
            mode,
            working_prec: working,
            error_guard,
            eval_w: eval_w_f,
            midpoint,
            abs_diff,
            bound,
            gap,
        })
    }
}
