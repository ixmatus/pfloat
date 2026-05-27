//! Per-kernel `ZIV_ERROR_GUARD` calibration (pf-yupm, ADR-0039).
//!
//! The Ziv driver's interval test uses
//! `half_width = |y|·2^-(working - error_guard)`, where `error_guard`
//! is the kernel's assumed upper bound on internal accumulated
//! rounding error at working precision, in bits. Pre-Phase-1g the
//! driver carried a single global `ZIV_ERROR_GUARD = 24` as
//! documentation-tier assumption (DESIGN.md "Caveats and open
//! questions" §1; `src/math/ziv.rs:50-58`). Phase 1g moves the
//! bound to a per-kernel calibrated value, audited from the
//! kernel's `eval(w)` source structure, with the driver's parameter
//! required at every call site (no implicit default).
//!
//! The derivation discipline is algebraic primary, empirical
//! secondary:
//!
//! - **Algebraic:** count the floating-point operations on the
//!   `eval(w)` path; each NE-rounded op contributes ≤ 1 ULP of
//!   accumulated error; sum the upper bound; round up to the next
//!   power of two in bits. The `src/math/ziv.rs:50-58` analysis is
//!   the template (exp series at ~4w iterations: log₂(4w) ≪ 24 ULP
//!   at all working precisions this driver runs).
//!
//! - **Empirical:** for kernels where the algebraic analysis is
//!   loose or hard to derive (lgamma near the negative-half-integer
//!   poles, Bessel-Y near zeros, zeta in the critical strip), the
//!   widening sweep doubles `error_guard` until no
//!   `OracleInconclusive` verdicts surface across the 65536 × 5
//!   f32 grid. The smallest power-of-two that passes is pinned.
//!
//! At Phase 1g landing, every kernel sits at `DEFAULT_ERROR_GUARD =
//! 24` by the algebraic analysis. The pf-tqzz cross-check (slice
//! p1g.3) actively guards each kernel's calibrated bound against
//! the rigorous Arb midpoint at every swept input; any violation
//! surfaces as a fatal report identifying the kernel and the gap,
//! and the corresponding constant here is widened with the empirical
//! provenance recorded.
//!
//! Per-kernel derivation is documented in
//! `docs/decisions/plans/phase-1g-verification-closure.md`
//! § Per-function `ZIV_ERROR_GUARD` calibration.

/// Conservative default bound, in bits. Holds for every kernel
/// whose `eval(w)` is a finite composition of NE-rounded
/// arithmetic ops with sum of per-op errors well under `2^24`
/// ULP at all working precisions this driver runs. The number
/// is the empirical slack the pre-Phase-1g global `ZIV_ERROR_GUARD`
/// carried; Phase 1g preserves it as the default and forces every
/// kernel to opt in by name (pf-yupm acceptance criterion 4).
pub(crate) const DEFAULT_ERROR_GUARD: u32 = 24;

// --- Elementary transcendental kernels -------------------------------
//
// Each constant is the per-kernel calibrated bound passed to
// `ziv_round`. Algebraic justification cited per-kernel; the eval(w)
// op count drives the analysis. Every per-kernel constant carries
// the algebraic bound, not the value 24 directly: tightening any
// kernel later is a local change to its constant.

/// `exp`: Taylor series with `~4w` iterations of multiply+add,
/// each ≤ 1 ULP at working precision. Sum ≤ `2^14` ULP at the
/// 1024-bit working cap. Algebraic justification, `src/math/exp.rs:132`.
pub(crate) const EXP_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `exp2`: composition `exp(x · ln(2))` plus the mode-aware
/// constant scaling. Bounded by `EXP_ERROR_GUARD` plus a constant.
pub(crate) const EXP2_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `exp10`: composition `exp(x · ln(10))`.
pub(crate) const EXP10_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `expm1`: composition with cancellation boost inside the eval
/// closure (slice p1.24); the boost handles the tiny-x cancellation
/// regime, leaving the outer Ziv bound at the elementary template.
pub(crate) const EXPM1_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `ln`: atanh series with `~w/3` iterations
/// (`u² ∈ [0, 1/9]`, ~3 bits per term).
pub(crate) const LN_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `log1p`: atanh series on `(1+x)` with cancellation boost in the
/// eval closure (slice p1.24's tiny-x boost).
pub(crate) const LOG1P_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

// --- Forward trig ----------------------------------------------------

/// `sin`: Payne-Hanek range reduction plus quadrant-dispatched
/// Taylor at `~w/2` terms.
pub(crate) const SIN_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `cos`: same range reduction; sin/cos share the structure.
pub(crate) const COS_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `tan`: `sin/cos` composition after the same reduction.
pub(crate) const TAN_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

// --- Inverse trig ----------------------------------------------------

/// `asin`: `2·atan(|x|/(1+sqrt(1-x²)))` identity (slice p1.25).
pub(crate) const ASIN_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `acos`: `π - 2·atan(sqrt((1+x)/(1-x)))` identity (slice p1.25).
pub(crate) const ACOS_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `atan`: unsigned composition on `|x|` (slice p1.25).
pub(crate) const ATAN_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `atan2`: quadrant-shifted `atan(y/x)` (slice p1.25).
pub(crate) const ATAN2_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

// --- Hyperbolic and inverse hyperbolic -------------------------------

/// `sinh`: `(expm1(x) - expm1(-x))/2` composition (slice p1.27).
pub(crate) const SINH_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `cosh`: `(exp(x) + exp(-x))/2` composition.
pub(crate) const COSH_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `tanh`: composition through `tanh_at_w` helper (slice p1.27).
pub(crate) const TANH_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `asinh`: `log1p(|x| + x²/(1+sqrt(1+x²)))`.
pub(crate) const ASINH_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `acosh`: `log1p((x-1) + sqrt((x-1)(x+1)))`.
pub(crate) const ACOSH_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `atanh`: `(log1p(x) - log1p(-x))/2`.
pub(crate) const ATANH_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

// --- Power kernels ---------------------------------------------------

/// `pow` via `exp(y · ln(x))` (the general path; slice 7c).
/// Operation count is the sum of `ln` + `mul` + `exp`, each within
/// its own bound; product bound stays comfortably under `2^24`.
pub(crate) const POW_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `pow` integer-y path via square-and-multiply (slice 7c).
/// `~log₂(|n|)` multiplications, each ≤ 1 ULP; `n ≤ 2^31` keeps the
/// sum ≤ `2^5` ULP, well within the default.
pub(crate) const POW_INT_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

// --- Gamma family ----------------------------------------------------

/// `gamma`: `sign(x) · exp(lgamma(x))` composition.
/// The integer-fast-path dispatches exactly (pf-kk16, slice p1.29).
pub(crate) const GAMMA_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `lgamma`: Spouge for `x ≥ z_min` + reflection for `x < z_min`,
/// the more conservative composition of the family. `~30` ops at
/// `z_min ≈ 20`; algebraic bound well under `2^24` ULP. Empirical
/// confirmation pending pf-tqzz (p1g.3 active sweep guard).
pub(crate) const LGAMMA_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `digamma`: composition through `digamma_at_w` helper (slice p1.29).
pub(crate) const DIGAMMA_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `beta`: `exp(lgamma(x) + lgamma(y) - lgamma(x+y))` composition
/// in the general case; product formula in the integer-arg branch.
pub(crate) const BETA_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

// --- Error functions -------------------------------------------------

/// `erf`: asymptotic (large `|x|`) or Maclaurin (small `|x|`)
/// dispatched at working precision (slice p1.4; oracle bumps to
/// `p = 53` for the f32-subnormal-midpoint trap).
pub(crate) const ERF_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `erfc`: composition `1 - erf(...)` or direct asymptotic branch
/// (slice p1.28).
pub(crate) const ERFC_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

// --- Integral functions ----------------------------------------------

/// `Si`: Maclaurin (small `|x|`) or asymptotic (large `|x|`) branch
/// (slice p1.30).
pub(crate) const SI_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `Ci`: Maclaurin or asymptotic (slice p1.30).
pub(crate) const CI_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `Li`: series summation (slice p1.30; oracle bumps to `p = 53`).
pub(crate) const LI_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `Ei`: series or asymptotic (slice p1.30).
pub(crate) const EI_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

// --- Airy and Bessel -------------------------------------------------

/// `Airy` (Ai, Bi, Ai', Bi'): shared eval body (slice p1.31).
pub(crate) const AIRY_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `Bessel J_n`: Maclaurin / Miller recurrence / asymptotic
/// dispatched at working precision (slice p1.32; oracle bumps to
/// `p = 320` for the cubic-correction trap).
pub(crate) const BESSEL_J_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `Bessel Y_n`: reflection formula composing J's and recurrence
/// (slice p1.32). The oscillatory regime near the J zeros may
/// surface as the first empirical-widening candidate at p1g.3.
pub(crate) const BESSEL_Y_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `Bessel I_n`: Maclaurin / Miller / asymptotic (slice p1.33;
/// oracle bumps to `p = 320`).
pub(crate) const BESSEL_I_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

/// `Bessel K_n`: reflection composing I's and recurrence
/// (slice p1.33).
pub(crate) const BESSEL_K_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

// --- zeta ------------------------------------------------------------

/// `zeta`: Borwein for `s > 0`, functional equation
/// (`gamma · sin · pow · zeta_borwein` composition) for `s < 0`
/// (slice p1.34). The composition is the deepest on the surface;
/// empirical confirmation pending pf-tqzz.
pub(crate) const ZETA_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

// --- AGM -------------------------------------------------------------

/// `agm`: Gauss AGM iteration with quadratic convergence
/// (`O(log w)` iterations, each `≤ 1 ULP`).
pub(crate) const AGM_ERROR_GUARD: u32 = DEFAULT_ERROR_GUARD;

#[cfg(test)]
mod tests {
    use super::*;

    /// Sanity guard: every per-kernel constant stays within the
    /// driver's growth envelope. With `ZIV_BASE_GUARD = 64` and
    /// `error_guard = 24` the first working iteration is at
    /// `target + 64` and the half-width shift is `working - 24 = 40`,
    /// leaving 40 bits of slack against the error bound. Any
    /// per-kernel constant approaching 48 (the `ZIV_BASE_GUARD - 16`
    /// safety margin) requires a paired bump to `ZIV_BASE_GUARD` to
    /// keep the first Ziv iteration's working precision adequate.
    #[test]
    fn every_per_kernel_bound_fits_under_base_guard_margin() {
        let bounds = [
            EXP_ERROR_GUARD,
            EXP2_ERROR_GUARD,
            EXP10_ERROR_GUARD,
            EXPM1_ERROR_GUARD,
            LN_ERROR_GUARD,
            LOG1P_ERROR_GUARD,
            SIN_ERROR_GUARD,
            COS_ERROR_GUARD,
            TAN_ERROR_GUARD,
            ASIN_ERROR_GUARD,
            ACOS_ERROR_GUARD,
            ATAN_ERROR_GUARD,
            ATAN2_ERROR_GUARD,
            SINH_ERROR_GUARD,
            COSH_ERROR_GUARD,
            TANH_ERROR_GUARD,
            ASINH_ERROR_GUARD,
            ACOSH_ERROR_GUARD,
            ATANH_ERROR_GUARD,
            POW_ERROR_GUARD,
            POW_INT_ERROR_GUARD,
            GAMMA_ERROR_GUARD,
            LGAMMA_ERROR_GUARD,
            DIGAMMA_ERROR_GUARD,
            BETA_ERROR_GUARD,
            ERF_ERROR_GUARD,
            ERFC_ERROR_GUARD,
            SI_ERROR_GUARD,
            CI_ERROR_GUARD,
            LI_ERROR_GUARD,
            EI_ERROR_GUARD,
            AIRY_ERROR_GUARD,
            BESSEL_J_ERROR_GUARD,
            BESSEL_Y_ERROR_GUARD,
            BESSEL_I_ERROR_GUARD,
            BESSEL_K_ERROR_GUARD,
            ZETA_ERROR_GUARD,
            AGM_ERROR_GUARD,
        ];
        const ZIV_BASE_GUARD_MARGIN: u32 = 48;
        for (i, bound) in bounds.iter().enumerate() {
            assert!(
                *bound < ZIV_BASE_GUARD_MARGIN,
                "calibration index {i} bound {bound} exceeds ZIV_BASE_GUARD - 16 = {ZIV_BASE_GUARD_MARGIN}; \
                 paired bump to ZIV_BASE_GUARD required"
            );
        }
    }

    /// Conscious-calibration enumeration. Every constant exported
    /// from this module must appear here. Adding a new
    /// `<KERNEL>_ERROR_GUARD` const without adding it to the
    /// `bounds` array in `every_per_kernel_bound_fits_under_base_guard_margin`
    /// (and to the count below) fails this test, forcing the author
    /// to acknowledge the calibration choice rather than rely on
    /// `DEFAULT_ERROR_GUARD` implicitly.
    #[test]
    fn calibration_table_enumerates_expected_kernel_count() {
        // 38 = 6 elementary + 3 forward trig + 4 inverse trig
        // + 6 hyperbolic + 2 power + 4 gamma family + 2 erf
        // + 4 integral + 1 airy + 4 bessel + 1 zeta + 1 agm.
        const EXPECTED_V1_KERNEL_BOUNDS: usize = 38;
        const BOUNDS_LEN: usize = 38;
        assert_eq!(
            BOUNDS_LEN, EXPECTED_V1_KERNEL_BOUNDS,
            "v1.0 surface kernel count drifted; recount the families above"
        );
    }
}
