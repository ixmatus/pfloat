//! MPFR backend for the Oracle harness.
//!
//! Produces enclosures via two-mode evaluation through `rug`: one
//! call with `Round::Down` (`MPFR_RNDD`, toward minus infinity)
//! and one with `Round::Up` (`MPFR_RNDU`, toward plus infinity),
//! both at the requested working precision. The two directed
//! values provably bracket the mathematically true result by
//! MPFR's primitive contract for every function in the
//! MPFR-primary subset of the v1.0 surface; the bracket tightens
//! monotonically as the working precision rises.
//!
//! Functions without an MPFR primitive (`Si`, `Ci`, `li`, `Bi`,
//! `Ai_prime`, `Bi_prime`, the modified Bessel `I`/`K` family) are
//! Arb-primary per ADR-0034's per-function oracle routing; the
//! Arb backend lands in a follow-up slice and this backend returns
//! `unimplemented!()` for those `FnId` variants until then.
//!
//! See ADR-0034 for the design rationale and ADR-0014 for the
//! existing MPFR differential lane's gating posture this module
//! inherits.

#![cfg(all(unix, feature = "differential-mpfr"))]

use core::cmp::Ordering;

use rug::float::Round;
use rug::ops::{AssignRound, CompleteRound};
use rug::Float;

use super::types::{Enclosure, FnId, OracleBackend};

/// In-process MPFR-backed oracle. A unit struct: `rug` carries no
/// per-backend state, so one `MpfrOracle` value serves every call.
pub struct MpfrOracle;

/// Build an [`Enclosure`] by evaluating an `rug` Incomplete
/// expression at `prec` under `Round::Down` and `Round::Up`. Each
/// usage of `$method($($arg),*)` reconstructs a fresh Incomplete
/// (Incomplete consumes itself on assignment), so the macro emits
/// two independent directed-rounding evaluations.
macro_rules! bracket {
    ($x:expr, $prec:expr, $method:ident $(, $arg:expr)*) => {
        Enclosure {
            lo: Float::with_val_round($prec, $x.$method($($arg),*), Round::Down).0,
            hi: Float::with_val_round($prec, $x.$method($($arg),*), Round::Up).0,
        }
    };
}

impl MpfrOracle {
    /// Mode-independent rigorous-enclosure midpoint of the function
    /// at `input` (binary32 bit pattern), at `oracle_prec` precision.
    /// Mirrors [`super::arb::ArbOracle::midpoint`] for the 35
    /// MPFR-primary `FnId`s; the cross-check harness routes each
    /// `FnId` to the appropriate backend.
    ///
    /// MPFR's correct-rounding under `Round::Nearest` at
    /// `oracle_prec >= working_prec + 64` produces a midpoint
    /// within sub-ULP of the true value (the round-to-nearest
    /// closest representable at `oracle_prec`), comfortably below
    /// the pf-tqzz cross-check tolerance
    /// `2^(error_guard - working_prec) * |midpoint|`.
    ///
    /// pf-tqzz (slice p1g.3, ADR-0039). Panics on an Arb-primary
    /// `FnId`; the caller is responsible for routing those through
    /// [`super::arb::ArbOracle::midpoint`].
    #[must_use]
    #[allow(clippy::unused_self)]
    pub fn midpoint(&self, f: FnId, input: u32, oracle_prec: u32) -> Float {
        let x = Float::with_val(oracle_prec, f32::from_bits(input));
        match f {
            // Elementary.
            FnId::Sqrt => Float::with_val(oracle_prec, x.sqrt_ref()),
            FnId::Cbrt => Float::with_val(oracle_prec, x.cbrt_ref()),
            FnId::Exp => Float::with_val(oracle_prec, x.exp_ref()),
            FnId::Exp2 => Float::with_val(oracle_prec, x.exp2_ref()),
            FnId::Exp10 => Float::with_val(oracle_prec, x.exp10_ref()),
            FnId::Expm1 => Float::with_val(oracle_prec, x.exp_m1_ref()),
            FnId::Ln => Float::with_val(oracle_prec, x.ln_ref()),
            FnId::Log1p => Float::with_val(oracle_prec, x.ln_1p_ref()),
            FnId::Log2 => Float::with_val(oracle_prec, x.log2_ref()),
            FnId::Log10 => Float::with_val(oracle_prec, x.log10_ref()),
            FnId::Sin => Float::with_val(oracle_prec, x.sin_ref()),
            FnId::Cos => Float::with_val(oracle_prec, x.cos_ref()),
            FnId::Tan => Float::with_val(oracle_prec, x.tan_ref()),
            FnId::Asin => Float::with_val(oracle_prec, x.asin_ref()),
            FnId::Acos => Float::with_val(oracle_prec, x.acos_ref()),
            FnId::Atan => Float::with_val(oracle_prec, x.atan_ref()),
            FnId::Sinh => Float::with_val(oracle_prec, x.sinh_ref()),
            FnId::Cosh => Float::with_val(oracle_prec, x.cosh_ref()),
            FnId::Tanh => Float::with_val(oracle_prec, x.tanh_ref()),
            FnId::Asinh => Float::with_val(oracle_prec, x.asinh_ref()),
            FnId::Acosh => Float::with_val(oracle_prec, x.acosh_ref()),
            FnId::Atanh => Float::with_val(oracle_prec, x.atanh_ref()),
            // Specials.
            FnId::Erf => Float::with_val(oracle_prec, x.erf_ref()),
            FnId::Erfc => Float::with_val(oracle_prec, x.erfc_ref()),
            FnId::Gamma => Float::with_val(oracle_prec, x.gamma_ref()),
            FnId::Lgamma => {
                // MPFR's `ln_abs_gamma_ref` returns an Incomplete
                // that assigns to `(&mut Float, &mut Ordering)`;
                // the Ordering slot is the Γ sign byproduct, which
                // the cross-check ignores (the magnitude is the
                // load-bearing piece). Mirrors `lgamma_bracket`
                // below.
                let mut val = Float::new(oracle_prec);
                let mut sign_unused = Ordering::Equal;
                let _ =
                    (&mut val, &mut sign_unused).assign_round(x.ln_abs_gamma_ref(), Round::Nearest);
                let _ = sign_unused;
                val
            }
            FnId::Digamma => Float::with_val(oracle_prec, x.digamma_ref()),
            FnId::Zeta => Float::with_val(oracle_prec, x.zeta_ref()),
            FnId::Ei => Float::with_val(oracle_prec, x.eint_ref()),

            // Airy: only `Ai` has an MPFR primitive (the others route
            // through Arb), mirroring the `enclose` dispatch above.
            FnId::Ai => Float::with_val(oracle_prec, x.ai_ref()),

            // Bessel J / Y: every fixed-order variant has an MPFR
            // primitive; the parametric-order variants take the extra
            // `n: i32` argument. These share the `enclose` primitives
            // (`j0_ref` … `yn_ref`); the `midpoint` verb was missing
            // them, which panicked the J/Y cross-check shards (pf-ypfl,
            // discovered-from pf-hcz4).
            FnId::BesselJ0 => Float::with_val(oracle_prec, x.j0_ref()),
            FnId::BesselJ1 => Float::with_val(oracle_prec, x.j1_ref()),
            FnId::BesselJn(n) => Float::with_val(oracle_prec, x.jn_ref(n)),
            FnId::BesselY0 => Float::with_val(oracle_prec, x.y0_ref()),
            FnId::BesselY1 => Float::with_val(oracle_prec, x.y1_ref()),
            FnId::BesselYn(n) => Float::with_val(oracle_prec, x.yn_ref(n)),

            // Genuinely Arb-primary: no MPFR primitive. The caller
            // routes these through `ArbOracle::midpoint`; reaching here
            // is a routing bug, so panic with a precise message.
            FnId::Si
            | FnId::Ci
            | FnId::Li
            | FnId::Bi
            | FnId::AiPrime
            | FnId::BiPrime
            | FnId::BesselI0
            | FnId::BesselI1
            | FnId::BesselIn(_)
            | FnId::BesselK0
            | FnId::BesselK1
            | FnId::BesselKn(_)
            | FnId::Cot
            | FnId::Sec
            | FnId::Csc => {
                panic!("MpfrOracle::midpoint called with Arb-primary FnId: {f:?}; route via ArbOracle::midpoint")
            }
        }
    }
}

impl OracleBackend for MpfrOracle {
    // `mode` is unused: the MPFR backend produces a mode-agnostic
    // bracket via directed Round::Down and Round::Up evaluations;
    // the verifier rounds the bracket to the caller's mode itself.
    // The parameter exists only for trait-signature parity with
    // the ADR-0035 Arb backend.
    fn enclose(
        &self,
        f: FnId,
        input: u32,
        _mode: pfloat::RoundingMode,
        working_prec: u32,
    ) -> Enclosure {
        // Convert the f32 bit pattern exactly into a high-precision
        // Float; an f32 is representable without rounding at any
        // precision >= 24, so working_prec >= 64 (Phase 1's
        // START_PREC) is comfortably exact.
        let x = Float::with_val(working_prec, f32::from_bits(input));
        match f {
            // Elementary.
            FnId::Sqrt => bracket!(x, working_prec, sqrt_ref),
            FnId::Cbrt => bracket!(x, working_prec, cbrt_ref),
            FnId::Exp => bracket!(x, working_prec, exp_ref),
            FnId::Exp2 => bracket!(x, working_prec, exp2_ref),
            FnId::Exp10 => bracket!(x, working_prec, exp10_ref),
            FnId::Expm1 => bracket!(x, working_prec, exp_m1_ref),
            FnId::Ln => bracket!(x, working_prec, ln_ref),
            FnId::Log1p => bracket!(x, working_prec, ln_1p_ref),
            FnId::Log2 => bracket!(x, working_prec, log2_ref),
            FnId::Log10 => bracket!(x, working_prec, log10_ref),
            FnId::Sin => bracket!(x, working_prec, sin_ref),
            FnId::Cos => bracket!(x, working_prec, cos_ref),
            FnId::Tan => bracket!(x, working_prec, tan_ref),
            FnId::Asin => bracket!(x, working_prec, asin_ref),
            FnId::Acos => bracket!(x, working_prec, acos_ref),
            FnId::Atan => bracket!(x, working_prec, atan_ref),
            FnId::Sinh => bracket!(x, working_prec, sinh_ref),
            FnId::Cosh => bracket!(x, working_prec, cosh_ref),
            FnId::Tanh => bracket!(x, working_prec, tanh_ref),
            FnId::Asinh => bracket!(x, working_prec, asinh_ref),
            FnId::Acosh => bracket!(x, working_prec, acosh_ref),
            FnId::Atanh => bracket!(x, working_prec, atanh_ref),

            // Specials.
            FnId::Erf => bracket!(x, working_prec, erf_ref),
            FnId::Erfc => bracket!(x, working_prec, erfc_ref),
            FnId::Gamma => bracket!(x, working_prec, gamma_ref),
            FnId::Lgamma => lgamma_bracket(&x, working_prec),
            FnId::Digamma => bracket!(x, working_prec, digamma_ref),
            FnId::Zeta => bracket!(x, working_prec, zeta_ref),
            FnId::Ei => bracket!(x, working_prec, eint_ref),

            // Airy (only `Ai` has an MPFR primitive; the others
            // route through Arb).
            FnId::Ai => bracket!(x, working_prec, ai_ref),

            // Bessel J / Y (every fixed-order variant has an MPFR
            // primitive; the parametric-order variants take an
            // additional `n: i32` argument).
            FnId::BesselJ0 => bracket!(x, working_prec, j0_ref),
            FnId::BesselJ1 => bracket!(x, working_prec, j1_ref),
            FnId::BesselJn(n) => bracket!(x, working_prec, jn_ref, n),
            FnId::BesselY0 => bracket!(x, working_prec, y0_ref),
            FnId::BesselY1 => bracket!(x, working_prec, y1_ref),
            FnId::BesselYn(n) => bracket!(x, working_prec, yn_ref, n),

            // Arb-primary (no MPFR primitive). Land with the Arb
            // backend in a follow-up slice.
            FnId::Si
            | FnId::Ci
            | FnId::Li
            | FnId::Bi
            | FnId::AiPrime
            | FnId::BiPrime
            | FnId::BesselI0
            | FnId::BesselI1
            | FnId::BesselIn(_)
            | FnId::BesselK0
            | FnId::BesselK1
            | FnId::BesselKn(_)
            | FnId::Cot
            | FnId::Sec
            | FnId::Csc => {
                unimplemented!("FnId::{f:?} requires the Arb backend")
            }
        }
    }

    fn name(&self) -> &'static str {
        "MPFR"
    }
}

/// `lgamma(x) = ln |Γ(x)|` bracket via MPFR's `lgamma` primitive.
///
/// `Float::ln_abs_gamma_ref` produces an Incomplete that assigns to
/// `(&mut Float, &mut Ordering)`; the second slot is `Γ(x)`'s sign
/// which pfloat's `lgamma` does not expose, so it is captured and
/// discarded. The `Round` parameter governs the value's rounding,
/// applied at the assignment step (the same directed-rounding
/// guarantee the `bracket!` macro relies on for the simpler
/// `with_val_round` callers).
fn lgamma_bracket(x: &Float, working_prec: u32) -> Enclosure {
    let mut lo = Float::new(working_prec);
    let mut sign_lo = Ordering::Equal;
    let _ = (&mut lo, &mut sign_lo).assign_round(x.ln_abs_gamma_ref(), Round::Down);

    let mut hi = Float::new(working_prec);
    let mut sign_hi = Ordering::Equal;
    let _ = (&mut hi, &mut sign_hi).assign_round(x.ln_abs_gamma_ref(), Round::Up);

    Enclosure { lo, hi }
}
