//! Core types of the Oracle harness.
//!
//! [`Enclosure`] is a proven bracket `[lo, hi]` of `f(input)` at the
//! backend's working precision. [`OracleBackend`] wraps any backend
//! that can produce enclosures. [`FnId`] identifies which function
//! the harness is verifying; backends match on it to dispatch to
//! the corresponding rug or python-flint call. [`Verdict`] is the
//! per-input result the verifier accumulates.
//!
//! See ADR-0034 for the architectural rationale and the alternatives
//! considered.

#![cfg(all(unix, feature = "differential-mpfr"))]

use pfloat::RoundingMode;
use rug::Float;

/// Proven bracket of the true value: `lo <= true <= hi`. Endpoints
/// are produced at the backend's working precision; the bracket
/// tightens as the working precision rises.
///
/// For the MPFR backend, `lo` is `f(input)` evaluated with
/// `MPFR_RNDD` (toward minus infinity) and `hi` is `f(input)`
/// evaluated with `MPFR_RNDU` (toward plus infinity). The two
/// directed values provably bracket the mathematically true result;
/// MPFR guarantees this for every primitive it ships.
///
/// For the Arb backend (next slice) the ball type already carries
/// midpoint plus radius natively, so the `(lo, hi)` shape is a
/// direct read of `(mid - rad, mid + rad)` rather than a two-mode
/// evaluation.
#[derive(Clone, Debug)]
pub struct Enclosure {
    pub lo: Float,
    pub hi: Float,
}

/// The outcome of an [`OracleBackend::enclose`] call.
///
/// Two outcomes must stay distinct, and conflating them is exactly
/// the honesty bug pf-41ou fixed:
///
/// - [`Self::Bracket`] carries a proven bracket. Its endpoints may be
///   finite, infinite, or NaN. NaN endpoints certify a *genuinely
///   NaN* true value (`ln(x)` for `x < 0`, `Ci(x)` for `x < 0`, a NaN
///   input): the oracle *knows* the answer is NaN, and the verifier
///   accepts a matching NaN from the kernel.
/// - [`Self::Inconclusive`] is the backend saying "I could not
///   certify a unique `f32`." An authoritative backend (the Arb
///   worker) returns this when its own internal precision loop
///   reaches its cap without the bracket collapsing to one `f32`
///   (the worker's `INC` reply), and the `differential-arb`-absent
///   fallback returns it because no oracle exists for the `FnId`.
///   The verifier maps it to [`Verdict::OracleInconclusive`], never
///   to `Ok`.
///
/// The old code encoded both as an [`Enclosure`] with NaN endpoints,
/// so an `INC` reply was read as a certified NaN and silently counted
/// as agreement whenever the kernel also returned NaN. The sum type
/// makes that state unrepresentable: `Inconclusive` is not a bracket
/// at all, so no rounding of it can ever certify.
#[derive(Clone, Debug)]
pub enum Enclosed {
    /// A proven bracket `[lo, hi]` of the true value.
    Bracket(Enclosure),
    /// The backend could not certify a unique `f32` (the worker's
    /// `INC` reply, or the absent-oracle fallback). Distinct from a
    /// [`Self::Bracket`] with NaN endpoints (a certified NaN value).
    Inconclusive,
}

/// A backend that can enclose a function value at a requested
/// working precision.
///
/// The trait is `Send` so the harness's per-function driver can
/// fan shards across cores in a follow-up slice. The methods are
/// `&self` so a single backend instance serves the entire sweep.
pub trait OracleBackend: Send + Sync {
    /// Enclose `f(input)` at `working_prec` bits under the caller's
    /// rounding mode. The returned bracket's endpoints are exact
    /// at the working precision.
    ///
    /// `mode` is a hint that authoritative backends (see
    /// [`OracleBackend::is_authoritative`]) use to produce a
    /// single-point enclosure at the mode's certified `f32`. Non
    /// authoritative backends ignore `mode`; the verifier does the
    /// mode-specific rounding from the returned bracket itself.
    /// Added at slice p1.8 / ADR-0035 to support the
    /// worker-reports-certified-f32-directly protocol while keeping
    /// the trait shape uniform for the older `MpfrOracle` enclose
    /// path.
    ///
    /// Returns [`Enclosed::Inconclusive`] when the backend cannot
    /// certify a unique `f32` (an authoritative worker's `INC` reply,
    /// or the absent-oracle fallback). This is deliberately distinct
    /// from an [`Enclosed::Bracket`] whose endpoints are NaN, which
    /// certifies a genuinely-NaN true value. The verifier routes the
    /// two differently: a bracket is rounded and compared, an
    /// inconclusive outcome becomes [`Verdict::OracleInconclusive`]
    /// (pf-41ou).
    fn enclose(&self, f: FnId, input: u32, mode: RoundingMode, working_prec: u32) -> Enclosed;

    /// Backend identifier for the status table's `oracle` column
    /// (e.g. `"MPFR"`, `"Arb"`).
    fn name(&self) -> &'static str;

    /// `true` when [`Self::enclose`] already runs the backend's own
    /// internal precision loop and further calls at higher
    /// `working_prec` would not refine the answer. The verifier
    /// short-circuits its Ziv-at-oracle loop for authoritative
    /// backends and calls `enclose` exactly once.
    ///
    /// The default returns `false`, matching the historical
    /// `MpfrOracle` shape where the bracket tightens as
    /// `working_prec` grows. The `ArbOracle` overrides this to
    /// return `true` under ADR-0035: the worker's own Ziv loop
    /// inside the subprocess decides the certified `f32` directly,
    /// and the Rust-side verifier just receives the answer.
    fn is_authoritative(&self) -> bool {
        false
    }
}

/// Identifier for each function the Phase 1 sweep verifies.
///
/// Matches the frozen unary surface enumerated in
/// `docs/v1.0-surface.md`. Backends dispatch on this enum to the
/// corresponding `rug` (MPFR backend) or `python-flint` (Arb
/// backend) call. The variants with an integer argument
/// (`Jn(order)`, `Yn(order)`, `In(order)`, `Kn(order)`) verify a
/// specific Bessel order; the harness iterates over a chosen set
/// of orders per the Phase 1 plan and records `N` in the status
/// table per function family.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum FnId {
    // Elementary.
    Sqrt,
    Cbrt,
    Exp,
    Exp2,
    Exp10,
    Expm1,
    Ln,
    Log1p,
    Log2,
    Log10,
    Sin,
    Cos,
    Tan,
    Cot,
    Sec,
    Csc,
    Asin,
    Acos,
    Atan,
    Sinh,
    Cosh,
    Tanh,
    Asinh,
    Acosh,
    Atanh,
    // Specials.
    Erf,
    Erfc,
    Gamma,
    Lgamma,
    Digamma,
    Zeta,
    Ei,
    Si,
    Ci,
    Li,
    // Airy.
    Ai,
    Bi,
    AiPrime,
    BiPrime,
    // Bessel J/Y (fixed and parametric orders).
    BesselJ0,
    BesselJ1,
    BesselJn(i32),
    BesselY0,
    BesselY1,
    BesselYn(i32),
    // Bessel I/K (fixed and parametric orders).
    BesselI0,
    BesselI1,
    BesselIn(i32),
    BesselK0,
    BesselK1,
    BesselKn(i32),
}

impl FnId {
    /// Display name for the status table's `function` column. The
    /// parametric variants name their order explicitly (e.g.
    /// `"Jn(5)"`).
    pub fn name(&self) -> &'static str {
        match self {
            Self::Sqrt => "sqrt",
            Self::Cbrt => "cbrt",
            Self::Exp => "exp",
            Self::Exp2 => "exp2",
            Self::Exp10 => "exp10",
            Self::Expm1 => "expm1",
            Self::Ln => "ln",
            Self::Log1p => "log1p",
            Self::Log2 => "log2",
            Self::Log10 => "log10",
            Self::Sin => "sin",
            Self::Cos => "cos",
            Self::Tan => "tan",
            Self::Cot => "cot",
            Self::Sec => "sec",
            Self::Csc => "csc",
            Self::Asin => "asin",
            Self::Acos => "acos",
            Self::Atan => "atan",
            Self::Sinh => "sinh",
            Self::Cosh => "cosh",
            Self::Tanh => "tanh",
            Self::Asinh => "asinh",
            Self::Acosh => "acosh",
            Self::Atanh => "atanh",
            Self::Erf => "erf",
            Self::Erfc => "erfc",
            Self::Gamma => "gamma",
            Self::Lgamma => "lgamma",
            Self::Digamma => "digamma",
            Self::Zeta => "zeta",
            Self::Ei => "Ei",
            Self::Si => "Si",
            Self::Ci => "Ci",
            Self::Li => "li",
            Self::Ai => "Ai",
            Self::Bi => "Bi",
            Self::AiPrime => "Ai_prime",
            Self::BiPrime => "Bi_prime",
            Self::BesselJ0 => "J0",
            Self::BesselJ1 => "J1",
            Self::BesselJn(_) => "Jn",
            Self::BesselY0 => "Y0",
            Self::BesselY1 => "Y1",
            Self::BesselYn(_) => "Yn",
            Self::BesselI0 => "I0",
            Self::BesselI1 => "I1",
            Self::BesselIn(_) => "In",
            Self::BesselK0 => "K0",
            Self::BesselK1 => "K1",
            Self::BesselKn(_) => "Kn",
        }
    }
}

/// Per-input verdict the harness accumulates over the sweep.
///
/// `Ok` means pfloat's `f32`-rounded output matched the oracle's
/// certified `f32` under `mode`. `Mismatch` carries the exact bits
/// of the disagreement; the harness appends these to the
/// per-function regression corpus.
/// `OracleInconclusive` means the oracle could not certify the
/// rounding before the precision cap (`MAX_PREC`); these go to a
/// worst-case-candidate file rather than the failure corpus.
/// `Panic` carries the panic message captured via
/// `std::panic::catch_unwind`; these go to a panic-regression file
/// that runs on every CI push regardless of the next per-release
/// sweep schedule.
#[derive(Clone, Debug)]
pub enum Verdict {
    Ok,
    Mismatch {
        input: u32,
        mode: RoundingMode,
        expected: u32,
        got: u32,
    },
    OracleInconclusive {
        input: u32,
        mode: RoundingMode,
    },
    Panic {
        input: u32,
        mode: RoundingMode,
        message: String,
    },
}
