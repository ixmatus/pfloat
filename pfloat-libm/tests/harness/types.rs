//! Core types of the pfloat-libm verification harness.
//!
//! [`Enclosure`] is a proven bracket `[lo, hi]` of `f(input)` produced
//! by MPFR at a working precision (the `lo` end rounded toward minus
//! infinity, the `hi` end toward plus infinity). [`LibmFnId`] names the
//! function under test; the oracle and the shell dispatch on it.
//! [`Verdict`] is the per-input result the driver accumulates.
//!
//! This mirrors `pfloat/tests/oracle/types.rs`, reimplemented here
//! because the dependency direction (`pfloat-libm -> pfloat`) forbids
//! importing pfloat's test code. The libm version is simpler: the shell
//! already returns a hardware float, so there is no `BigFloat` bridge,
//! and the surface is the 27 v0.1 functions rather than pfloat's full
//! special-function set. See ADR-0058.

#![cfg(all(unix, feature = "differential-mpfr"))]

use pfloat_libm::RoundingMode;
use rug::Float;

/// Proven bracket of the true value: `lo <= f(x) <= hi`. The endpoints
/// are `f(x)` evaluated with `MPFR_RNDD` (toward minus infinity) and
/// `MPFR_RNDU` (toward plus infinity) at the oracle's working
/// precision; MPFR guarantees the bracket for every primitive it ships,
/// and the bracket tightens as the precision rises.
#[derive(Clone, Debug)]
pub struct Enclosure {
    pub lo: Float,
    pub hi: Float,
}

/// Hardware width a sweep targets. Recorded in the status row so a
/// reader can tell an `f32` exhaustive run from an `f64` differential.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Width {
    F32,
    F64,
}

impl Width {
    pub fn name(self) -> &'static str {
        match self {
            Width::F32 => "f32",
            Width::F64 => "f64",
        }
    }
}

/// The 27 functions the v0.1 shell exposes. The 25 unary variants are
/// unit; the two binary variants carry the information the oracle and
/// the shell need to rebuild the call from the `LibmFnId` alone:
/// `Rootn` carries its integer order, and `Hypot`'s partner travels
/// out of band in [`LibmArg`] (so the unary sweep machinery, which
/// iterates one input axis, can drive a binary function against a
/// fixed partner).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum LibmFnId {
    // exponentials
    Exp,
    Exp2,
    Exp10,
    Expm1,
    // logarithms
    Ln,
    Log2,
    Log10,
    Log1p,
    // roots
    Sqrt,
    Cbrt,
    // circular
    Sin,
    Cos,
    Tan,
    Cot,
    Sec,
    Csc,
    // inverse circular
    Asin,
    Acos,
    Atan,
    // hyperbolic
    Sinh,
    Cosh,
    Tanh,
    // inverse hyperbolic
    Asinh,
    Acosh,
    Atanh,
    // binary
    Hypot,
    Rootn(i32),
}

/// Out-of-band second operand for the binary functions, carried
/// alongside a swept input so the unary sweep machinery can drive a
/// binary function against a fixed partner. `Hypot`'s partner `y` is
/// stored as its hardware bit pattern (f32 in the low 32 bits, f64 in
/// all 64). `Rootn` carries its order in the variant, so its arg is
/// [`LibmArg::None`].
#[derive(Clone, Copy, Debug)]
pub enum LibmArg {
    None,
    HypotY(u64),
}

impl LibmFnId {
    /// Canonical name for the status table's `function` column.
    pub fn name(self) -> &'static str {
        match self {
            LibmFnId::Exp => "exp",
            LibmFnId::Exp2 => "exp2",
            LibmFnId::Exp10 => "exp10",
            LibmFnId::Expm1 => "expm1",
            LibmFnId::Ln => "ln",
            LibmFnId::Log2 => "log2",
            LibmFnId::Log10 => "log10",
            LibmFnId::Log1p => "log1p",
            LibmFnId::Sqrt => "sqrt",
            LibmFnId::Cbrt => "cbrt",
            LibmFnId::Sin => "sin",
            LibmFnId::Cos => "cos",
            LibmFnId::Tan => "tan",
            LibmFnId::Cot => "cot",
            LibmFnId::Sec => "sec",
            LibmFnId::Csc => "csc",
            LibmFnId::Asin => "asin",
            LibmFnId::Acos => "acos",
            LibmFnId::Atan => "atan",
            LibmFnId::Sinh => "sinh",
            LibmFnId::Cosh => "cosh",
            LibmFnId::Tanh => "tanh",
            LibmFnId::Asinh => "asinh",
            LibmFnId::Acosh => "acosh",
            LibmFnId::Atanh => "atanh",
            LibmFnId::Hypot => "hypot",
            LibmFnId::Rootn(_) => "rootn",
        }
    }

    pub fn is_binary(self) -> bool {
        matches!(self, LibmFnId::Hypot | LibmFnId::Rootn(_))
    }

    /// The 25 unary functions: the exhaustive-`f32`-sweep set. Binary
    /// functions are excluded (they cannot be exhausted over one axis).
    pub const UNARY: &'static [LibmFnId] = &[
        LibmFnId::Exp,
        LibmFnId::Exp2,
        LibmFnId::Exp10,
        LibmFnId::Expm1,
        LibmFnId::Ln,
        LibmFnId::Log2,
        LibmFnId::Log10,
        LibmFnId::Log1p,
        LibmFnId::Sqrt,
        LibmFnId::Cbrt,
        LibmFnId::Sin,
        LibmFnId::Cos,
        LibmFnId::Tan,
        LibmFnId::Cot,
        LibmFnId::Sec,
        LibmFnId::Csc,
        LibmFnId::Asin,
        LibmFnId::Acos,
        LibmFnId::Atan,
        LibmFnId::Sinh,
        LibmFnId::Cosh,
        LibmFnId::Tanh,
        LibmFnId::Asinh,
        LibmFnId::Acosh,
        LibmFnId::Atanh,
    ];
}

/// Which IEEE 754 status flag a [`Verdict::FlagMismatch`] is about. The
/// harness gates `INVALID` and `DIV_BY_ZERO` hard (an independent
/// expectation, see `status_gate`). `INEXACT` is gated for the exp/log
/// family and sin/cos (pf-njs5, ADR-0060) against an enclosure-derived
/// expectation; `OVERFLOW`/`UNDERFLOW` remain ungated.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FlagKind {
    Invalid,
    DivByZero,
    Inexact,
}

impl FlagKind {
    pub fn name(self) -> &'static str {
        match self {
            FlagKind::Invalid => "INVALID",
            FlagKind::DivByZero => "DIV_BY_ZERO",
            FlagKind::Inexact => "INEXACT",
        }
    }
}

/// Per-input verdict the driver accumulates. The hardware bit patterns
/// are widened to `u64` (f32 occupies the low 32 bits) so one type
/// serves both widths; the driver tags the width.
///
/// `ValueMismatch` is the headline failure: the shell's bits disagree
/// with the oracle's certified bits. `FlagMismatch` is a hard
/// `INVALID`/`DIV_BY_ZERO` disagreement. `OracleInconclusive` means the
/// enclosure still straddled an `f{32,64}` rounding boundary at the
/// precision cap (a measure-zero hard-to-round input); it is recorded,
/// not a failure, mirroring pfloat's own posture.
#[derive(Clone, Debug)]
pub enum Verdict {
    Ok,
    ValueMismatch {
        input: u64,
        mode: RoundingMode,
        expected: u64,
        got: u64,
    },
    FlagMismatch {
        input: u64,
        mode: RoundingMode,
        flag: FlagKind,
        expected: bool,
        got: bool,
    },
    OracleInconclusive {
        input: u64,
        mode: RoundingMode,
    },
}
