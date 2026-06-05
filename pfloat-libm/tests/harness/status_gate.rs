//! The status-flag gate policy.
//!
//! VALUE bit-exactness is always the hard gate (see `verify`). On the
//! IEEE status flags this module decides the policy. Under the default
//! [`StatusGate::ValueAndDomainHard`], `INVALID` and `DIV_BY_ZERO` are
//! hard, checked against an independent expectation, and `INEXACT` is
//! hard for the exp/log family and sin/cos (pf-njs5, ADR-0060), checked
//! against an expectation derived from the oracle enclosure rather than
//! the shell. `OVERFLOW` and `UNDERFLOW` remain ungated.
//!
//! pfloat's `INEXACT` over-reported on composed-exact results
//! (`log10(1000)=3`, `exp10(2)=100`, `exp2(10)=1024`) and under-reported
//! on sub-working-precision residuals (`exp(2^-1074)=1.0`); ADR-0060
//! closes both for the gated functions via a pre-Ziv exact-input
//! dispatch plus a transcendental force, so the flag is now reliable
//! enough to gate. The function set ([`inexact_is_gated`]) is exactly
//! the kernels fixed there; the rest stay ungated pending pf-uqd1.
//!
//! The gated expectations are derived independently of the shell:
//!
//! - `INVALID` comes straight from the oracle enclosure: a finite or
//!   infinite (non-NaN) input whose true value is NaN is a domain error
//!   (`ln` of a negative, `asin` out of range, trig of an infinity,
//!   ...). A NaN input producing NaN is propagation, not `INVALID`.
//!   This needs no per-function domain table; the enclosure being NaN
//!   *is* the domain-error witness.
//!
//! - `DIV_BY_ZERO` comes from a tiny exact-pole set. A pole raises the
//!   flag only at an exactly representable input (irrational poles such
//!   as `cot`/`sec` at `pi`/`pi/2` are never hit by a float, so the
//!   value there is finite). Deriving it from "enclosure is infinite"
//!   instead would misfire on the overflow regime, where MPFR's own
//!   exponent range overflows `exp(huge)` to infinity though the true
//!   value is finite. The exact-input list sidesteps that.

#![cfg(all(unix, feature = "differential-mpfr"))]

use pfloat_libm::Status;

use super::types::{Enclosure, FlagKind, LibmFnId};

/// How strictly the harness treats the IEEE status flags.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Default)]
pub enum StatusGate {
    /// VALUE + `INVALID` + `DIV_BY_ZERO` are hard; the rest are
    /// ungated. The default.
    #[default]
    ValueAndDomainHard,
    /// Only VALUE is hard; flags are not checked.
    ValueOnly,
}

/// `INVALID` is expected when a non-NaN input has a NaN true value (a
/// domain error). Both enclosure endpoints are NaN for such inputs.
pub fn expected_invalid(input_is_nan: bool, enc: &Enclosure) -> bool {
    !input_is_nan && enc.lo.is_nan() && enc.hi.is_nan()
}

/// `DIV_BY_ZERO` is expected at the exactly representable poles of the
/// v0.1 surface. `value` is the hardware input widened to `f64`
/// (`0.0 == -0.0` in Rust, so the `±0` poles match either sign).
pub fn expected_div_by_zero(f: LibmFnId, value: f64) -> bool {
    match f {
        // log family: pole at 0.
        LibmFnId::Ln | LibmFnId::Log2 | LibmFnId::Log10 => value == 0.0,
        // log1p(x) = ln(1 + x): pole at x = -1.
        LibmFnId::Log1p => value == -1.0,
        // atanh: poles at ±1.
        LibmFnId::Atanh => value == 1.0 || value == -1.0,
        // cot, csc: pole at 0 (the only representable zero of sin).
        // sec/tan/sin/cos have no representable pole.
        LibmFnId::Cot | LibmFnId::Csc => value == 0.0,
        // rootn(0, negative n) = 1 / 0^|n| = +inf.
        LibmFnId::Rootn(n) => value == 0.0 && n < 0,
        _ => false,
    }
}

/// Whether `INEXACT` is hard-gated for `f`. The set is the kernels whose
/// exact-output set over the dyadic inputs is fully characterized and
/// whose flag was corrected under ADR-0060 (the exp/log family and
/// sin/cos) and ADR-0063 (the rest of the elementary transcendental
/// surface: the remaining trig, the reciprocal trig, inverse trig,
/// hyperbolic, inverse hyperbolic, and expm1/log1p). Algebraic kernels
/// (sqrt/cbrt/hypot/rootn) are not transcendental and stay ungated.
///
/// The expectation itself is `INEXACT` ⇔ the true result is not exactly
/// representable in the hardware type, witnessed by the oracle
/// enclosure (not the shell). That test is width-generic, so it lives
/// at the `verify` call site where the `Hw` width and oracle precision
/// are in scope; this predicate is the policy half.
pub fn inexact_is_gated(f: LibmFnId) -> bool {
    matches!(
        f,
        // exp/log family + sin/cos (ADR-0060).
        LibmFnId::Exp
            | LibmFnId::Exp2
            | LibmFnId::Exp10
            | LibmFnId::Ln
            | LibmFnId::Log2
            | LibmFnId::Log10
            | LibmFnId::Sin
            | LibmFnId::Cos
            // The rest of the elementary transcendentals (ADR-0063).
            | LibmFnId::Tan
            | LibmFnId::Cot
            | LibmFnId::Sec
            | LibmFnId::Csc
            | LibmFnId::Asin
            | LibmFnId::Acos
            | LibmFnId::Atan
            | LibmFnId::Sinh
            | LibmFnId::Cosh
            | LibmFnId::Tanh
            | LibmFnId::Asinh
            | LibmFnId::Acosh
            | LibmFnId::Atanh
            | LibmFnId::Expm1
            | LibmFnId::Log1p
    )
}

/// Check the gated flags. Returns the first disagreement as
/// `(flag, expected, got)`, or `None` when the gated flags agree (or
/// the gate is [`StatusGate::ValueOnly`]). Only called after VALUE has
/// matched. `INEXACT` is gated separately at the `verify` call site
/// (see [`inexact_is_gated`]); it needs the hardware width and oracle
/// precision this function does not carry.
pub fn check_flags(
    gate: StatusGate,
    f: LibmFnId,
    input_is_nan: bool,
    value_f64: f64,
    enc: &Enclosure,
    got: Status,
) -> Option<(FlagKind, bool, bool)> {
    if let StatusGate::ValueOnly = gate {
        return None;
    }
    let exp_invalid = expected_invalid(input_is_nan, enc);
    if got.invalid() != exp_invalid {
        return Some((FlagKind::Invalid, exp_invalid, got.invalid()));
    }
    let exp_dbz = expected_div_by_zero(f, value_f64);
    if got.div_by_zero() != exp_dbz {
        return Some((FlagKind::DivByZero, exp_dbz, got.div_by_zero()));
    }
    None
}
