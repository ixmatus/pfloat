//! The status-flag gate policy.
//!
//! VALUE bit-exactness is always the hard gate (see `verify`). On the
//! IEEE status flags this module decides the policy. Under the default
//! [`StatusGate::ValueAndDomainHard`], `INVALID` and `DIV_BY_ZERO` are
//! also hard, checked against an independent expectation; `INEXACT`,
//! `OVERFLOW`, and `UNDERFLOW` are not gated, because the directed-pair
//! shell conservatively over-reports `INEXACT` on composed-exact
//! results (`log10(1000)=3`, `exp10(2)=100`, `exp2(10)=1024`; pf-njs5).
//!
//! The two gated expectations are derived independently of the shell:
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

/// Check the gated flags. Returns the first disagreement as
/// `(flag, expected, got)`, or `None` when the gated flags agree (or
/// the gate is [`StatusGate::ValueOnly`]). Only called after VALUE has
/// matched.
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
