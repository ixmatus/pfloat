//! The MPFR enclosure oracle.
//!
//! For one `LibmFnId` and one input value (a `rug::Float` holding the
//! hardware input exactly), [`enclose`] returns the directed pair
//! `[f(x) toward -inf, f(x) toward +inf]` at the requested working
//! precision. MPFR guarantees correct rounding for every primitive used
//! here, so the pair brackets the true value, and the bracket tightens
//! as the precision grows (the Ziv-at-oracle argument in `verify`).
//!
//! rug 1.30 (MPFR 4.2.2) exposes a primitive for every v0.1 function,
//! including cot/sec/csc (`cot_ref`/`sec_ref`/`csc_ref`), cbrt, hypot,
//! and the IEEE 754-2019 n-th root (`root_ref` for n > 0, `root_i_ref`
//! for signed n). So the harness is MPFR-only, with no Arb or
//! python-flint dependency. This corrects pfloat's own oracle, which
//! routes cot/sec/csc through Arb on the (mistaken) belief that MPFR
//! lacks them (ADR-0058).

#![cfg(all(unix, feature = "differential-mpfr"))]

use rug::float::{Round, Special};
use rug::Float;

use super::types::{Enclosure, LibmFnId};

/// Evaluate `$x.$m(args)` rounded toward `-inf` and toward `+inf` at
/// precision `$p`, returning the bracket.
macro_rules! bracket {
    ($x:expr, $p:expr, $m:ident $(, $arg:expr)*) => {
        Enclosure {
            lo: Float::with_val_round($p, $x.$m($($arg),*), Round::Down).0,
            hi: Float::with_val_round($p, $x.$m($($arg),*), Round::Up).0,
        }
    };
}

/// A bracket whose endpoints are both NaN. Used for domain errors MPFR
/// cannot express directly (`rootn` with order 0); the verifier reads a
/// both-NaN bracket as a certified NaN result and the gate derives an
/// expected `INVALID`.
fn nan_bracket(p: u32) -> Enclosure {
    Enclosure {
        lo: Float::with_val(p, Special::Nan),
        hi: Float::with_val(p, Special::Nan),
    }
}

/// Enclose `f(input)` at `prec` bits. `input` is the exact `rug::Float`
/// of the hardware argument; `partner` supplies `hypot`'s second
/// coordinate (already lifted exactly at `prec`).
pub fn enclose(f: LibmFnId, input: &Float, partner: Option<&Float>, prec: u32) -> Enclosure {
    match f {
        LibmFnId::Exp => bracket!(input, prec, exp_ref),
        LibmFnId::Exp2 => bracket!(input, prec, exp2_ref),
        LibmFnId::Exp10 => bracket!(input, prec, exp10_ref),
        LibmFnId::Expm1 => bracket!(input, prec, exp_m1_ref),
        LibmFnId::Ln => bracket!(input, prec, ln_ref),
        LibmFnId::Log2 => bracket!(input, prec, log2_ref),
        LibmFnId::Log10 => bracket!(input, prec, log10_ref),
        LibmFnId::Log1p => bracket!(input, prec, ln_1p_ref),
        LibmFnId::Sqrt => bracket!(input, prec, sqrt_ref),
        LibmFnId::Cbrt => bracket!(input, prec, cbrt_ref),
        LibmFnId::Sin => bracket!(input, prec, sin_ref),
        LibmFnId::Cos => bracket!(input, prec, cos_ref),
        LibmFnId::Tan => bracket!(input, prec, tan_ref),
        LibmFnId::Cot => bracket!(input, prec, cot_ref),
        LibmFnId::Sec => bracket!(input, prec, sec_ref),
        LibmFnId::Csc => bracket!(input, prec, csc_ref),
        LibmFnId::Asin => bracket!(input, prec, asin_ref),
        LibmFnId::Acos => bracket!(input, prec, acos_ref),
        LibmFnId::Atan => bracket!(input, prec, atan_ref),
        LibmFnId::Sinh => bracket!(input, prec, sinh_ref),
        LibmFnId::Cosh => bracket!(input, prec, cosh_ref),
        LibmFnId::Tanh => bracket!(input, prec, tanh_ref),
        LibmFnId::Asinh => bracket!(input, prec, asinh_ref),
        LibmFnId::Acosh => bracket!(input, prec, acosh_ref),
        LibmFnId::Atanh => bracket!(input, prec, atanh_ref),
        LibmFnId::Hypot => {
            let y = partner.expect("hypot requires a partner operand");
            bracket!(input, prec, hypot_ref, y)
        }
        LibmFnId::Rootn(n) => {
            if n == 0 {
                // IEEE 754-2019 §9.2: rootn(x, 0) is a domain error.
                nan_bracket(prec)
            } else if n > 0 {
                bracket!(input, prec, root_ref, n as u32)
            } else {
                bracket!(input, prec, root_i_ref, n)
            }
        }
    }
}

/// Oracle name for the status table's `oracle` column.
pub const ORACLE_NAME: &str = "MPFR";
