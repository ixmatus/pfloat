//! Verification core: Ziv-at-oracle precision loop.
//!
//! Given an [`OracleBackend`] and a function under test, the
//! verifier rounds the oracle's enclosure endpoints to `f32` under
//! the caller's rounding mode. If both endpoints round to the same
//! `f32` value, every point in the bracket (including the true
//! mathematical value) rounds there too and the correctly-rounded
//! `f32` is determined. If they straddle a rounding boundary, the
//! oracle's working precision doubles and the bracket tightens;
//! this is the Ziv interval-test argument from ADR-0022 applied to
//! the oracle's evaluation rather than to a pfloat working-precision
//! pass.
//!
//! The loop starts at `START_PREC = 64` (comfortably above f32's
//! 24 mantissa bits) and caps at `MAX_PREC = 1024`. Inputs that
//! exhaust the cap return [`Verdict::OracleInconclusive`]; they are
//! captured to a worst-case-candidate file by the per-function
//! driver, not the failure corpus.
//!
//! The verifier does not invoke a pfloat kernel directly. The
//! [`Kernel`] callback is the seam between the oracle side and
//! the pfloat-side dispatch table; the per-function driver supplies
//! the closure that calls pfloat's actual `FnId`-corresponding
//! kernel.

#![cfg(all(unix, feature = "differential-mpfr"))]

use pfloat::RoundingMode;

use super::convert::round_f32;
use super::types::{Enclosure, FnId, OracleBackend, Verdict};

/// First Ziv-at-oracle working precision. The oracle evaluates its
/// first enclosure at this precision; common-case inputs (mode-
/// uniform overflow / underflow, non-hard-to-round normals) certify
/// here without further iteration.
pub const START_PREC: u32 = 64;

/// Maximum Ziv-at-oracle working precision. On the measure-zero
/// inputs whose true value lies arbitrarily close to an f32 rounding
/// boundary, the loop caps here and the verdict is
/// [`Verdict::OracleInconclusive`] rather than `Ok` or `Mismatch`.
/// The honest caveat: a future sweep can raise this cap for a
/// specific function if the inconclusive set warrants it.
pub const MAX_PREC: u32 = 1024;

/// Closure that evaluates pfloat's kernel for a given `FnId` at
/// `f32` precision under the requested mode, returning the result's
/// f32 bit pattern. The verifier holds a reference; the per-function
/// driver builds the closure that dispatches to the right pfloat
/// API.
pub type Kernel<'a> = dyn Fn(FnId, u32, RoundingMode) -> u32 + 'a;

/// Returns the unique `f32` both endpoints of `enc` round to under
/// `mode`, or `None` when they straddle a rounding boundary (or
/// when either endpoint is NaN, treated as inconclusive for the
/// f32-rounding check).
pub fn certified_round_f32(enc: &Enclosure, mode: RoundingMode) -> Option<f32> {
    let lo_r = round_f32(&enc.lo, mode)?;
    let hi_r = round_f32(&enc.hi, mode)?;
    if lo_r.to_bits() == hi_r.to_bits() {
        Some(lo_r)
    } else {
        None
    }
}

/// Verify one `(function, input, mode)` triple. Runs the Ziv-at-
/// oracle loop until the enclosure certifies a unique f32 or the
/// precision cap is reached. Compares the certified `f32` to
/// pfloat's `f32` result from `kernel`. Returns the verdict.
///
/// `kernel` is the per-function pfloat dispatch the driver provides.
/// `f32` inputs that map to NaN under pfloat get reported via
/// `Verdict::Mismatch` (the `got` field carries pfloat's NaN bit
/// pattern); NaN-vs-NaN bit equality is meaningful here even though
/// IEEE says `NaN != NaN` for arithmetic comparison.
pub fn verify_input(
    oracle: &dyn OracleBackend,
    f: FnId,
    input: u32,
    mode: RoundingMode,
    kernel: &Kernel<'_>,
) -> Verdict {
    let mut prec = START_PREC;
    loop {
        let enc = oracle.enclose(f, input, prec);
        if let Some(expected) = certified_round_f32(&enc, mode) {
            let got = kernel(f, input, mode);
            return if got == expected.to_bits() {
                Verdict::Ok
            } else {
                Verdict::Mismatch {
                    input,
                    mode,
                    expected: expected.to_bits(),
                    got,
                }
            };
        }
        if prec >= MAX_PREC {
            return Verdict::OracleInconclusive { input, mode };
        }
        prec = (prec * 2).min(MAX_PREC);
    }
}
