//! Verification core: the Ziv-at-oracle precision loop, width-generic.
//!
//! For one `(function, input, mode)` triple the verifier rounds the
//! oracle's enclosure endpoints to the hardware width. When both
//! endpoints round to the same value, every point in the bracket
//! (including the true value) rounds there too, so the correctly
//! rounded result is determined; compare it bit-for-bit to the shell's
//! output. When the endpoints straddle a rounding boundary, double the
//! oracle's working precision and tighten the bracket. The loop caps at
//! `MAX_PREC`; inputs that exhaust it are [`Verdict::OracleInconclusive`]
//! (a measure-zero hard-to-round input, recorded not failed).
//!
//! Mirrors `pfloat/tests/oracle/verify.rs`, but generic over [`Hw`] so
//! one loop serves f32 and f64, and it sees the shell's `Status` so the
//! flag gate can run after the value gate.

#![cfg(all(unix, feature = "differential-mpfr"))]

use pfloat_libm::RoundingMode;

use super::hw::Hw;
use super::oracle::enclose;
use super::status_gate::{check_flags, StatusGate};
use super::types::{Enclosure, LibmArg, LibmFnId, Verdict};

/// First Ziv-at-oracle working precision (above both f32's 24 and f64's
/// 53 mantissa bits).
pub const START_PREC: u32 = 64;

/// Precision cap. Hard-to-round inputs that still straddle here yield
/// [`Verdict::OracleInconclusive`].
pub const MAX_PREC: u32 = 1024;

/// The unique hardware float both endpoints of `enc` round to under
/// `mode`, or `None` when they straddle. A both-NaN bracket (the true
/// value is undefined) certifies a NaN result.
pub fn certified_round<H: Hw>(enc: &Enclosure, mode: RoundingMode) -> Option<H> {
    if enc.lo.is_nan() && enc.hi.is_nan() {
        return Some(H::nan());
    }
    let lo = H::round(&enc.lo, mode)?;
    let hi = H::round(&enc.hi, mode)?;
    if H::to_bits(lo) == H::to_bits(hi) {
        Some(lo)
    } else {
        None
    }
}

/// Verify one `(f, input, arg, mode)` triple under the given gate.
pub fn verify_input<H: Hw>(
    f: LibmFnId,
    input: H::Bits,
    arg: LibmArg,
    mode: RoundingMode,
    gate: StatusGate,
) -> Verdict {
    let mut prec = START_PREC;
    loop {
        let xf = H::lift(input, prec);
        let partner = match arg {
            LibmArg::HypotY(yb) => Some(H::lift(H::partner_bits(yb), prec)),
            LibmArg::None => None,
        };
        let enc = enclose(f, &xf, partner.as_ref(), prec);

        if let Some(expected) = certified_round::<H>(&enc, mode) {
            let (got_bits, got_status) = H::shell(f, input, arg, mode);
            let got = H::from_bits(got_bits);

            // VALUE gate: NaN-aware (any NaN matches any NaN), ±0
            // distinguished by bit pattern.
            let nan_match = H::is_nan(expected) && H::is_nan(got);
            if !(nan_match || H::to_bits(expected) == got_bits) {
                return Verdict::ValueMismatch {
                    input: H::bits_to_u64(input),
                    mode,
                    expected: H::bits_to_u64(H::to_bits(expected)),
                    got: H::bits_to_u64(got_bits),
                };
            }

            // FLAG gate: INVALID / DIV_BY_ZERO (run only after value
            // matched).
            let xv = H::from_bits(input);
            if let Some((flag, exp, g)) =
                check_flags(gate, f, H::is_nan(xv), H::as_f64(xv), &enc, got_status)
            {
                return Verdict::FlagMismatch {
                    input: H::bits_to_u64(input),
                    mode,
                    flag,
                    expected: exp,
                    got: g,
                };
            }
            return Verdict::Ok;
        }

        if prec >= MAX_PREC {
            return Verdict::OracleInconclusive {
                input: H::bits_to_u64(input),
                mode,
            };
        }
        prec = (prec * 2).min(MAX_PREC);
    }
}
