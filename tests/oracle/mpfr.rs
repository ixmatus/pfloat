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
//! See ADR-0034 for the design rationale and ADR-0014 for the
//! existing MPFR differential lane's gating posture this module
//! inherits.

#![cfg(all(unix, feature = "differential-mpfr"))]

use rug::float::Round;
use rug::Float;

use super::types::{Enclosure, FnId, OracleBackend};

/// In-process MPFR-backed oracle. A unit struct: `rug` carries no
/// per-backend state, so one `MpfrOracle` value serves every call.
pub struct MpfrOracle;

impl OracleBackend for MpfrOracle {
    fn enclose(&self, f: FnId, input: u32, working_prec: u32) -> Enclosure {
        // Convert the f32 bit pattern exactly into a high-precision
        // Float; an f32 is representable without rounding at any
        // precision >= 24, so working_prec >= 64 (Phase 1's
        // START_PREC) is comfortably exact.
        let x = Float::with_val(working_prec, f32::from_bits(input));
        match f {
            FnId::Sqrt => sqrt_bracket(&x, working_prec),
            // Remaining MPFR-primary dispatch entries land in
            // commit p1.3.3.
            _ => unimplemented!("FnId::{:?} dispatch lands in slice p1.3.3", f),
        }
    }

    fn name(&self) -> &'static str {
        "MPFR"
    }
}

/// `[sqrt(x) with RNDD, sqrt(x) with RNDU]` at `working_prec`.
///
/// The `*_ref` accessor returns an `Incomplete` placeholder so the
/// directed rounding gets applied at the `with_val_round` assignment
/// rather than at an intermediate to-nearest step.
fn sqrt_bracket(x: &Float, working_prec: u32) -> Enclosure {
    let (lo, _) = Float::with_val_round(working_prec, x.sqrt_ref(), Round::Down);
    let (hi, _) = Float::with_val_round(working_prec, x.sqrt_ref(), Round::Up);
    Enclosure { lo, hi }
}
