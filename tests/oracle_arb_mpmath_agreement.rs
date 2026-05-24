//! ADR-0035 Tier 2 cross-check: for the smoke-gate anchor inputs of
//! every non-parametric Arb-primary `FnId`, the Arb worker and the
//! mpmath worker MUST agree on the certified `f32` bit pattern
//! under NE.
//!
//! Disagreement between two independent multi-precision libraries
//! (FLINT/Arb and mpmath share no code lineage) is vanishingly
//! unlikely to be coincident, so a divergence here is a strong
//! signal that one of the two has a silent defect that the
//! single-oracle protocol (slice p1.7's pf-6a4e episode showed
//! was possible) would have hidden.
//!
//! Cost: ~3 seconds debug for 10 functions x 16 inputs x 2 oracles.
//! Not in the per-push CI lane (still gated on `differential-arb`
//! which requires the Python venv); runs at slice-close cadence
//! alongside the full sweep.

#![cfg(all(unix, feature = "differential-arb"))]

#[path = "oracle/mod.rs"]
mod oracle;

use oracle::{ArbOracle, Enclosure, FnId, MpmathOracle, OracleBackend};
use pfloat::RoundingMode;
use rug::float::Round;

const NE: RoundingMode = RoundingMode::NearestEven;

/// 16 consecutive f32 inputs anchored at each function's natural
/// domain. Mirrors `tests/oracle_arb_smoke.rs::smoke_inputs`; if
/// that anchor set changes, this should too.
fn smoke_inputs(f: FnId) -> std::ops::Range<u32> {
    let anchor = match f {
        FnId::BesselI0
        | FnId::BesselI1
        | FnId::BesselIn(_)
        | FnId::BesselK0
        | FnId::BesselK1
        | FnId::BesselKn(_)
        | FnId::Bi
        | FnId::AiPrime
        | FnId::BiPrime => 0x3f80_0000u32, // 1.0
        FnId::Si | FnId::Ci => 0x3f00_0000u32, // 0.5
        FnId::Li => 0x4000_0000u32,            // 2.0
        _ => panic!("smoke_inputs called with non-Arb-primary FnId: {f:?}"),
    };
    anchor..(anchor + 16)
}

const ARB_PRIMARY_FNIDS: &[FnId] = &[
    FnId::Si,
    FnId::Ci,
    FnId::Li,
    FnId::Bi,
    FnId::AiPrime,
    FnId::BiPrime,
    FnId::BesselI0,
    FnId::BesselI1,
    FnId::BesselK0,
    FnId::BesselK1,
];

/// Extract the f32 bit pattern that an authoritative backend's
/// single-point [`Enclosure`] certifies. Authoritative backends
/// (per ADR-0035) return an enclosure with both endpoints equal to
/// the certified `f32`; the bit pattern is recovered via
/// `to_f32_round(Round::Nearest)` on either endpoint (the value is
/// exactly representable at f32 precision so the rounding is
/// trivial). NaN endpoints decode as a sentinel `f32::NAN` bit
/// pattern (worker returned `INC`); the cross-check treats two NaN
/// sentinels as agreeing.
fn extract_certified_f32(enc: &Enclosure) -> u32 {
    if enc.lo.is_nan() && enc.hi.is_nan() {
        return f32::NAN.to_bits();
    }
    let lo_f32 = enc.lo.to_f32_round(Round::Nearest);
    let hi_f32 = enc.hi.to_f32_round(Round::Nearest);
    // Single-point enclosures from the authoritative worker have
    // lo == hi; the rounding gives the same f32. We assert this
    // invariant to catch any worker that emits a multi-point
    // enclosure unexpectedly.
    assert_eq!(
        lo_f32.to_bits(),
        hi_f32.to_bits(),
        "authoritative enclosure endpoints disagree on f32: lo={lo_f32}, hi={hi_f32}"
    );
    lo_f32.to_bits()
}

#[test]
fn arb_and_mpmath_workers_agree_on_smoke_inputs() {
    let arb = ArbOracle::new()
        .expect("ArbOracle::new (requires the python-flint venv; run scripts/setup_arb_oracle.sh)");
    let mpm = MpmathOracle::new()
        .expect("MpmathOracle::new (requires mpmath in the same venv; install via pip)");

    let mut divergences: Vec<(FnId, u32, u32, u32)> = Vec::new();
    let mut total_checks = 0u32;

    for &f in ARB_PRIMARY_FNIDS {
        for input in smoke_inputs(f) {
            let enc_arb = arb.enclose(f, input, NE, 64);
            let enc_mpm = mpm.enclose(f, input, NE, 64);
            let bits_arb = extract_certified_f32(&enc_arb);
            let bits_mpm = extract_certified_f32(&enc_mpm);
            total_checks += 1;
            if bits_arb != bits_mpm {
                divergences.push((f, input, bits_arb, bits_mpm));
                if divergences.len() <= 8 {
                    eprintln!(
                        "[arb-mpmath] DIVERGENCE: {f:?} input={input:#010x} \
                         arb={bits_arb:#010x} mpmath={bits_mpm:#010x}"
                    );
                }
            }
        }
    }

    eprintln!(
        "[arb-mpmath] checked {total_checks} (input, fn) pairs across {} functions; {} divergences",
        ARB_PRIMARY_FNIDS.len(),
        divergences.len()
    );

    assert!(
        divergences.is_empty(),
        "Arb and mpmath disagreed on {} of {} smoke inputs; see eprintln above",
        divergences.len(),
        total_checks
    );
}
