//! Smoke tests for the Oracle harness type definitions.
//!
//! Confirms `Enclosure`, `OracleBackend`, `FnId`, and `Verdict` from
//! `tests/oracle/types.rs` compile in a test crate and carry the
//! shape the verifier and runners depend on.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod oracle;

use oracle::{Enclosure, FnId, OracleBackend, Verdict};
use pfloat::RoundingMode;
use rug::Float;

#[test]
fn enclosure_carries_two_endpoints_at_working_precision() {
    let p = 64;
    let lo = Float::with_val(p, 1.0);
    let hi = Float::with_val(p, 2.0);
    let enc = Enclosure {
        lo: lo.clone(),
        hi: hi.clone(),
    };
    assert!(enc.lo <= enc.hi);
    assert_eq!(enc.lo.prec(), p);
    assert_eq!(enc.hi.prec(), p);
}

#[test]
fn fnid_names_match_v1_0_surface_doc() {
    // A spot-check across each variant kind; the full list is
    // structural and the docs/v1.0-surface.md document is the
    // authoritative inventory.
    assert_eq!(FnId::Sqrt.name(), "sqrt");
    assert_eq!(FnId::Exp.name(), "exp");
    assert_eq!(FnId::Ln.name(), "ln");
    assert_eq!(FnId::Sin.name(), "sin");
    assert_eq!(FnId::Gamma.name(), "gamma");
    assert_eq!(FnId::Lgamma.name(), "lgamma");
    assert_eq!(FnId::Ai.name(), "Ai");
    assert_eq!(FnId::AiPrime.name(), "Ai_prime");
    assert_eq!(FnId::BesselJ0.name(), "J0");
    assert_eq!(FnId::BesselJn(5).name(), "Jn");
    assert_eq!(FnId::BesselI0.name(), "I0");
    assert_eq!(FnId::BesselKn(-2).name(), "Kn");
    assert_eq!(FnId::Ei.name(), "Ei");
    assert_eq!(FnId::Li.name(), "li");
}

#[test]
fn fnid_parametric_orders_distinguish_in_eq() {
    assert_ne!(FnId::BesselJn(3), FnId::BesselJn(5));
    assert_eq!(FnId::BesselJn(3), FnId::BesselJn(3));
    assert_ne!(FnId::BesselJ0, FnId::BesselJn(0));
}

#[test]
fn verdict_variants_carry_diagnostic_payloads() {
    let ok = Verdict::Ok;
    let mismatch = Verdict::Mismatch {
        input: 0x3f800000,
        mode: RoundingMode::NearestEven,
        expected: 0x3f800000,
        got: 0x3f800001,
    };
    let inconclusive = Verdict::OracleInconclusive {
        input: 0x7f7fffff,
        mode: RoundingMode::NearestEven,
    };
    let panic = Verdict::Panic {
        input: 0xffc00000,
        mode: RoundingMode::NearestEven,
        message: "unreachable kernel branch".into(),
    };
    // Debug renders without panicking.
    let _ = format!("{ok:?} {mismatch:?} {inconclusive:?} {panic:?}");
}

/// Minimal `OracleBackend` impl so the trait shape is exercised by a
/// concrete type in tests; the real MPFR backend lands in commit 3.
struct StubBackend;

impl OracleBackend for StubBackend {
    fn enclose(&self, _f: FnId, _input: u32, working_prec: u32) -> Enclosure {
        let lo = Float::with_val(working_prec, 0.0);
        let hi = Float::with_val(working_prec, 0.0);
        Enclosure { lo, hi }
    }

    fn name(&self) -> &'static str {
        "stub"
    }
}

#[test]
fn oracle_backend_trait_object_dispatches_through_dyn() {
    let stub: &dyn OracleBackend = &StubBackend;
    let enc = stub.enclose(FnId::Exp, 0x3f800000, 64);
    assert_eq!(enc.lo.prec(), 64);
    assert_eq!(stub.name(), "stub");
}
