//! Smoke test for the ADR-0039 (pf-tqzz) MIDPOINT verb.
//!
//! Confirms the Arb worker's MIDPOINT response parses correctly on
//! the Rust side, and that the returned [`rug::Float`] midpoint is
//! within a tiny tolerance of the true mathematical value for an
//! Arb-primary kernel (Si) at a representative input.
//!
//! The full pf-tqzz cross-check sweep (every kernel × every input ×
//! every mode) is the per-release harness deferred to a follow-up
//! sub-slice; this smoke validates the wire format end-to-end so
//! the follow-up's only remaining work is per-kernel scaffolding.
//!
//! The test is silently no-op when the Arb venv is not configured
//! (the `pfloat-arb-oracle` venv at `${HOME}/.cache/...` is not
//! provisioned by the default CI build per ADR-0034's LGPL-isolation
//! posture; production developers and the per-release CI lane run
//! `scripts/setup_arb_oracle.sh` first).

#![cfg(all(unix, feature = "differential-arb"))]

mod oracle;

use oracle::arb::ArbOracle;
use oracle::FnId;

/// Si(1) ≈ 0.9460830703671830149... (NIST DLMF 6.7.1 / Abramowitz
/// & Stegun 5.1.4 / standard reference value).
const SI_AT_ONE_REFERENCE_F64: f64 = 0.946_083_070_367_183;

#[test]
fn arb_midpoint_si_at_one_returns_close_to_reference() {
    let oracle = match ArbOracle::new() {
        Ok(o) => o,
        Err(e) => {
            eprintln!(
                "skipping arb_midpoint smoke (Arb venv unavailable): {e}\n\
                 (run scripts/setup_arb_oracle.sh to provision)"
            );
            return;
        }
    };

    let one_bits = 1.0f32.to_bits();
    // Request a generous oracle precision (128 bits) so the midpoint
    // sits well within f64-comparison tolerance of the true value.
    let mid = oracle
        .midpoint(FnId::Si, one_bits, 128)
        .expect("Si(1) midpoint at p=128");

    // The MIDPOINT verb returns the mode-independent midpoint at the
    // requested precision; comparing the f64 projection against the
    // reference value with a generous 1e-13 tolerance leaves plenty
    // of room for the reference rounding and the Arb ball radius.
    let mid_f64 = mid.to_f64();
    let gap = (mid_f64 - SI_AT_ONE_REFERENCE_F64).abs();
    assert!(
        gap < 1e-13,
        "Si(1) midpoint {mid_f64} vs reference {SI_AT_ONE_REFERENCE_F64}: gap {gap}"
    );
}

#[test]
fn arb_midpoint_si_at_zero_returns_zero() {
    let oracle = match ArbOracle::new() {
        Ok(o) => o,
        Err(_) => return,
    };

    let zero_bits = 0.0f32.to_bits();
    // Si(0) = 0 exactly; the worker should return the zero-encoded
    // midpoint (OK + 0 0 wire form).
    let mid = oracle
        .midpoint(FnId::Si, zero_bits, 128)
        .expect("Si(0) midpoint at p=128");
    assert!(
        mid.is_zero(),
        "Si(0) midpoint must be exact zero, got {mid}"
    );
}

#[test]
fn arb_midpoint_bessel_k0_at_one_is_finite_positive() {
    let oracle = match ArbOracle::new() {
        Ok(o) => o,
        Err(_) => return,
    };

    // K_0(1) ≈ 0.4210244382... (NIST DLMF 10.32.9 standard reference).
    // The smoke confirms the BesselK0 wire path returns a finite
    // positive midpoint at an arbitrary input the worker dispatches
    // through `arb.bessel_k`.
    let one_bits = 1.0f32.to_bits();
    let mid = oracle
        .midpoint(FnId::BesselK0, one_bits, 128)
        .expect("K_0(1) midpoint at p=128");
    let mid_f64 = mid.to_f64();
    assert!(
        mid_f64 > 0.42 && mid_f64 < 0.43,
        "K_0(1) midpoint {mid_f64} not in expected (0.42, 0.43)"
    );
}
