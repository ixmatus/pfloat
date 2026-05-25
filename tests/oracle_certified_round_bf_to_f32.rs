//! Exhaustive property test for
//! [`oracle::convert::certified_round_bf_to_f32`] (Phase 1f slice
//! p1.23, ADR-0038).
//!
//! The mode-aware `certified_round_bf_to_f32` is the load-bearing
//! piece that lifts the silent NE-only constraint of the existing
//! Display + Rust `f32` parser bridge (`bf_to_f32_bits`, per
//! `feedback_bf_to_f32_directed_mode`). The two bridges must agree
//! bit-exact under `NearestEven` on every binary32 input — at
//! `p = 24` the BigFloat lands on the f32 grid and both routes are
//! exact re-encodes, so any disagreement signals a bridge bug.
//!
//! This is the 2^32-input property test. Wall-clock budget ~2
//! minutes on a single thread (the rug FFI dominates); marked
//! `#[ignore]` so the per-push CI gate stays Python-free and
//! within its compute budget per ADR-0035. Run via:
//!
//! ```bash
//! cargo test --features differential-mpfr \
//!     --test oracle_certified_round_bf_to_f32 -- --ignored
//! ```
//!
//! The four directed modes are exercised via the pinned-corpus +
//! sweep coverage; full 2^32 directed-mode sweeps would mean four
//! more 2-minute passes here and are deferred to the per-family
//! slice gates where they fall out of the per-mode status TOML
//! migration's sweep verdict.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod oracle;

use oracle::{bf24_of_bits, bf_to_f32_bits, certified_round_bf_to_f32};
use pfloat::RoundingMode;

#[test]
#[ignore = "2^32 exhaustive sweep; ~2 minutes wall-clock"]
fn certified_round_bf_to_f32_agrees_with_bf_to_f32_bits_under_ne_at_p24() {
    let mut mismatch_count = 0u64;
    let mut first_mismatch: Option<(u32, u32, u32)> = None;
    for bits in 0u32..=u32::MAX {
        let bf = bf24_of_bits(bits);
        let existing = bf_to_f32_bits(&bf);
        let Some(certified) = certified_round_bf_to_f32(&bf, RoundingMode::NearestEven) else {
            // NaN: existing bridge produces a NaN bit pattern via
            // parse; certified returns None per the IEEE NaN
            // convention. Skip the bit-pattern comparison for NaN
            // inputs (no unique f32 they round to).
            assert_eq!(
                (existing >> 23) & 0xFF,
                0xFF,
                "certified None on non-NaN existing pattern 0x{existing:08x}"
            );
            assert_ne!(
                existing & 0x007F_FFFF,
                0,
                "certified None on infinity (mant=0) existing pattern 0x{existing:08x}"
            );
            continue;
        };
        if certified != existing {
            mismatch_count += 1;
            if first_mismatch.is_none() {
                first_mismatch = Some((bits, existing, certified));
            }
        }
    }
    assert_eq!(
        mismatch_count, 0,
        "{mismatch_count} mismatches over 2^32 inputs; first: {first_mismatch:?}"
    );
}

#[test]
fn certified_round_bf_to_f32_handles_signaling_inputs_at_p24() {
    // A small sanity gate that runs in the default test suite (not
    // ignored) so the bridge gets exercised on every cargo test
    // pass: cover the six binary32 representative classes — zero,
    // subnormal, smallest normal, midrange normal, largest normal,
    // infinity. NaN tested separately because the f32 parser may
    // emit any quiet-NaN payload, so the bit pattern comparison
    // would compare implementation choices, not correctness.
    let representative_bits: &[u32] = &[
        0x0000_0000, // +0
        0x0000_0001, // smallest +subnormal
        0x007F_FFFF, // largest +subnormal
        0x0080_0000, // smallest +normal
        0x3F80_0000, // +1.0
        0x7F7F_FFFF, // largest +normal
        0x7F80_0000, // +inf
        0x8000_0000, // -0
    ];
    for &bits in representative_bits {
        let bf = bf24_of_bits(bits);
        let existing = bf_to_f32_bits(&bf);
        let certified = certified_round_bf_to_f32(&bf, RoundingMode::NearestEven)
            .expect("non-NaN input: certified is Some");
        assert_eq!(
            certified, existing,
            "p=24 NE bridge disagreement at 0x{bits:08x}: existing=0x{existing:08x}, \
             certified=0x{certified:08x}"
        );
    }
}
