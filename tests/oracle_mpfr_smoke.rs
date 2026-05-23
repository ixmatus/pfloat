//! Smoke tests for the MPFR backend's enclosure shape.
//!
//! Confirms the directed-rounding bracket from `RNDD` + `RNDU`
//! actually straddles MPFR's RNDN evaluation for a handful of
//! `f32` inputs, and that the bracket tightens monotonically as
//! the working precision rises. These are the load-bearing
//! properties the verifier's `certified_round_f32` step depends
//! on; if either ever broke, the entire harness loses its
//! correctness argument.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod oracle;

use oracle::{Enclosure, FnId, MpfrOracle, OracleBackend};
use rug::float::Round;
use rug::Float;

#[test]
fn mpfr_backend_name_is_mpfr() {
    let o = MpfrOracle;
    assert_eq!(o.name(), "MPFR");
}

#[test]
fn sqrt_enclosure_straddles_nearest_value() {
    let o = MpfrOracle;
    // sqrt(2) at p=64 NRD/NRU brackets the irrational, the
    // canonical hard-to-round example.
    let two_bits = 2.0_f32.to_bits();
    let Enclosure { lo, hi } = o.enclose(FnId::Sqrt, two_bits, 64);
    let middle = {
        let x = Float::with_val(64, 2.0_f32);
        Float::with_val_round(64, x.sqrt_ref(), Round::Nearest).0
    };
    assert!(lo <= middle, "lo {lo} > middle {middle}");
    assert!(middle <= hi, "middle {middle} > hi {hi}");
}

#[test]
fn sqrt_enclosure_exact_value_collapses_to_zero_width() {
    let o = MpfrOracle;
    // sqrt(4) = 2 exactly at any precision >= 2 bits.
    let four_bits = 4.0_f32.to_bits();
    let Enclosure { lo, hi } = o.enclose(FnId::Sqrt, four_bits, 64);
    assert_eq!(lo, hi, "exact sqrt should produce zero-width bracket");
    let two = Float::with_val(64, 2.0);
    assert_eq!(lo, two);
}

#[test]
fn sqrt_enclosure_tightens_with_precision() {
    let o = MpfrOracle;
    let two_bits = 2.0_f32.to_bits();
    let Enclosure {
        lo: lo_64,
        hi: hi_64,
    } = o.enclose(FnId::Sqrt, two_bits, 64);
    let Enclosure {
        lo: lo_128,
        hi: hi_128,
    } = o.enclose(FnId::Sqrt, two_bits, 128);
    // The 128-bit bracket must lie inside the 64-bit one (or be
    // tighter on at least one side). The directed rounding modes
    // guarantee monotone refinement.
    assert!(
        lo_64 <= lo_128,
        "lo did not tighten: 64-bit {lo_64} > 128-bit {lo_128}"
    );
    assert!(
        hi_128 <= hi_64,
        "hi did not tighten: 128-bit {hi_128} > 64-bit {hi_64}"
    );
    // And the 128-bit width is strictly less than the 64-bit width
    // (sqrt(2) is irrational, so both precisions miss but the
    // higher one misses by less).
    let width_64 = Float::with_val(128, &hi_64 - &lo_64);
    let width_128 = Float::with_val(128, &hi_128 - &lo_128);
    assert!(
        width_128 < width_64,
        "bracket did not tighten: width64={width_64}, width128={width_128}"
    );
}

#[test]
fn sqrt_enclosure_handles_zero_input() {
    let o = MpfrOracle;
    let zero_bits = 0.0_f32.to_bits();
    let Enclosure { lo, hi } = o.enclose(FnId::Sqrt, zero_bits, 64);
    let zero = Float::with_val(64, 0.0);
    assert_eq!(lo, zero);
    assert_eq!(hi, zero);
}

#[test]
fn sqrt_enclosure_handles_infinity_input() {
    let o = MpfrOracle;
    let inf_bits = f32::INFINITY.to_bits();
    let Enclosure { lo, hi } = o.enclose(FnId::Sqrt, inf_bits, 64);
    assert!(lo.is_infinite() && lo.is_sign_positive());
    assert!(hi.is_infinite() && hi.is_sign_positive());
}
