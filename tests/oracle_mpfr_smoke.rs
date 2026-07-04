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

use oracle::{Enclosed, Enclosure, FnId, MpfrOracle, OracleBackend};
use rug::float::Round;
use rug::Float;

use pfloat::RoundingMode;

const NE: RoundingMode = RoundingMode::NearestEven;

/// Unwrap the bracket from an [`Enclosed`]. The MPFR backend never
/// abstains (it always produces a directed-rounding bracket), so
/// `Inconclusive` here is unreachable; the helper keeps the smoke
/// tests' destructures on `Enclosure` after the pf-41ou sum-type
/// split.
fn bracket(e: Enclosed) -> Enclosure {
    match e {
        Enclosed::Bracket(b) => b,
        Enclosed::Inconclusive => unreachable!("MPFR backend always brackets"),
    }
}

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
    let Enclosure { lo, hi } = bracket(o.enclose(FnId::Sqrt, two_bits, NE, 64));
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
    let Enclosure { lo, hi } = bracket(o.enclose(FnId::Sqrt, four_bits, NE, 64));
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
    } = bracket(o.enclose(FnId::Sqrt, two_bits, NE, 64));
    let Enclosure {
        lo: lo_128,
        hi: hi_128,
    } = bracket(o.enclose(FnId::Sqrt, two_bits, NE, 128));
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
    let Enclosure { lo, hi } = bracket(o.enclose(FnId::Sqrt, zero_bits, NE, 64));
    let zero = Float::with_val(64, 0.0);
    assert_eq!(lo, zero);
    assert_eq!(hi, zero);
}

#[test]
fn sqrt_enclosure_handles_infinity_input() {
    let o = MpfrOracle;
    let inf_bits = f32::INFINITY.to_bits();
    let Enclosure { lo, hi } = bracket(o.enclose(FnId::Sqrt, inf_bits, NE, 64));
    assert!(lo.is_infinite() && lo.is_sign_positive());
    assert!(hi.is_infinite() && hi.is_sign_positive());
}

/// Every MPFR-primary `FnId` produces a non-empty bracket
/// (`lo <= hi`) on a representative input. The directed-rounding
/// guarantee should hold uniformly across the surface.
#[test]
fn mpfr_primary_enclosures_are_well_formed() {
    let o = MpfrOracle;
    let p = 64;
    // (FnId, input as f32) pairs chosen to be in each function's
    // domain. Negative-input domain restrictions:
    //   sqrt, ln, log1p, log2, log10, gamma (poles), Ei (x > 0),
    //   acosh (x >= 1), atanh (|x| < 1), asin/acos (|x| <= 1).
    // 0.5 falls in every domain except acosh; we use 1.5 for that.
    let half = 0.5_f32.to_bits();
    let two = 2.0_f32.to_bits();
    let small = 0.25_f32.to_bits();
    let cases: Vec<(FnId, u32)> = vec![
        (FnId::Sqrt, half),
        (FnId::Exp, half),
        (FnId::Exp2, half),
        (FnId::Exp10, half),
        (FnId::Expm1, half),
        (FnId::Ln, half),
        (FnId::Log1p, half),
        (FnId::Log2, half),
        (FnId::Log10, half),
        (FnId::Sin, half),
        (FnId::Cos, half),
        (FnId::Tan, half),
        (FnId::Asin, half),
        (FnId::Acos, half),
        (FnId::Atan, half),
        (FnId::Sinh, half),
        (FnId::Cosh, half),
        (FnId::Tanh, half),
        (FnId::Asinh, half),
        (FnId::Acosh, two),
        (FnId::Atanh, half),
        (FnId::Erf, half),
        (FnId::Erfc, half),
        (FnId::Gamma, half),
        (FnId::Lgamma, half),
        (FnId::Digamma, half),
        (FnId::Zeta, two),
        (FnId::Ei, half),
        (FnId::Ai, half),
        (FnId::BesselJ0, half),
        (FnId::BesselJ1, half),
        (FnId::BesselJn(3), half),
        (FnId::BesselY0, half),
        (FnId::BesselY1, half),
        (FnId::BesselYn(3), half),
    ];
    for (id, x_bits) in cases {
        let Enclosure { lo, hi } = bracket(o.enclose(id, x_bits, NE, p));
        // NaN endpoints can occur only if MPFR's primitive returns
        // NaN, which is not expected for any of these in-domain
        // inputs.
        assert!(
            !lo.is_nan() && !hi.is_nan(),
            "NaN endpoint for {id:?} at f32 {x_bits:#010x}: lo={lo}, hi={hi}"
        );
        // Bracket well-formedness: lo <= hi at every working
        // precision; the directed-rounding contract guarantees this.
        assert!(
            lo <= hi,
            "lo > hi for {id:?} at f32 {x_bits:#010x}: lo={lo}, hi={hi}"
        );
    }
    let _ = small; // placeholder for future small-input cases
}

/// The Arb-primary `FnId` variants must `unimplemented!()` rather
/// than silently fall through; the MPFR backend should not pretend
/// to cover functions MPFR has no primitive for.
#[test]
#[should_panic(expected = "requires the Arb backend")]
fn arb_only_fnids_panic_under_mpfr_backend() {
    let o = MpfrOracle;
    let _ = o.enclose(FnId::BesselI0, 0.5_f32.to_bits(), NE, 64);
}
