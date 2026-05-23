//! Smoke tests for the verification core.
//!
//! Covers `certified_round_f32` (oracle bracket → unique f32 under
//! a rounding mode, or None), the f32 conversion bridges
//! (`bf24_of_bits`, `bf_to_f32_bits`, `round_f32`), and the
//! `verify_input` Ziv-at-oracle loop against the live MPFR backend.
//! Pfloat's `sqrt` kernel serves as the kernel-under-test for the
//! end-to-end check: sqrt is the only function the MPFR backend
//! has wired at slice p1.3.2; the broader-surface verification
//! lands at the smoke gate in slice p1.3.5 once the full pfloat
//! kernel dispatch (slice p1.3.4 in the original numbering, now
//! folded forward) ships.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod oracle;

use oracle::{
    bf24_of_bits, bf_to_f32_bits, certified_round_f32, round_f32, verify_input, Enclosure, FnId,
    Kernel, MpfrOracle, Verdict, MAX_PREC, START_PREC,
};
use pfloat::RoundingMode;
use rug::float::Round;
use rug::Float;

const NE: RoundingMode = RoundingMode::NearestEven;

// --- f32 ↔ BigFloat conversion bridges ---

#[test]
fn bf24_of_bits_round_trips_positive_normal() {
    let one_bits = 1.0_f32.to_bits();
    let bf = bf24_of_bits(one_bits);
    assert_eq!(bf.precision(), 24);
    assert_eq!(bf_to_f32_bits(&bf), one_bits);
}

#[test]
fn bf24_of_bits_round_trips_negative_normal() {
    let neg_pi_bits = (-core::f32::consts::PI).to_bits();
    let bf = bf24_of_bits(neg_pi_bits);
    assert_eq!(bf_to_f32_bits(&bf), neg_pi_bits);
}

#[test]
fn bf24_of_bits_round_trips_subnormal_min() {
    // f32 subnormal min = 0x00000001 = 2^-149.
    let sub_min = 1u32;
    let bf = bf24_of_bits(sub_min);
    assert_eq!(bf_to_f32_bits(&bf), sub_min);
}

#[test]
fn bf24_of_bits_round_trips_subnormal_max() {
    // f32 subnormal max = 0x007FFFFF.
    let sub_max = 0x007F_FFFF;
    let bf = bf24_of_bits(sub_max);
    assert_eq!(bf_to_f32_bits(&bf), sub_max);
}

#[test]
fn bf24_of_bits_handles_signed_zero() {
    assert_eq!(bf_to_f32_bits(&bf24_of_bits(0)), 0);
    assert_eq!(bf_to_f32_bits(&bf24_of_bits(0x8000_0000)), 0x8000_0000);
}

#[test]
fn bf24_of_bits_handles_infinity() {
    let pos_inf = f32::INFINITY.to_bits();
    let neg_inf = f32::NEG_INFINITY.to_bits();
    assert_eq!(bf_to_f32_bits(&bf24_of_bits(pos_inf)), pos_inf);
    assert_eq!(bf_to_f32_bits(&bf24_of_bits(neg_inf)), neg_inf);
}

// --- round_f32 (rug Float → f32 under mode) ---

#[test]
fn round_f32_returns_none_for_nan() {
    let nan = Float::with_val(64, rug::float::Special::Nan);
    assert!(round_f32(&nan, NE).is_none());
}

#[test]
fn round_f32_nearest_even_matches_rug_to_nearest() {
    // A non-trivial value at p=64: sqrt(2). `Float::with_val_round`
    // with the sqrt incomplete evaluates at p=64 under nearest.
    let two = Float::with_val(64, 2.0_f32);
    let (sqrt2, _) = Float::with_val_round(64, two.sqrt_ref(), Round::Nearest);
    let r_ne = round_f32(&sqrt2, NE).unwrap();
    let expected = sqrt2.to_f32_round(Round::Nearest);
    assert_eq!(r_ne.to_bits(), expected.to_bits());
}

#[test]
fn round_f32_directed_modes_match_rug() {
    let two = Float::with_val(64, 2.0_f32);
    let (sqrt2, _) = Float::with_val_round(64, two.sqrt_ref(), Round::Nearest);
    assert_eq!(
        round_f32(&sqrt2, RoundingMode::TowardZero)
            .unwrap()
            .to_bits(),
        sqrt2.to_f32_round(Round::Zero).to_bits()
    );
    assert_eq!(
        round_f32(&sqrt2, RoundingMode::TowardPositive)
            .unwrap()
            .to_bits(),
        sqrt2.to_f32_round(Round::Up).to_bits()
    );
    assert_eq!(
        round_f32(&sqrt2, RoundingMode::TowardNegative)
            .unwrap()
            .to_bits(),
        sqrt2.to_f32_round(Round::Down).to_bits()
    );
}

#[test]
fn round_f32_nearest_away_differs_from_nearest_even_on_ties() {
    // 1 ULP at f32 above 1.0 = 2^-23. Half-ULP tie = 1 + 2^-24.
    // NE rounds to even (1.0); NA rounds away (next f32).
    let prec = 64;
    let one = Float::with_val(prec, 1.0);
    let half_ulp = Float::with_val(prec, 1u32) >> 24u32;
    let tie = Float::with_val(prec, &one + &half_ulp);
    let ne = round_f32(&tie, NE).unwrap();
    let na = round_f32(&tie, RoundingMode::NearestAway).unwrap();
    assert_eq!(ne.to_bits(), 1.0_f32.to_bits()); // even
    assert!(na.to_bits() > ne.to_bits()); // away from zero
}

// --- certified_round_f32 ---

#[test]
fn certified_returns_some_when_both_endpoints_round_same() {
    let p = 64;
    // Both endpoints equal -> trivially certified.
    let lo = Float::with_val(p, 1.0_f32);
    let hi = Float::with_val(p, 1.0_f32);
    let enc = Enclosure { lo, hi };
    let r = certified_round_f32(&enc, NE).unwrap();
    assert_eq!(r.to_bits(), 1.0_f32.to_bits());
}

#[test]
fn certified_returns_none_when_endpoints_straddle_f32_boundary() {
    let p = 64;
    let one = Float::with_val(p, 1.0_f32);
    // A bracket spanning f32 ULPs at 1.0: lo just below 1.0, hi
    // just above 1.0 + 1 ULP. The NE-rounded values differ at f32.
    let ulp_f32 = Float::with_val(p, 1u32) >> 23u32; // 2^-23 (f32 ULP at 1.0)
    let lo = Float::with_val(p, &one - &ulp_f32);
    let hi = Float::with_val(p, &one + &ulp_f32);
    let enc = Enclosure { lo, hi };
    assert!(certified_round_f32(&enc, NE).is_none());
}

#[test]
fn certified_returns_none_when_endpoint_is_nan() {
    let p = 64;
    let nan = Float::with_val(p, rug::float::Special::Nan);
    let one = Float::with_val(p, 1.0_f32);
    let enc = Enclosure {
        lo: nan.clone(),
        hi: one.clone(),
    };
    assert!(certified_round_f32(&enc, NE).is_none());
}

// --- verify_input end-to-end on sqrt ---

fn sqrt_kernel(_f: FnId, input: u32, mode: RoundingMode) -> u32 {
    // pfloat's sqrt at p=24 (binary32). Builds the f32 input via
    // the bit-exact bridge so subnormal inputs are preserved.
    let x = bf24_of_bits(input);
    let (result, _) = x.sqrt(mode);
    bf_to_f32_bits(&result)
}

#[test]
fn verify_sqrt_at_simple_input_returns_ok() {
    let o = MpfrOracle;
    let k: &Kernel = &sqrt_kernel;
    let v = verify_input(&o, FnId::Sqrt, 4.0_f32.to_bits(), NE, k);
    assert!(matches!(v, Verdict::Ok), "got {v:?}");
}

#[test]
fn verify_sqrt_at_zero_returns_ok() {
    let o = MpfrOracle;
    let k: &Kernel = &sqrt_kernel;
    let v = verify_input(&o, FnId::Sqrt, 0.0_f32.to_bits(), NE, k);
    assert!(matches!(v, Verdict::Ok), "got {v:?}");
}

#[test]
fn verify_sqrt_at_infinity_returns_ok() {
    let o = MpfrOracle;
    let k: &Kernel = &sqrt_kernel;
    let v = verify_input(&o, FnId::Sqrt, f32::INFINITY.to_bits(), NE, k);
    assert!(matches!(v, Verdict::Ok), "got {v:?}");
}

/// A non-trivial input where the oracle's first guard suffices
/// (the value is not arbitrarily close to a rounding boundary).
#[test]
fn verify_sqrt_at_two_returns_ok_under_all_modes() {
    let o = MpfrOracle;
    let k: &Kernel = &sqrt_kernel;
    for mode in [
        RoundingMode::NearestEven,
        RoundingMode::NearestAway,
        RoundingMode::TowardZero,
        RoundingMode::TowardPositive,
        RoundingMode::TowardNegative,
    ] {
        let v = verify_input(&o, FnId::Sqrt, 2.0_f32.to_bits(), mode, k);
        assert!(matches!(v, Verdict::Ok), "mode={mode:?}: {v:?}");
    }
}

/// A small dense sweep over the first 1024 f32 normals: every input
/// either certifies as `Ok` or returns `OracleInconclusive` (NaN
/// from sqrt of negative would surface as `Mismatch` only if pfloat
/// and the oracle disagreed on the NaN-bit-pattern; here every input
/// is non-negative). No `Mismatch`, no panic.
#[test]
fn verify_sqrt_dense_sweep_first_1024_normals_returns_ok() {
    let o = MpfrOracle;
    let k: &Kernel = &sqrt_kernel;
    let mut ok = 0u32;
    let mut inconclusive = 0u32;
    // f32 bit patterns 0x3f800000..0x3f800400: 1024 inputs starting
    // at 1.0 going up by 1 ULP each.
    for bits in 0x3f80_0000u32..0x3f80_0400 {
        match verify_input(&o, FnId::Sqrt, bits, NE, k) {
            Verdict::Ok => ok += 1,
            Verdict::OracleInconclusive { .. } => inconclusive += 1,
            other => panic!("unexpected verdict at {bits:#010x}: {other:?}"),
        }
    }
    assert_eq!(ok + inconclusive, 1024);
    // sqrt is one of the easy ones; expect no inconclusives in
    // this range.
    assert_eq!(
        inconclusive, 0,
        "sqrt sweep had {inconclusive} inconclusives"
    );
}

// --- constants are reasonable ---

#[test]
fn precision_bounds_are_sane() {
    assert!(START_PREC >= 24); // comfortably above f32 mantissa width
    assert!(MAX_PREC >= START_PREC * 8); // at least three doublings of headroom
    assert!(MAX_PREC <= 8192); // bounded so a pathological input fails fast
}
