//! serde round-trip and deserialize-validation tests for the `serde`
//! feature (ADR-0068).
//!
//! The impls do not branch on `is_human_readable`, so a JSON round-trip
//! exercises them fully; the rejection tests confirm the deserialize
//! trust boundary rejects malformed input rather than coercing it.

// `FixedFloat<PREC>`'s `[(); limbs_for(PREC)]` bound needs the same
// nightly feature the library declares.
#![feature(generic_const_exprs)]
#![allow(incomplete_features)]
#![cfg(all(feature = "serde", feature = "fixed", feature = "std"))]

use serde::de::DeserializeOwned;
use serde::Serialize;

use pfloat::{BigFloat, BuildError, FixedFloat, IeeeClass, ParseError, RoundingMode, Sign, Status};

/// Serialize then deserialize through JSON; the result must equal the
/// input bit-for-bit (`BigFloat`'s `PartialEq` is structural).
fn rt<T>(x: &T)
where
    T: Serialize + DeserializeOwned + PartialEq + core::fmt::Debug,
{
    let s = serde_json::to_string(x).expect("serialize");
    let y: T = serde_json::from_str(&s).expect("deserialize");
    assert_eq!(*x, y, "round-trip mismatch via {s}");
}

fn big(n: i64, p: u32) -> BigFloat {
    BigFloat::try_from_i64_exact(n, p).expect("precision >= 1")
}

#[test]
fn bigfloat_specials_round_trip() {
    rt(&BigFloat::try_new_zero(Sign::Positive, 53).unwrap());
    rt(&BigFloat::try_new_zero(Sign::Negative, 53).unwrap());
    rt(&BigFloat::try_new_infinity(Sign::Positive, 113).unwrap());
    rt(&BigFloat::try_new_infinity(Sign::Negative, 113).unwrap());
    // Quiet and signaling NaNs, with and without a nonzero payload.
    rt(&BigFloat::try_new_quiet_nan(Sign::Positive, 53, &[]).unwrap());
    rt(&BigFloat::try_new_quiet_nan(Sign::Negative, 53, &[0xABCD]).unwrap());
    rt(&BigFloat::try_new_signaling_nan(Sign::Positive, 128, &[0x1, 0x2]).unwrap());
}

#[test]
fn bigfloat_normals_round_trip() {
    // Integers and a non-dyadic quotient across precisions that span the
    // storage-padding edge cases: pad>0 (53, 65, 113), pad==0 (64, 128).
    for &(n, p) in &[
        (42i64, 53u32),
        (-7, 113),
        (3, 64),
        (3, 65),
        (5, 128),
        (-1, 256),
    ] {
        rt(&big(n, p));
    }
    let three = big(3, 256);
    let one = big(1, 256);
    let (third, _) = one.div(&three, RoundingMode::NearestEven);
    rt(&third);
}

#[test]
fn fixedfloat_round_trips() {
    let v = FixedFloat::<53>::try_from_big_exact(big(42, 53)).unwrap();
    rt(&v);
    let w = FixedFloat::<128>::try_from_big_exact(big(-9, 128)).unwrap();
    rt(&w);
}

#[test]
fn simple_types_round_trip() {
    rt(&Sign::Positive);
    rt(&Sign::Negative);
    for &m in &[
        RoundingMode::NearestEven,
        RoundingMode::NearestAway,
        RoundingMode::TowardZero,
        RoundingMode::TowardPositive,
        RoundingMode::TowardNegative,
    ] {
        rt(&m);
    }
    rt(&Status::OK);
    rt(&Status::INVALID);
    rt(&Status::INEXACT);
    rt(&IeeeClass::PositiveNormal);
    rt(&IeeeClass::QuietNaN);
    rt(&BuildError::PrecisionZero);
    rt(&BuildError::ValueExceedsPrecision {
        value_bits: 60,
        requested: 53,
    });
    rt(&ParseError::Empty);
    rt(&ParseError::InvalidExponent);
}

// ----- deserialize trust boundary: malformed input is rejected -----

fn rejected(json: &str) {
    let r: Result<BigFloat, _> = serde_json::from_str(json);
    assert!(r.is_err(), "expected rejection, got {r:?} for {json}");
}

#[test]
fn deserialize_rejects_zero_precision() {
    rejected(r#"{"precision":0,"class":{"Zero":{"sign":"Positive"}}}"#);
}

#[test]
fn deserialize_rejects_wrong_mantissa_limb_count() {
    // precision 53 needs exactly one limb; two is malformed.
    rejected(
        r#"{"precision":53,"class":{"Normal":{"sign":"Positive","exponent":0,"mantissa":[9223372036854775808,0]}}}"#,
    );
}

#[test]
fn deserialize_rejects_unnormalized_mantissa() {
    // Top bit of the most-significant limb clear.
    rejected(
        r#"{"precision":53,"class":{"Normal":{"sign":"Positive","exponent":0,"mantissa":[1]}}}"#,
    );
}

#[test]
fn deserialize_rejects_padding_bits_set() {
    // precision 53 leaves 11 padding bits; bit 0 is inside them.
    rejected(
        r#"{"precision":53,"class":{"Normal":{"sign":"Positive","exponent":0,"mantissa":[9223372036854775809]}}}"#,
    );
}

#[test]
fn deserialize_rejects_wrong_nan_payload_length() {
    // precision 53 needs a one-limb payload; empty is malformed.
    rejected(r#"{"precision":53,"class":{"Nan":{"quiet":true,"sign":"Positive","payload":[]}}}"#);
}

#[test]
fn deserialize_accepts_canonical_normal() {
    // 1.0 at precision 53: top bit set, exponent 0, padding clear.
    let json = r#"{"precision":53,"class":{"Normal":{"sign":"Positive","exponent":0,"mantissa":[9223372036854775808]}}}"#;
    let v: BigFloat = serde_json::from_str(json).expect("canonical value accepted");
    assert_eq!(v, big(1, 53));
}

#[test]
fn fixedfloat_rejects_precision_mismatch() {
    // A value serialized at precision 64 must not deserialize into
    // FixedFloat<53>.
    let s = serde_json::to_string(&big(3, 64)).unwrap();
    let r: Result<FixedFloat<53>, _> = serde_json::from_str(&s);
    assert!(r.is_err(), "precision 64 must not load into FixedFloat<53>");
}
