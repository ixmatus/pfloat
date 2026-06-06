//! Fuzz target: arbitrary bytes through serde JSON deserialization of
//! [`pfloat::BigFloat`].
//!
//! Deserialization is a trust boundary (ADR-0068): it revalidates the
//! canonical form (precision, mantissa limb count, top bit set, padding
//! bits clear, NaN payload length) and must never panic, only `Err` on
//! malformed input. Coverage-guided fuzzing learns the wire shape and
//! drives the validation in `big_from_repr`. On success the value must
//! be canonical, which is checked by re-serializing and re-deserializing
//! to an identical value: a non-canonical survivor would diverge here.
//!
//! Per ADR-0013 the corpus is not checked in; libFuzzer evolves its
//! own. Counterexamples get promoted to `tests/serde_roundtrip.rs`.

#![no_main]

use libfuzzer_sys::fuzz_target;

use pfloat::BigFloat;

fuzz_target!(|data: &[u8]| {
    let Ok(bf) = serde_json::from_slice::<BigFloat>(data) else {
        return;
    };
    // A deserialized value is always canonical: it re-serializes and
    // re-deserializes to itself.
    let bytes = serde_json::to_vec(&bf).expect("serialize a deserialized value");
    let bf2 = serde_json::from_slice::<BigFloat>(&bytes).expect("re-deserialize canonical output");
    assert_eq!(bf, bf2);
});
