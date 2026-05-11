//! Fuzz target: BigFloat::Display round-trip on integer inputs.
//!
//! Given an i64, build a BigFloat, format it, re-parse, and check
//! numerical equality. Panic-freedom plus the strong round-trip
//! invariant. Complements the `parse` target which goes the other
//! direction (raw bytes → parse → reformat → re-parse).

#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;

use pfloat::{BigFloat, RoundingMode};

#[derive(Arbitrary, Debug)]
struct Input {
    value: i64,
    prec: u32,
}

fuzz_target!(|input: Input| {
    let prec = (input.prec % 256).max(1);
    let Ok(v) = BigFloat::try_from_i64_exact(input.value, prec) else {
        return;
    };
    let rendered = format!("{v}");
    let Ok((reparsed, _)) =
        BigFloat::parse_str(&rendered, prec, RoundingMode::NearestEven)
    else {
        // If parse rejects pfloat's own Display output, that is a
        // bug. Surface as a panic.
        panic!("parse rejected Display output: {rendered}");
    };
    let (cmp, _) = v.partial_cmp(&reparsed);
    assert_eq!(
        cmp,
        Some(core::cmp::Ordering::Equal),
        "round-trip mismatch for {v} → '{rendered}' → {reparsed}"
    );
});
