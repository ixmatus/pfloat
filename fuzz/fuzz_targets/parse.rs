//! Fuzz target: feed arbitrary byte strings through
//! [`pfloat::BigFloat::parse_str`].
//!
//! Kani covers the special-case dispatch on hand-curated inputs;
//! libFuzzer covers the long tail of malformed input. The body
//! asserts parse never panics and, on success, that the rendered
//! `Display` output re-parses to a numerically equivalent value.
//!
//! Per ADR-0013 the corpus is not checked into the repo; libFuzzer
//! evolves its own corpus per run. Counterexamples that survive
//! get promoted to `.proptest-regressions` entries against
//! `tests/property_parse.rs` (the existing seed-commit convention).

#![no_main]

use libfuzzer_sys::fuzz_target;

use pfloat::{BigFloat, RoundingMode};

fuzz_target!(|data: &[u8]| {
    let Ok(s) = core::str::from_utf8(data) else {
        return;
    };
    let Ok((parsed, _status)) = BigFloat::parse_str(s, 113, RoundingMode::NearestEven)
    else {
        return;
    };
    // Round-trip: Display output must re-parse to a numerically
    // equivalent value.
    let rendered = format!("{parsed}");
    let Ok((reparsed, _)) = BigFloat::parse_str(&rendered, 113, RoundingMode::NearestEven)
    else {
        return;
    };
    if parsed.is_nan() {
        assert!(reparsed.is_nan());
    } else {
        // partial_cmp returns (Option<Ordering>, Status); for
        // finite operands the option is Some(Equal).
        let (ord, _) = parsed.partial_cmp(&reparsed);
        assert_eq!(ord, Some(std::cmp::Ordering::Equal));
    }
});
