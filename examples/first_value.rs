//! Tutorial 1: your first correctly-rounded value to N digits.
//!
//! Run with: `cargo run --example first_value`
//!
//! Companion to `docs/guides/01-first-value.md`. Every line here is the code
//! the guide walks through.

use pfloat::{BigFloat, RoundingMode};

const NE: RoundingMode = RoundingMode::NearestEven;

fn main() {
    // Precision is a bit count, chosen up front. 200 bits holds about 60
    // decimal digits of headroom.
    let two = BigFloat::try_from_i64_exact(2, 200).unwrap();

    // Every kernel returns the value paired with a Status of IEEE sticky flags.
    let (root, status) = two.sqrt(NE);

    // sqrt(2) is irrational, so the result is inexact at any finite precision.
    assert!(status.inexact());

    // Output digits are chosen separately from the working precision. Ask for
    // 50 significant digits, correctly rounded to nearest even.
    let fifty = root.to_decimal_string(50, NE);
    println!("sqrt(2) to 50 digits = {fifty}");

    // Or let the formatter pick the shortest string that round-trips back to
    // the same value at this precision.
    let shortest = root.to_shortest_decimal_string();
    println!("sqrt(2) shortest     = {shortest}");

    // The first 50 digits of sqrt(2) are known; pin them so the example is
    // self-checking.
    assert!(fifty.starts_with("1.4142135623730950488016887242096980785696718753769"));
}
