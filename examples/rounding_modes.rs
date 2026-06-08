//! Tutorial 3: choosing a rounding mode, and why directed modes matter.
//!
//! Run with: `cargo run --example rounding_modes`
//!
//! Companion to `docs/guides/03-rounding-modes.md`. Every line here is the code
//! the guide walks through. Default features (std + fmt + big) are enough; the
//! irrational value is sqrt(2), so no exp-log or trig feature is needed.

use pfloat::{BigFloat, RoundingMode};

fn main() {
    // Compute sqrt(2) once, with headroom, under the default mode. The working
    // precision (256 bits) is a separate choice from the rounding mode we apply
    // when we read the value back out to a concrete grid.
    let two = BigFloat::try_from_i64_exact(2, 256).unwrap();
    let (root, status) = two.sqrt(RoundingMode::NearestEven);

    // sqrt(2) is irrational, so it never lands exactly on a finite grid: the
    // mode genuinely changes the answer.
    assert!(status.inexact());

    // Apply each of the five modes by rounding the high precision value onto the
    // f64 grid. f64 is a fixed grid, so "the next representable value" is well
    // defined and we can compare bit patterns directly.
    let (ne, _) = root.to_f64_round(RoundingMode::NearestEven);
    let (na, _) = root.to_f64_round(RoundingMode::NearestAway);
    let (tz, _) = root.to_f64_round(RoundingMode::TowardZero);
    let (tp, _) = root.to_f64_round(RoundingMode::TowardPositive);
    let (tn, _) = root.to_f64_round(RoundingMode::TowardNegative);

    println!("NearestEven    = {ne:.17}  bits={:#018x}", ne.to_bits());
    println!("NearestAway    = {na:.17}  bits={:#018x}", na.to_bits());
    println!("TowardZero     = {tz:.17}  bits={:#018x}", tz.to_bits());
    println!("TowardPositive = {tp:.17}  bits={:#018x}", tp.to_bits());
    println!("TowardNegative = {tn:.17}  bits={:#018x}", tn.to_bits());

    // The directed pair brackets the true value: rounding down never exceeds
    // rounding to nearest, which never exceeds rounding up. The true sqrt(2)
    // lies in the closed interval [tn, tp].
    assert!(tn <= ne);
    assert!(ne <= tp);

    // The bracket is as tight as the grid allows: TowardNegative and
    // TowardPositive are adjacent f64 values, exactly one ULP apart. Their bit
    // patterns differ by one because consecutive finite f64 values of the same
    // sign have consecutive bit patterns.
    assert_eq!(tp.to_bits() - tn.to_bits(), 1);

    // sqrt(2) is positive, so "toward zero" rounds in the same direction as
    // "toward negative infinity": both step down to the smaller neighbor.
    assert_eq!(tz.to_bits(), tn.to_bits());

    // For this value the nearest neighbor happens to be the upper one, so both
    // nearest variants agree with TowardPositive. NearestEven is the IEEE 754
    // default; you reach for a directed mode when you need a guaranteed side.
    assert_eq!(ne.to_bits(), tp.to_bits());
    assert_eq!(na.to_bits(), ne.to_bits());

    // The same bracket is visible in decimal once you ask for enough digits to
    // pass the point where the two directed results diverge.
    let down = root.to_decimal_string(40, RoundingMode::TowardNegative);
    let up = root.to_decimal_string(40, RoundingMode::TowardPositive);
    println!("down (40 sig) = {down}");
    println!("up   (40 sig) = {up}");

    // Both decimal strings begin with the digits everyone agrees on; they part
    // ways only deep in the tail, which is the whole point of a tight bracket.
    assert!(down.starts_with("1.41421356237309504880168872420"));
    assert!(up.starts_with("1.41421356237309504880168872420"));
    assert_ne!(down, up);
}
