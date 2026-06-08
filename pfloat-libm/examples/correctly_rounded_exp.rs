//! A correctly-rounded `f32` `exp`, contrasted with `std`.
//!
//! Run with: `cargo run --example correctly_rounded_exp`
//!
//! pfloat-libm computes `exp` at high precision through pfloat's kernel and
//! rounds the result to `f32` only once an enclosure proves the rounding,
//! so the last bit is always correct. The directed `_round` form lets you
//! pick the IEEE rounding mode and read the sticky status flags.

use pfloat_libm::{f32 as lm, RoundingMode};

fn main() {
    let x: f32 = 1.5;

    // The bare call rounds to nearest even, correctly rounded.
    let r = lm::exp(x);
    println!("correctly-rounded exp(1.5) = {r}");

    // Round the true exp(1.5) toward +inf and toward -inf and read the flags.
    let (up, status) = lm::exp_round(x, RoundingMode::TowardPositive);
    let (down, _) = lm::exp_round(x, RoundingMode::TowardNegative);

    // exp(1.5) is transcendental, so it lands strictly between two f32 grid
    // points: the upward and downward roundings differ by exactly one ULP,
    // and the result is inexact.
    assert!(up > down);
    assert_eq!(up.to_bits() - down.to_bits(), 1);
    assert!(status.inexact());
    println!("toward +inf = {up}, toward -inf = {down} (one ULP apart)");

    // Saturation fast-path: an argument far past the f32 overflow point
    // returns +inf with OVERFLOW and INEXACT without paying the cost of
    // argument reduction.
    let (big, st) = lm::exp_round(1000.0, RoundingMode::NearestEven);
    assert!(big.is_infinite() && st.overflow() && st.inexact());
    println!("exp(1000.0) saturates to {big} (overflow flag set)");

    // For contrast: std's f32::exp is fast but not guaranteed correct in the
    // last bit, so it may differ from the correctly-rounded value above.
    println!("std exp(1.5) for contrast = {}", x.exp());
}
