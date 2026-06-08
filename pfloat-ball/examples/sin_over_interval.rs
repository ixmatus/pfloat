//! A rigorous enclosure of `sin` over an input interval.
//!
//! Run with: `cargo run --example sin_over_interval --features trig`
//!
//! A `Ball` denotes a closed interval, and a ball operation returns a sound
//! superset of the image: `sin` applied to EVERY point of the input interval
//! is provably inside the result ball. This example encloses `sin` over a
//! tight interval around `1.0` and checks that the true `sin` at an interior
//! point lands inside the returned bounds. The radius, not a status flag, is
//! the accuracy channel for a ball (the spec's Law 5).

use pfloat::{BigFloat, RoundingMode};
use pfloat_ball::{Ball, Mag};

const NE: RoundingMode = RoundingMode::NearestEven;

fn main() {
    // A tight input interval: 1.0 plus or minus 2^-20, at 128-bit precision.
    let mid = BigFloat::try_from_i64_exact(1, 128).unwrap();
    let x = Ball::new(mid.clone(), Mag::from_pow2(-20)).unwrap();

    // Enclose sin over the whole interval.
    let (s, _status) = x.sin();
    let lower = s.lower();
    let upper = s.upper();
    println!("sin([1 - 2^-20, 1 + 2^-20]) is enclosed by");
    println!("  [{lower},");
    println!("   {upper}]");

    // The enclosure must contain sin(t) for every t in the interval. Check an
    // interior point that is not the midpoint: t = 1 + 2^-22.
    let (eps, _) = BigFloat::try_from_i64_exact(1, 128)
        .unwrap()
        .scale_by_pow2(-22);
    let (t, _) = mid.add(&eps, NE);
    let (sin_t, _) = t.sin(NE);
    let at_or_above_lower = lower.partial_cmp(&sin_t).0 != Some(core::cmp::Ordering::Greater);
    let at_or_below_upper = upper.partial_cmp(&sin_t).0 != Some(core::cmp::Ordering::Less);
    assert!(
        at_or_above_lower && at_or_below_upper,
        "sin(1 + 2^-22) must lie inside the enclosure"
    );
    println!("sin(1 + 2^-22) is provably inside the enclosure");

    // The radius surfaces the accuracy directly: how many leading bits of the
    // result the enclosure pins down.
    println!("certified accuracy: {} bits", s.rel_accuracy_bits());
}
