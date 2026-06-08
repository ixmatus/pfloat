//! Using pfloat-ball as a rigorous self-oracle for a pfloat scalar.
//!
//! Run with: `cargo run --example ball_oracle`
//!
//! The two crates share one scalar engine, so the same value computed as a
//! `pfloat` scalar and as a `pfloat-ball` enclosure must agree: a correctly
//! rounded scalar always lies between the ball's lower and upper endpoints.
//! This is the cross-crate composition pattern, and it needs no extra feature
//! (`Ball::sqrt` is part of the always-available arithmetic surface).

use pfloat::{BigFloat, RoundingMode};
use pfloat_ball::Ball;

fn main() {
    // A pfloat scalar: sqrt(2), correctly rounded to nearest even at 200 bits.
    let two = BigFloat::try_from_i64_exact(2, 200).unwrap();
    let (scalar_sqrt2, _status) = two.sqrt(RoundingMode::NearestEven);

    // The same computation as a rigorous enclosure of an exact (point) ball.
    let ball_two = Ball::point(two).unwrap();
    let (ball_sqrt2, _flags) = ball_two.sqrt();

    // The correctly rounded scalar must lie inside [lower, upper].
    use core::cmp::Ordering::Greater;
    let lo = ball_sqrt2.lower();
    let hi = ball_sqrt2.upper();
    let scalar_at_or_above_lo = lo.partial_cmp(&scalar_sqrt2).0 != Some(Greater);
    let scalar_at_or_below_hi = scalar_sqrt2.partial_cmp(&hi).0 != Some(Greater);
    assert!(
        scalar_at_or_above_lo && scalar_at_or_below_hi,
        "the scalar must lie inside the rigorous ball enclosure"
    );

    println!("scalar sqrt(2)      = {scalar_sqrt2}");
    println!("ball enclosure lo   = {lo}");
    println!("ball enclosure hi   = {hi}");
    println!("the scalar lies inside the rigorous enclosure");
}
