//! Reaching a target accuracy without guessing a working precision.
//!
//! Run with: `cargo run --example target_accuracy`
//! Companion guide: docs/guides/06-target-accuracy.md
//!
//! You state the relative accuracy you want, in bits, and
//! `refine_to_accuracy` grows the working precision until the ball's
//! certified radius is small enough to deliver it. The computation is
//! written once as a function of the precision `p`; the driver supplies
//! whatever `p` the target demands. `sqrt` is part of the always-available
//! arithmetic surface, so this example needs no extra feature.

use pfloat::BigFloat;
use pfloat_ball::{refine_to_accuracy, Ball};

fn main() {
    // The accuracy we want: 200 bits of certified relative accuracy on
    // sqrt(2). We never name a working precision; the driver finds one.
    let target_bits: i64 = 200;

    let (ball, _status) = refine_to_accuracy(
        target_bits,
        // Start low on purpose, to show the driver climbing.
        32,
        // A ceiling, so an unreachable target cannot loop forever.
        4096,
        |p| {
            // Build the input at the precision the driver handed us, then
            // run the rigorous kernel. The radius shrinks as `p` grows.
            let two = Ball::point(BigFloat::try_from_i64_exact(2, p).unwrap()).unwrap();
            two.sqrt()
        },
    );

    let certified = ball.rel_accuracy_bits();
    println!("requested accuracy : {target_bits} bits");
    println!("certified accuracy : {certified} bits");
    println!("midpoint           : {}", ball.midpoint());
    println!("lower endpoint     : {}", ball.lower());
    println!("upper endpoint     : {}", ball.upper());

    // The driver delivered at least what we asked for.
    assert!(
        certified >= target_bits,
        "the ball must certify at least the requested accuracy"
    );

    // sqrt(2) is irrational, so the enclosure is not a point: it carries a
    // positive radius. The accuracy lives in that radius, not in a flag.
    assert!(!ball.is_exact(), "an irrational root has a positive radius");

    // The enclosure is sound: the true sqrt(2), computed independently at a
    // much higher precision, lies between the endpoints. We bracket it with
    // its square rather than transcribing digits: lower^2 <= 2 <= upper^2.
    use core::cmp::Ordering::Greater;
    use pfloat::RoundingMode::NearestEven;
    let two_hi = BigFloat::try_from_i64_exact(2, 4096).unwrap();
    let lo = ball.lower();
    let hi = ball.upper();
    let (lo_sq, _) = lo.mul(&lo, NearestEven);
    let (hi_sq, _) = hi.mul(&hi, NearestEven);
    // lower^2 <= 2: lo_sq is not greater than two.
    assert!(lo_sq.partial_cmp(&two_hi).0 != Some(Greater));
    // 2 <= upper^2: two is not greater than hi_sq.
    assert!(two_hi.partial_cmp(&hi_sq).0 != Some(Greater));

    println!("the enclosure brackets the true sqrt(2)");
}
