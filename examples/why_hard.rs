//! Tutorial 9: why correct rounding and interval soundness are hard.
//!
//! Run with: `cargo run --example why_hard`
//!
//! Companion to `docs/guides/09-why-hard.md`. The square root of two is
//! irrational, so at any finite output format it falls strictly between two
//! representable grid points. This example pins that gap with a directed pair:
//! `to_f64_round(TowardNegative)` and `to_f64_round(TowardPositive)` return the
//! two `f64` grid points that bracket the true value, exactly one unit in the
//! last place apart. It then shows why the nearest grid point alone is not a
//! bound, which is the seed of the table maker's dilemma and the reason a sound
//! enclosure must round its radius outward.

use pfloat::{BigFloat, RoundingMode};

const NE: RoundingMode = RoundingMode::NearestEven;
const TN: RoundingMode = RoundingMode::TowardNegative;
const TP: RoundingMode = RoundingMode::TowardPositive;

fn main() {
    // Compute sqrt(2) with generous headroom: 300 bits is far more than the 53
    // an f64 holds, so the high precision value resolves the true result well
    // past the f64 grid spacing.
    let two = BigFloat::try_from_i64_exact(2, 300).unwrap();
    let (root, status) = two.sqrt(NE);

    // Irrational result: no finite format represents it exactly.
    assert!(status.inexact());

    // The directed pair. Rounding the SAME high precision value toward minus
    // infinity and toward plus infinity yields the two f64 grid points that the
    // true value lies between. Both conversions report inexact, because the true
    // value sits strictly inside the open interval (lo, hi).
    let (lo, lo_status) = root.to_f64_round(TN);
    let (hi, hi_status) = root.to_f64_round(TP);
    assert!(lo_status.inexact());
    assert!(hi_status.inexact());

    // lo is strictly below hi, and they are adjacent: exactly one ULP apart.
    // Adjacency is the whole point. The true value has nowhere to hide between
    // them, yet neither grid point equals it.
    assert!(lo < hi);
    let ulp_gap = hi.to_bits() - lo.to_bits();
    assert_eq!(ulp_gap, 1, "the bracket must be exactly one ULP wide");

    // Nearest even lands on one of the two grid points (here, the upper one),
    // and it agrees with the hardware sqrt. Nearest gives you the closer grid
    // point; it does not give you a bound.
    let (ne, _) = root.to_f64_round(NE);
    assert!(ne == lo || ne == hi);
    assert_eq!(ne, 2.0_f64.sqrt());

    println!("sqrt(2) is bracketed by two adjacent f64 grid points:");
    println!(
        "  lo (round down) = {lo:.20}  bits = {:#018x}",
        lo.to_bits()
    );
    println!(
        "  hi (round up)   = {hi:.20}  bits = {:#018x}",
        hi.to_bits()
    );
    println!("  the two differ by exactly {ulp_gap} ULP, and the true value lies between them");

    // Why a midpoint is not a bound. Take the midpoint of the two grid points in
    // f64. The exact mathematical midpoint sits exactly halfway between two
    // adjacent grid points, is itself not representable, and rounds back onto a
    // grid point (here, lo). A single rounded number cannot carry "the answer is
    // somewhere in this gap"; it can only be one point, and one point is not an
    // interval.
    let naive_midpoint = f64::midpoint(lo, hi);
    assert!(naive_midpoint == lo || naive_midpoint == hi);
    println!(
        "  the f64 average of lo and hi collapses back onto a grid point: {}",
        if naive_midpoint == lo { "lo" } else { "hi" }
    );

    // To carry the gap you need two numbers (the pair lo and hi), or a center
    // plus a nonnegative radius. The radius is the topic of guide 4 (pfloat-ball)
    // and guide 6 (refine_to_accuracy); for soundness it must round OUTWARD, so
    // the stored bound is never tighter than the truth. The directed pair above
    // is the simplest sound enclosure: round one end down, the other end up, and
    // the true value is guaranteed inside.
    let true_value_below_hi = root.to_f64_round(TP).0 == hi;
    let true_value_above_lo = root.to_f64_round(TN).0 == lo;
    assert!(true_value_below_hi && true_value_above_lo);
    println!("the directed pair [lo, hi] is a sound enclosure of sqrt(2)");
}
