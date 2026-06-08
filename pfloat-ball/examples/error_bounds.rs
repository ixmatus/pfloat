//! Rigorous error bounds with pfloat-ball: a guaranteed enclosure, a
//! radius as the accuracy channel, and exact-in-exact-out.
//!
//! Run with: `cargo run --example error_bounds`
//! Companion guide: docs/guides/04-error-bounds.md
//!
//! Every assertion here pins a claim the guide makes, so a clean run is a
//! proof the walk through is faithful. The example uses only the always
//! available arithmetic surface (`sqrt`, `mul`, `add`), so it needs no
//! `trig` or `exp-log` feature.

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode};
use pfloat_ball::{Ball, Mag};

const NE: RoundingMode = RoundingMode::NearestEven;

fn bf(n: i64, precision: u32) -> BigFloat {
    BigFloat::try_from_i64_exact(n, precision).unwrap()
}

/// `true` when `lo <= x <= hi`, the containment test the FTIA promises.
fn within(lo: &BigFloat, x: &BigFloat, hi: &BigFloat) -> bool {
    lo.partial_cmp(x).0 != Some(Ordering::Greater) && x.partial_cmp(hi).0 != Some(Ordering::Greater)
}

fn main() {
    // A ball you build directly: [1 +/- 2^-50]. The radius is a Mag, an
    // unsigned magnitude that can only round up, so an unsound (inward)
    // radius is unrepresentable.
    let measured = Ball::new(bf(1, 64), Mag::from_pow2(-50)).unwrap();

    // Law 5: the radius is the accuracy channel. rel_accuracy_bits reads
    // it as log2(|mid| / rad): mid = 1 (exponent 0), rad = 2^-50, so the
    // certified relative accuracy is exactly 50 bits.
    assert_eq!(measured.rel_accuracy_bits(), 50);
    assert!(!measured.is_exact());

    // Law 1 (FTIA): every point the ball denotes lies in [lower, upper].
    // The endpoints are exact directed roundings (Law 4), so 1 itself,
    // the midpoint, sits inside the guaranteed enclosure.
    let lo = measured.lower();
    let hi = measured.upper();
    assert!(within(&lo, &bf(1, 64), &hi));

    // A square root over an exact (point) input. sqrt(2) is irrational,
    // so the output ball has a positive radius: it is a guaranteed
    // enclosure of the true sqrt(2), never a claim narrower than the
    // truth.
    let two = Ball::point(bf(2, 200)).unwrap();
    let (root_two, _status) = two.sqrt();
    assert!(
        !root_two.is_exact(),
        "sqrt(2) cannot be exact at finite precision"
    );

    // The FTIA guarantee, checked against a known interior value: a
    // separately and correctly rounded sqrt(2) must lie inside the ball.
    let scalar_root_two = bf(2, 200).sqrt(NE).0;
    let lo = root_two.lower();
    let hi = root_two.upper();
    assert!(
        within(&lo, &scalar_root_two, &hi),
        "the rigorous enclosure must contain the true sqrt(2)"
    );

    // The accuracy channel survives an arithmetic op. Squaring the
    // enclosure of sqrt(2) gives a ball that still contains 2, with a
    // radius that records the propagated error.
    let (squared, _status) = root_two.mul(&root_two);
    assert!(within(&squared.lower(), &bf(2, 200), &squared.upper()));
    println!(
        "sqrt(2)^2 enclosure carries {} bits of certified accuracy",
        squared.rel_accuracy_bits()
    );

    // Law 3 (exact-in-exact-out): a square root of a perfect square is
    // exact, so the directed pair coincides and the output ball has a
    // zero radius. sqrt of the point ball {4} is the point ball {2}.
    let four = Ball::point(bf(4, 64)).unwrap();
    let (root_four, _status) = four.sqrt();
    assert!(root_four.is_exact(), "sqrt(4) = 2 exactly: an exact ball");
    assert_eq!(root_four.rel_accuracy_bits(), i64::MAX);
    // Both endpoints collapse onto 2.
    assert!(within(&root_four.lower(), &bf(2, 64), &root_four.upper()));

    // An exact op on exact balls stays exact: 3 * 7 = 21 with no rounding.
    let three = Ball::point(bf(3, 64)).unwrap();
    let seven = Ball::point(bf(7, 64)).unwrap();
    let (twenty_one, _status) = three.mul(&seven);
    assert!(twenty_one.is_exact(), "3 * 7 = 21 is exact in, exact out");
    assert!(within(
        &twenty_one.lower(),
        &bf(21, 64),
        &twenty_one.upper()
    ));

    println!("all enclosure laws hold: FTIA containment, radius accuracy, exact-in-exact-out");
}
