//! Test your own numerical code against a rigorous ball oracle.
//!
//! Run with: `cargo run --example oracle_test`
//! Companion guide: docs/guides/07-oracle-test.md
//!
//! A `Ball` gives you a pass/fail oracle with a guaranteed bracket. For each
//! input we build the rigorous enclosure `Ball::point(input).sqrt()`, which
//! provably contains the true square root. We round its lower endpoint down
//! and its upper endpoint up onto the `f64` grid, giving a bracket of `f64`
//! values that is guaranteed to straddle the truth. A result inside the
//! bracket is certified correct to the last bit; a result outside is a proven
//! bug. The harness reports the smallest failing input, so the first thing you
//! see is the easiest case to debug.

use pfloat::{BigFloat, RoundingMode};
use pfloat_ball::Ball;

const TOWARD_NEG: RoundingMode = RoundingMode::TowardNegative;
const TOWARD_POS: RoundingMode = RoundingMode::TowardPositive;

/// Working precision for the oracle. Generous headroom above the 53 bits of an
/// `f64` keeps the enclosure tighter than one `f64` step, so the outward
/// rounded bracket spans at most the two `f64` values straddling the truth.
const ORACLE_BITS: u32 = 200;

/// The function under test: the platform's `f64` square root. Swap this for
/// your own kernel; the harness does not care how it is implemented.
fn function_under_test(x: f64) -> f64 {
    x.sqrt()
}

/// A deliberately broken kernel: one Newton step from a poor seed. It is exact
/// at the fixed points 0 and 1 but wrong everywhere else, so the oracle should
/// catch it and report 2 as the smallest failing input.
fn buggy_sqrt(x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    let guess = x; // a poor seed: far off for x > 1
    0.5 * (guess + x / guess) // a single Newton step, not converged
}

/// The rigorous `f64` bracket for `sqrt(n)`: round the ball's lower endpoint
/// toward negative infinity and its upper endpoint toward positive infinity,
/// so both rounding steps push outward and the returned pair is guaranteed to
/// enclose the true square root.
fn sqrt_bracket(n: i64) -> (f64, f64) {
    let x = BigFloat::try_from_i64_exact(n, ORACLE_BITS).unwrap();
    let (ball, _flags) = Ball::point(x).unwrap().sqrt();
    let (lo, _) = ball.lower().to_f64_round(TOWARD_NEG);
    let (hi, _) = ball.upper().to_f64_round(TOWARD_POS);
    (lo, hi)
}

/// Run `kernel` against the ball oracle over `inputs`. Returns the smallest
/// failing input, or `None` if every input passed.
fn first_failure(kernel: fn(f64) -> f64, inputs: &[i64]) -> Option<i64> {
    let mut smallest: Option<i64> = None;
    for &n in inputs {
        let (lo, hi) = sqrt_bracket(n);
        let got = kernel(n as f64);
        let inside = lo <= got && got <= hi;
        if !inside {
            // Track the smallest failing input, not just the first one seen.
            smallest = Some(match smallest {
                Some(prev) if prev <= n => prev,
                _ => n,
            });
        }
    }
    smallest
}

fn main() {
    let inputs: [i64; 8] = [0, 1, 2, 3, 4, 5, 9, 16];

    // Pass 1: the real kernel. The platform sqrt is correctly rounded, so
    // every result must lie inside the rigorous bracket.
    let result = first_failure(function_under_test, &inputs);
    assert_eq!(
        result, None,
        "f64::sqrt must pass the ball oracle on every input"
    );
    println!("function_under_test: PASS on {} inputs", inputs.len());

    // Show the bracket for one irrational case. The two endpoints are adjacent
    // f64 values straddling the true sqrt(2), and our result sits between them.
    let (lo2, hi2) = sqrt_bracket(2);
    let got2 = function_under_test(2.0);
    println!("  sqrt(2): bracket [{lo2:.17}, {hi2:.17}]");
    println!("           result   {got2:.17} lies inside");
    assert!(lo2 <= got2 && got2 <= hi2);

    // A perfect square gives an exact ball: the bracket collapses to a single
    // f64 and an exact result lands on it.
    let (lo9, hi9) = sqrt_bracket(9);
    assert_eq!(lo9, hi9, "sqrt(9) is exact: the bracket is a single f64");
    assert_eq!(lo9, 3.0);
    println!("  sqrt(9): bracket collapses to the single f64 {lo9}");

    // Pass 2: the broken kernel. The oracle must catch it and report the
    // smallest input where it leaves the bracket.
    let bug = first_failure(buggy_sqrt, &inputs);
    match bug {
        Some(n) => println!("buggy_sqrt: FAIL, smallest failing input = {n}"),
        None => panic!("the broken kernel should have failed the oracle"),
    }
    // One Newton step is exact at the fixed points 0 and 1, so the smallest
    // input it gets wrong is 2.
    assert_eq!(bug, Some(2), "buggy_sqrt first fails at input 2");

    println!("the ball oracle gives a proven pass/fail bracket, not an estimate");
}
