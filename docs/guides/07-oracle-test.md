# How to: test your own numerical code against a trusted oracle

Type: how to (task oriented). Time: fifteen minutes. Runnable companion:
`examples/oracle_test.rs` in the `pfloat-ball` crate.

You have an `f64` routine and you want to know whether it is right. A `Ball`
turns that question into a pass/fail test with a guaranteed answer: it encloses
the true result in a proven bracket, and you check whether your value lands
inside. No reference table, no second implementation to trust, no tolerance you
picked by feel. This guide wires a ball oracle around `f64::sqrt`, runs it over
a list of inputs, and reports the smallest input that fails.

## Goal

Build a rigorous `f64` bracket for each input, run the code under test, and
flag any result that falls outside the bracket, surfacing the smallest failing
input first.

## 1. Name the function under test

The oracle does not care how your routine is built. It takes an `f64` in and an
`f64` out and treats the body as opaque. Here it is the platform square root;
substitute your own kernel and nothing else changes.

```rust
fn function_under_test(x: f64) -> f64 {
    x.sqrt()
}
```

## 2. Build the rigorous bracket from a ball

`Ball::point(input)` is the exact input as a zero radius ball. Applying
`sqrt()` returns a ball that provably contains the true square root: the
midpoint is the rounded value and the radius bounds the residual. The bracket
is part of the always available arithmetic surface, so the default features
suffice and no rounding mode enters the ball call.

To compare against an `f64`, push both ball endpoints onto the `f64` grid in
the outward direction. Round the lower endpoint toward negative infinity and
the upper endpoint toward positive infinity; both steps move away from the
truth, so the resulting pair of `f64` values is guaranteed to straddle the true
square root.

```rust
const TOWARD_NEG: RoundingMode = RoundingMode::TowardNegative;
const TOWARD_POS: RoundingMode = RoundingMode::TowardPositive;
const ORACLE_BITS: u32 = 200;

fn sqrt_bracket(n: i64) -> (f64, f64) {
    let x = BigFloat::try_from_i64_exact(n, ORACLE_BITS).unwrap();
    let (ball, _flags) = Ball::point(x).unwrap().sqrt();
    let (lo, _) = ball.lower().to_f64_round(TOWARD_NEG);
    let (hi, _) = ball.upper().to_f64_round(TOWARD_POS);
    (lo, hi)
}
```

The working precision sits well above the 53 bits of an `f64`, so the ball is
tighter than one `f64` step. The outward rounded bracket therefore spans at
most the two `f64` values on either side of the truth. For a perfect square the
true root is itself an `f64`, both endpoints round to it, and the bracket
collapses to a single value.

## 3. Check the result and report the smallest failure

A result passes when it lands in the closed bracket. Because the bracket is two
plain `f64` values, the test is two `f64` comparisons. Sweep the inputs and
keep the smallest one that fails, which is the case you most want to look at
first.

```rust
fn first_failure(kernel: fn(f64) -> f64, inputs: &[i64]) -> Option<i64> {
    let mut smallest: Option<i64> = None;
    for &n in inputs {
        let (lo, hi) = sqrt_bracket(n);
        let got = kernel(n as f64);
        if !(lo <= got && got <= hi) {
            smallest = Some(match smallest {
                Some(prev) if prev <= n => prev,
                _ => n,
            });
        }
    }
    smallest
}
```

Run it over the inputs and a clean sweep returns `None`:

```rust
let inputs = [0_i64, 1, 2, 3, 4, 5, 9, 16];
let result = first_failure(function_under_test, &inputs);
assert_eq!(result, None);
println!("function_under_test: PASS on {} inputs", inputs.len());
```

`f64::sqrt` is correctly rounded, so it lands inside every bracket. For an
irrational case the bracket holds the two adjacent `f64` values that straddle
the truth, and the result is one of them. For a perfect square the bracket is a
single `f64` and the result lands on it.

## 4. Watch the oracle catch a real bug

A pass on a correct kernel proves nothing about the oracle's teeth. Point it at
a kernel you know is wrong. One Newton step from a poor seed is exact at the
fixed points 0 and 1 and wrong everywhere else, so the oracle should report 2
as the smallest failing input.

```rust
fn buggy_sqrt(x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    let guess = x; // a poor seed
    0.5 * (guess + x / guess) // one Newton step, not converged
}

let bug = first_failure(buggy_sqrt, &inputs);
assert_eq!(bug, Some(2));
```

At input 2 the broken kernel returns 1.5; the bracket is roughly
`[1.41421356237309492, 1.41421356237309515]`. The result sits far outside, the
oracle flags it, and the smallest failing input points you straight at the
easiest case to debug.

## The shape of the technique

- A ball is a proven bracket, not an estimate. The true result is inside the
  enclosure by construction, so a result outside it is a bug with a proof, not
  a value that merely looks off.
- Round the endpoints outward to the target type. Lower toward negative
  infinity, upper toward positive infinity; the bracket can only widen, never
  exclude the truth, and the comparison stays in the target type.
- Report the smallest failing input. A sweep that returns the minimal counter
  example hands you the simplest reproducer rather than whichever case the loop
  happened to reach first.

## Next

- Guide 06 reaches a target accuracy without guessing a precision, using
  `refine_to_accuracy`; pair it with this oracle to grow the bracket until it
  certifies the bit you care about.
- Guide 04 introduces the `Ball` type and the certified enclosure from the
  ground up, if the ball mechanics here moved fast.
- Guide 09 explains why a midpoint is not a bound and why a sound radius must
  round outward, the reasoning under the outward rounding step above.