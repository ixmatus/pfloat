# How to: reach a target accuracy without guessing a precision

Type: how to (task oriented). Time: ten minutes. Runnable companion:
`examples/target_accuracy.rs` in the `pfloat-ball` crate.

You want a result good to a stated number of bits, and you do not want to
guess a working precision that delivers it. Too low and the answer is short
of the mark; too high and you pay for accuracy you never asked for. This
guide states the accuracy up front and lets `refine_to_accuracy` grow the
precision until the ball certifies it.

## Goal

Compute the square root of two to at least 200 bits of certified relative
accuracy, naming the accuracy you want rather than the precision you hope
will reach it. `sqrt` is part of the always available arithmetic surface,
so no extra feature is needed.

## 1. Separate the two numbers

A ball carries two numbers that are easy to conflate. The working
precision is the effort: the bit width the kernel computes at. The
certified accuracy is the result: how tightly the radius pins the true
value down, read off the ball with `rel_accuracy_bits`. The whole point of
this guide is that you set the second and let the driver find the first.

```rust
use pfloat_ball::{refine_to_accuracy, Ball};
use pfloat::BigFloat;

let target_bits: i64 = 200;
```

## 2. Write the computation as a function of precision

`refine_to_accuracy` calls your closure with a precision `p` and expects a
ball back. Build the inputs at `p`, run the rigorous kernel, return the
result and its status. Write it once; the driver decides what `p` to pass.

```rust
let (ball, _status) = refine_to_accuracy(
    target_bits,
    32,   // start precision: deliberately low, to watch the driver climb
    4096, // max precision: a ceiling so an unreachable target stops
    |p| {
        let two = Ball::point(BigFloat::try_from_i64_exact(2, p).unwrap()).unwrap();
        two.sqrt()
    },
);
```

The driver evaluates the closure at the start precision, reads
`rel_accuracy_bits` off the result, and stops the moment it meets the
target. Otherwise it grows the precision geometrically (times 1.5, at least
plus 32 bits) and tries again, so the loop terminates in a logarithmic
number of evaluations. The `max_precision` ceiling is the escape hatch: a
genuinely entire result can never reach the target, and the driver returns
the `max_precision` ball rather than spinning forever.

## 3. Read the certified accuracy

The returned ball meets or beats the target. Read it back and confirm.

```rust
let certified = ball.rel_accuracy_bits();
println!("requested accuracy : {target_bits} bits");
println!("certified accuracy : {certified} bits");
assert!(certified >= target_bits);
```

The certified figure usually overshoots a little, because precision grows
in discrete steps and the step that crosses the target lands past it. The
companion prints 215 bits for a 200 bit request: the first precision the
driver tried that cleared the bar delivered fifteen bits of headroom.

## 4. Remember the accuracy lives in the radius, not a flag

The square root of two is irrational, so the enclosure is not a point: it
has a positive radius, and that radius is where the accuracy lives. On a
ball a small positive radius is the normal correct outcome of an inexact
operation, not a failure. You read accuracy with `rel_accuracy_bits`, not
by checking a status flag.

```rust
assert!(!ball.is_exact()); // an irrational root carries a positive radius
```

To prove the enclosure is sound without transcribing a long digit string,
bracket the true value with its square: square the endpoints and check that
`lower` squared does not exceed two and two does not exceed `upper`
squared. The companion runs exactly this check and it passes.

## What you reached

- You stated an accuracy in bits and the driver found a precision that
  delivers it. The computation was written once as a function of `p`; the
  effort was the driver's to choose.
- The certified accuracy comes off the ball with `rel_accuracy_bits`, and
  it meets or slightly overshoots the target because precision grows in
  discrete steps.
- A `max_precision` ceiling keeps an unreachable target (an entire result)
  from looping forever; the driver returns the best it found.

## Next

- Guide 04 introduces the ball type and the midpoint and radius model that
  `rel_accuracy_bits` reads.
- Guide 05 covers the elementary functions on a ball (`exp`, `ln`, `sin`),
  which compose into the same `refine_to_accuracy` closure when you enable
  the `exp-log` and `trig` features.
