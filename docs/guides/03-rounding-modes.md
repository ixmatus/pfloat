# Tutorial: choosing a rounding mode, and why directed modes matter

Type: tutorial (learning oriented). Time: ten minutes. Runnable companion:
`examples/rounding_modes.rs` in the `pfloat` crate.

By the end you will apply all five rounding modes to one irrational result, the
square root of two, and read the bit patterns they produce. You will see that
the two directed modes bracket the true value from below and above, that the
bracket is as tight as the grid allows, and that nearest is the default you keep
until you need a guaranteed side. Default features are enough; sqrt(2) needs no
extra feature flag.

## Step 1: compute the value once, at the working precision

A rounding mode is the choice you make when a real result has to land on a finite
grid. So first compute the value with headroom, then round it five different
ways. Two hundred fifty six bits pins sqrt(2) far past the f64 grid we will round
onto.

```rust
use pfloat::{BigFloat, RoundingMode};

let two = BigFloat::try_from_i64_exact(2, 256).unwrap();
let (root, status) = two.sqrt(RoundingMode::NearestEven);
assert!(status.inexact());
```

The `INEXACT` flag confirms what makes this tutorial interesting: sqrt(2) is
irrational, so it never sits exactly on a finite grid, and the mode you pick
genuinely changes the answer. An exact result (a perfect square) would round the
same way under every mode, and there would be nothing to choose.

## Step 2: apply all five modes

`RoundingMode` enumerates the five IEEE 754-2019 rounding attributes. Round the
high precision `root` onto the f64 grid under each one. f64 is a fixed grid, so
"the next representable value" is well defined and the bit patterns are directly
comparable.

```rust
let (ne, _) = root.to_f64_round(RoundingMode::NearestEven);
let (na, _) = root.to_f64_round(RoundingMode::NearestAway);
let (tz, _) = root.to_f64_round(RoundingMode::TowardZero);
let (tp, _) = root.to_f64_round(RoundingMode::TowardPositive);
let (tn, _) = root.to_f64_round(RoundingMode::TowardNegative);
```

`NearestEven` rounds to the closest grid value, breaking ties to the even last
bit; it is the IEEE 754 default and the one you keep most of the time.
`NearestAway` is the same except ties go away from zero. The three remaining
modes are the directed ones: `TowardZero` truncates, `TowardPositive` rounds up
toward positive infinity, and `TowardNegative` rounds down toward negative
infinity.

## Step 3: read the directed bracket

The two directed modes are the reason this tutorial exists. `TowardNegative`
gives the largest grid value at or below the true result, and `TowardPositive`
gives the smallest grid value at or above it. Together they bracket sqrt(2).

```rust
assert!(tn <= ne);
assert!(ne <= tp);
```

The true value of sqrt(2) lies in the closed interval `[tn, tp]`. This is a
guarantee, not an estimate: no finite computation can place the result outside
that bracket, because each endpoint was rounded in the direction that cannot
overshoot it. When you need a value you can prove the answer is below, ask for
`TowardPositive`; when you need one you can prove it is above, ask for
`TowardNegative`. Nearest gives you the closest value but no side guarantee.

## Step 4: see how tight the bracket is

For an irrational like sqrt(2), the directed endpoints land on two adjacent f64
values, one ULP apart. "ULP" is the unit in the last place, the gap between
consecutive grid values. Consecutive finite f64 values of the same sign have
consecutive bit patterns, so the bracket width shows up as a difference of one
in the raw bits.

```rust
assert_eq!(tp.to_bits() - tn.to_bits(), 1);
```

That is the tightest a bracket can be on this grid. You cannot enclose an
irrational value in anything narrower than the two grid points straddling it.

```rust
// sqrt(2) is positive, so toward zero rounds the same direction as toward
// negative infinity: both step down to the smaller neighbor.
assert_eq!(tz.to_bits(), tn.to_bits());
```

`TowardZero` agreed with `TowardNegative` here only because sqrt(2) is positive;
toward zero means down for a positive value and up for a negative one. For a
negative result the two would swap, which is exactly why both of the directed
infinity modes exist alongside truncation.

## Step 5: the same bracket in decimal

The bit level bracket is also visible in decimal once you print enough digits to
reach the place where the two directed results diverge.

```rust
let down = root.to_decimal_string(40, RoundingMode::TowardNegative);
let up = root.to_decimal_string(40, RoundingMode::TowardPositive);
assert!(down.starts_with("1.41421356237309504880168872420"));
assert!(up.starts_with("1.41421356237309504880168872420"));
assert_ne!(down, up);
```

Both strings agree on the leading digits everyone can prove and part ways only
deep in the tail. That shared prefix is the set of digits the bracket pins down;
the divergence marks where the uncertainty begins. A tighter working precision
would push the divergence further right.

## The mental model

Three ideas carry the rounding mode story.

- A rounding mode is a late choice, applied when a real value lands on a finite
  grid; it changes the answer only when the result is inexact. `NearestEven` is
  the IEEE 754 default you keep until you have a reason not to.
- The two directed modes give a guaranteed bracket. `TowardNegative` is a proven
  lower bound and `TowardPositive` a proven upper bound, so the true value sits
  in `[tn, tp]` with no overshoot on either side.
- The bracket is as tight as the grid: for an irrational result the directed
  endpoints are one ULP apart, the narrowest enclosure two grid points allow.

## Next

- Tutorial 1 covers precision as a bit count and the `(value, Status)` pair, if
  you skipped straight here.
- Tutorial 2 reads the full status surface and the thread local flag accumulator.
- Tutorial 4 moves from a one ULP bracket on a single value to a rigorous
  enclosure that tracks error through a whole computation, with `pfloat-ball`.
