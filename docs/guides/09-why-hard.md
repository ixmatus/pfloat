# Tutorial: why correct rounding and interval soundness are hard

Type: tutorial (learning oriented, an explanation). Time: fifteen minutes.
Runnable companion: `examples/why_hard.rs` in the `pfloat` crate.

This guide explains the two problems the rest of the family exists to solve.
The first is the table maker's dilemma: deciding the last bit of a correctly
rounded transcendental can demand far more working precision than the answer
itself occupies. The second is soundness: a guaranteed enclosure of a value is
not a single number, and the bound that holds it must round outward, away from
the truth, never toward it. You will pin a transcendental between two adjacent
`f64` grid points with a directed pair, watch a naive midpoint fail to be a
bound, and leave knowing why `pfloat-ball` stores a center and a radius and why
that radius always rounds up.

## Step 1: a transcendental lands between two grid points

A finite floating format is a grid. The representable `f64` values are spaced
points on the real line, and the spacing widens as the magnitude grows. The
square root of two is irrational, so it is not one of those points. It falls
strictly inside the gap between two of them.

```rust
use pfloat::{BigFloat, RoundingMode};

const NE: RoundingMode = RoundingMode::NearestEven;

let two = BigFloat::try_from_i64_exact(2, 300).unwrap();
let (root, status) = two.sqrt(NE);
assert!(status.inexact());
```

Three hundred bits of working precision resolve the true value far past the
fifty three bits an `f64` carries. The `INEXACT` flag confirms what irrationality
already tells you: no finite format holds this value exactly. The question correct
rounding has to answer is which grid point to report, and that question is harder
than it looks.

## Step 2: bracket the value with a directed pair

Round the same high precision value two ways. Rounding toward minus infinity
gives the grid point just below the true value; rounding toward plus infinity
gives the grid point just above. Together they bracket it.

```rust
const TN: RoundingMode = RoundingMode::TowardNegative;
const TP: RoundingMode = RoundingMode::TowardPositive;

let (lo, lo_status) = root.to_f64_round(TN);
let (hi, hi_status) = root.to_f64_round(TP);
assert!(lo_status.inexact());
assert!(hi_status.inexact());

assert!(lo < hi);
let ulp_gap = hi.to_bits() - lo.to_bits();
assert_eq!(ulp_gap, 1, "the bracket must be exactly one ULP wide");
```

The two grid points are adjacent: their bit patterns differ by one, so they are
exactly one unit in the last place apart. The true value sits somewhere in the
open interval between them, and both conversions report `INEXACT` because neither
endpoint is the value itself. This pair, one end rounded down and the other
rounded up, is the simplest sound enclosure there is.

## Step 3: the table maker's dilemma

Correct rounding to nearest has to decide which of `lo` and `hi` is closer.

```rust
let (ne, _) = root.to_f64_round(NE);
assert!(ne == lo || ne == hi);
assert_eq!(ne, 2.0_f64.sqrt());
```

For the square root of two the decision is easy, because the true value is not
near the halfway line between the grid points. The hard cases are the ones where
it is. A transcendental can land so close to the exact midpoint of two grid
points that you cannot tell which side it falls on until you have computed it to
many more bits than the result occupies. There is no general formula for how
many; you discover the bound by computing wider and wider until the answer stops
moving. That open ended demand for precision is the table maker's dilemma, named
for the eighteenth century compilers of logarithm tables who first hit it. Ziv's
strategy answers it operationally: compute at a trial precision, test whether the
result is far enough from a rounding boundary to be certain, and if it is not,
raise the precision and retry. The kernels in `pfloat` run that loop so a caller
never sees a wrong last bit.

The dilemma is why correct rounding is hard. The cost is unbounded in the worst
case, the worst case is rare and unpredictable, and getting it wrong is silent:
a single wrong bit looks exactly like a right one until something downstream
disagrees.

## Step 4: a midpoint is not a bound

Suppose you keep only the nearest grid point. You have thrown away the gap. The
nearest point is one number, and one number cannot say "the answer is somewhere
in here". Watch what happens if you try to recover the gap by averaging the two
endpoints.

```rust
let naive_midpoint = f64::midpoint(lo, hi);
assert!(naive_midpoint == lo || naive_midpoint == hi);
```

The exact mathematical midpoint of two adjacent grid points sits exactly halfway
between them, which is itself not representable, so it rounds back onto a grid
point. The average of `lo` and `hi` is `lo` again. You are back to a single
point with no width. A correctly rounded value answers "what is the closest
representable number"; it does not answer "what interval is the truth guaranteed
to lie in". Those are different questions, and the second one needs more than a
point to answer.

## Step 5: an enclosure rounds outward

To carry the gap you need two numbers. The directed pair `[lo, hi]` is one such
representation. `pfloat-ball` uses the other: a center and a nonnegative radius,
a `Ball` whose midpoint is a `BigFloat` and whose radius is a `Mag`. The interval
it denotes is the midpoint plus or minus the radius, and the rule that makes the
ball sound is the direction the radius rounds.

The radius must round outward, which here means up. A `Mag` is a nonnegative
magnitude that rounds up by construction, so every operation that widens a ball
overestimates the error rather than underestimating it. If a radius ever rounded
down, the stored interval could be narrower than the truth, and a value the ball
claims to enclose could in fact sit just outside it. Outward rounding buys
soundness: the interval is sometimes looser than necessary, never tighter than
correct. The directed pair obeys the same rule from both ends, rounding the lower
endpoint down and the upper endpoint up, so the bracket can only grow, never
shrink below the truth.

```rust
let true_value_below_hi = root.to_f64_round(TP).0 == hi;
let true_value_above_lo = root.to_f64_round(TN).0 == lo;
assert!(true_value_below_hi && true_value_above_lo);
```

That asymmetry, tighten by computing harder but widen the bound rather than risk
excluding the truth, is the whole discipline of rigorous numerics. Correct
rounding spends precision to nail the last bit of a point. Interval soundness
spends a little width to guarantee containment. The two problems are why the
family splits into a correctly rounded scalar layer (`pfloat`) and a verified
enclosure layer (`pfloat-ball`).

## What you learned

- A transcendental falls strictly between two grid points. A directed pair, the
  same value rounded down and rounded up, brackets it; for the square root of two
  the two endpoints are exactly one ULP apart.
- Deciding the last bit to nearest can demand far more precision than the result
  occupies, with no fixed bound in advance. That is the table maker's dilemma,
  and Ziv's compute, test, retry loop is how the kernels answer it.
- A single rounded number is a point, not an interval. A sound enclosure carries
  width, and the bound that carries it rounds outward (a `Mag` rounds up by
  construction) so the interval is never tighter than the truth.

## Next

- Guide 3 walks the five rounding modes in full, including when a directed mode
  is the right tool rather than nearest.
- Guide 4 builds a `Ball` with `pfloat-ball` and reads the certified enclosure
  the radius defines.
- Guide 6 reaches a target accuracy with `refine_to_accuracy`, growing the
  working precision until the enclosure is tight enough, the constructive answer
  to the dilemma this guide describes.