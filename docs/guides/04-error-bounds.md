# Tutorial: rigorous error bounds with pfloat-ball

Type: tutorial (learning oriented). Time: fifteen minutes. Runnable companion:
`examples/error_bounds.rs` in the `pfloat-ball` crate.

By the end you will build a `Ball`, run a square root and an arithmetic op over
it, and read a result that is not one rounded number but a guaranteed enclosure:
an interval that provably contains the true answer. You will learn the three
laws that make this rigorous, that the radius is the accuracy channel rather
than an afterthought, and that an exact input gives an exact output.

A `Ball` differs from the correctly rounded `BigFloat` of tutorial 1 in what it
claims. A `BigFloat` is the single closest representable value, paired with a
flag telling you it rounded. A `Ball` is a midpoint and a radius, `[mid ± rad]`,
denoting the closed interval `[mid − rad, mid + rad]`, and it promises the true
mathematical result lies somewhere inside. Reach for a ball when you need a
proof of containment, not a best guess.

## Step 1: build a ball and read its accuracy

A ball pairs a full precision midpoint with a radius. The radius is a `Mag`, an
unsigned magnitude that rounds only upward by construction, so an unsound
(inward) radius is a value you cannot even write. Here is a ball denoting a value
known to about fifty bits.

```rust
use pfloat::{BigFloat, RoundingMode};
use pfloat_ball::{Ball, Mag};

let measured = Ball::new(bf(1, 64), Mag::from_pow2(-50)).unwrap();

assert_eq!(measured.rel_accuracy_bits(), 50);
assert!(!measured.is_exact());
```

`rel_accuracy_bits()` reads the radius as roughly `log2(|mid| / rad)`. The
midpoint is one (binary exponent zero) and the radius is `2^-50`, so the
certified relative accuracy is exactly fifty bits. This is Law 5 of the
enclosure spec: the radius is the primary accuracy channel. On a ball, `INEXACT`
is the normal correct outcome, so you read the radius, not the status flag, to
learn how good the answer is.

## Step 2: enclose a square root, and check the FTIA guarantee

Build an exact ball from a point, then take its square root. Because the square
root of two is irrational, the output ball carries a positive radius. That is
the point: the ball does not pretend to a value it cannot represent; it brackets
the truth.

```rust
let two = Ball::point(bf(2, 200)).unwrap();
let (root_two, _status) = two.sqrt();
assert!(!root_two.is_exact());

let scalar_root_two = bf(2, 200).sqrt(RoundingMode::NearestEven).0;
let lo = root_two.lower();
let hi = root_two.upper();
assert!(within(&lo, &scalar_root_two, &hi));
```

`lower()` and `upper()` give the endpoints `mid − rad` and `mid + rad`, each
rounded outward (Law 4: ball to endpoints is exact, the tightest representable
bracket). The assertion is the Fundamental Theorem of Interval Arithmetic
(Law 1) made concrete: a separately and correctly rounded square root of two, a
known interior value, must lie inside `[lower, upper]`. Soundness runs one way
only. The radius may over estimate and the enclosure stays correct; it must
never under estimate, because an under estimating radius turns every downstream
bound into a falsehood.

## Step 3: carry the bound through an arithmetic op

The accuracy channel survives composition. Square the enclosure of the square
root of two and the result still contains two, with a radius that records the
error propagated through the multiply.

```rust
let (squared, _status) = root_two.mul(&root_two);
assert!(within(&squared.lower(), &bf(2, 200), &squared.upper()));
println!("certified accuracy: {} bits", squared.rel_accuracy_bits());
```

Every ball op computes its midpoint with pfloat's correctly rounded kernel and
sets the radius from the rounding error those kernels already report (Law 2, the
directed pair route). You never reason about the error budget by hand. You ask
the result for `rel_accuracy_bits()` and it tells you how many bits it still
vouches for.

## Step 4: exact in gives exact out

A square root of a perfect square is exact, so there is nothing to round. The
directed pair coincides, the radius is zero, and the ball denotes a single
point. This is Law 3: exactness in produces exactness out, with no spurious
slack.

```rust
let four = Ball::point(bf(4, 64)).unwrap();
let (root_four, _status) = four.sqrt();
assert!(root_four.is_exact());
assert_eq!(root_four.rel_accuracy_bits(), i64::MAX);

let three = Ball::point(bf(3, 64)).unwrap();
let seven = Ball::point(bf(7, 64)).unwrap();
let (twenty_one, _status) = three.mul(&seven);
assert!(twenty_one.is_exact());
```

The square root of the point ball `{4}` is the point ball `{2}`: `is_exact()`
holds and `rel_accuracy_bits()` returns `i64::MAX`, the marker for infinite
relative accuracy. The same holds for `3 * 7 = 21`, which rounds nowhere. An
exact ball is the rigorous counterpart of a `BigFloat` whose `inexact()` flag
stayed clear.

## What you learned

Three laws carry through the rest of the rigorous tower.

- A ball is a guaranteed enclosure. `[mid ± rad]` contains the true result, and
  `within(lower, value, upper)` is the FTIA containment you can check (Law 1).
  Soundness is one directional: the radius may widen, never narrow.
- The radius is the accuracy channel. `rel_accuracy_bits()` reports certified
  relative accuracy off the radius, and a small positive radius is the normal
  correct outcome of an inexact op, not a failure (Law 5).
- Exact in gives exact out. A perfect square root or an exact product yields a
  zero radius ball, `is_exact()` true, accuracy `i64::MAX` (Law 3).

## Next

- Tutorial 5 drives precision automatically with `refine_to_accuracy`, which
  re evaluates a computation at growing working precision until the certified
  radius reaches a target bit count.
- Tutorial 6 reaches the ball elementary functions (`exp`, `ln`, `sin`),
  enabled by the `exp-log` and `trig` features.
- Tutorial 1 covers the scalar `BigFloat` surface a ball midpoint is built from,
  if you skipped straight to enclosures.
