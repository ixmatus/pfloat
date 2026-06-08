# Tutorial: your first correctly rounded value to N digits

Type: tutorial (learning oriented). Time: ten minutes. Runnable companion:
`examples/first_value.rs` in the `pfloat` crate.

By the end you will compute the square root of two to fifty correct decimal
digits, know why the result is paired with a status, and hold the one mental
model the rest of the family builds on: precision is a bit count you choose up
front, and output digits are a separate choice you make at the end.

## Step 1: choose a precision

A `BigFloat` carries its precision at runtime, as a number of bits. You pick it
when you build the value. Two hundred bits holds about sixty decimal digits, so
it leaves headroom above the fifty digits we will print.

```rust
use pfloat::{BigFloat, RoundingMode};

const NE: RoundingMode = RoundingMode::NearestEven;

let two = BigFloat::try_from_i64_exact(2, 200).unwrap();
```

`try_from_i64_exact(2, 200)` builds the integer two as a two hundred bit value.
It is fallible because not every integer fits exactly in every precision; two
does, so the `unwrap` is safe here.

## Step 2: compute, and read the status

Every kernel returns the value paired with a `Status`, the IEEE 754-2019 sticky
exception flags the operation raised. You read the value and the flags together.

```rust
let (root, status) = two.sqrt(NE);
assert!(status.inexact());
```

The square root of two is irrational, so no finite precision represents it
exactly: the kernel rounds, and the `INEXACT` flag is set. This is the normal
outcome for a transcendental or an irrational result. The flag is how you tell
an exact answer (a perfect square, say) from a rounded one, per call, without
any global state. The rounding mode you passed, `NearestEven`, is the IEEE
default; step into tutorial 3 to see the other four.

## Step 3: print to the digits you want

The working precision and the output digits are two separate choices. You
computed at two hundred bits; now ask for exactly fifty significant digits,
correctly rounded.

```rust
let fifty = root.to_decimal_string(50, NE);
println!("sqrt(2) to 50 digits = {fifty}");
// 1.4142135623730950488016887242096980785696718753769
```

`to_decimal_string(digits, mode)` rounds the high precision value to the
requested number of significant digits under the mode you pass. Because you
computed with headroom above fifty digits, every one of those fifty is correct:
the rounding at the end lands on the true value, not on accumulated error.

## Step 4: let the formatter pick the shortest form

Often you do not want a fixed digit count; you want the shortest string that
reads back as the same value. That is what `to_shortest_decimal_string` gives,
using the Steele and White shortest decimal algorithm.

```rust
let shortest = root.to_shortest_decimal_string();
println!("sqrt(2) shortest = {shortest}");
// 1.41421356237309504880168872420969807856967187537694807317668
```

The shortest form runs to about sixty digits here, because that is how many the
two hundred bit value pins down. Compute at a higher precision and the shortest
form grows; compute at a lower one and it shrinks. The string always round
trips: parse it back and you recover the same `BigFloat`.

## The mental model

Three ideas carry through the rest of the family.

- Precision is a bit count, chosen when you build a value, and it sets how much
  the result pins down. Pick it with headroom above the digits you need.
- Every operation returns `(value, Status)`. The status is the per call record
  of what the operation rounded or flagged; `inexact()` is the one you read
  most.
- Output digits are a separate, late choice. `to_decimal_string(n, mode)` for a
  fixed count, `to_shortest_decimal_string` for the shortest round tripping
  form.

## Next

- Tutorial 2 reads the full status surface: `OVERFLOW`, `UNDERFLOW`, `INVALID`,
  `DIV_BY_ZERO`, and the thread local accumulator.
- Tutorial 3 walks the five rounding modes and shows when a directed mode is the
  right tool.
- Tutorial 4 moves from a single rounded value to a rigorous enclosure with
  `pfloat-ball`, when you need a guaranteed bound rather than a correctly
  rounded point.
