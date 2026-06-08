# Tutorial: read and act on the status flags

Type: tutorial (learning oriented). Time: ten minutes. Runnable companion:
`examples/status_flags.rs` in the `pfloat` crate.

By the end you will read the full status surface that tutorial 1 only glimpsed.
You will see an inexact result and an exact one side by side, force an
`OVERFLOW`, and use the thread local accumulator that gives you IEEE 754's
"global flags" model when you want it. The one idea to carry away: a status is
both a per call return value and, under `std`, a sticky thread local you can
union across a whole computation.

## Step 1: read an inexact result

Tutorial 1 met `inexact()` once. Here is the pattern again, so the contrast in
the next step lands. Every kernel returns the value paired with a `Status`, the
five IEEE 754-2019 sticky exception flags the call raised.

```rust
use pfloat::{flags, BigFloat, RoundingMode};

const NE: RoundingMode = RoundingMode::NearestEven;

let two = BigFloat::try_from_i64_exact(2, 200).unwrap();
let (_root2, s_inexact) = two.sqrt(NE);
assert!(s_inexact.inexact());
assert!(!s_inexact.overflow());
assert!(!s_inexact.invalid());
```

The square root of two is irrational, so the kernel rounds and sets `INEXACT`.
Nothing else fired: no overflow, no invalid operation. The status predicates are
all `const fn` returning `bool`, so you read them with no allocation and no
ceremony: `inexact()`, `overflow()`, `underflow()`, `invalid()`,
`div_by_zero()`, and `is_ok()`.

## Step 2: read an exact result

Now compute a square root that lands exactly, and watch the flags stay clear.
The square root of four is two, representable at any precision, so the kernel
returns it untouched.

```rust
let four = BigFloat::try_from_i64_exact(4, 200).unwrap();
let (_root4, s_exact) = four.sqrt(NE);
assert!(!s_exact.inexact());
assert!(s_exact.is_ok());
```

`is_ok()` is true exactly when no flag is set. This is how you tell an exact
answer from a rounded one without inspecting digits: the same `sqrt` call that
rounded two reports a clean status for four. The flag is the witness, computed by
the kernel, not reconstructed by you.

## Step 3: force an overflow

pfloat carries the exponent in an `i64`, so `OVERFLOW` is hard to reach by
accident: the exponent range dwarfs any normal computation. You reach it on
purpose by pushing a value past `i64::MAX`. `scale_by_pow2(k)` multiplies by two
to the `k`, adjusting only the exponent; scaling a value whose exponent is
already positive by `i64::MAX` saturates and raises the flag.

```rust
let three = BigFloat::try_from_i64_exact(3, 53).unwrap();
let (saturated, s_over) = three.scale_by_pow2(i64::MAX);
assert!(s_over.overflow());
assert!(!s_over.underflow());
assert!(!saturated.is_infinite());
```

Two details matter. The result is finite, not infinity: pfloat has no maximum
exponent to clamp against, so it saturates to the largest finite value it can
represent and tells you with the flag. And `overflow()` and `underflow()` are
distinct predicates; scaling the other direction by `i64::MIN` would set
`UNDERFLOW` instead. You act on the one you asked about.

## Step 4: accumulate flags in the thread local

The per call return is the only flag transport under `no_std`. Under `std`,
every flag producing operation also unions its status into a per thread sticky
cell, the IEEE 754 "global flags" mental model: run a whole computation, then ask
once whether anything inexact or invalid happened along the way.

The accessors live in the `pfloat::flags` module. Reading `src/status.rs`
confirms the names: `clear()` resets the cell and hands back its previous value,
`test()` reads without clearing, `set()` overwrites, and `raise()` unions a
status in (the operation the kernels use internally, by bitwise OR). This guide
uses `clear` and `test`.

```rust
let _ = flags::clear();

let _ = two.sqrt(NE);  // raises INEXACT
let _ = four.sqrt(NE); // raises nothing
let accumulated = flags::test();
assert!(accumulated.inexact());
```

You cleared the cell, ran an inexact op and an exact one, then read the union of
both. `INEXACT` survived because the accumulator unions; the exact second call
did not unset it. To start a fresh measurement, clear again and read back the
cleared set.

```rust
let cleared = flags::clear();
assert!(cleared.inexact());
let after = flags::test();
assert!(after.is_ok());
```

`clear()` returns what it cleared, so you lose nothing by resetting: inspect the
returned status, and the cell is already clean for the next computation. Each
thread owns its own cell, so a spawned worker starts at `OK` and never sees
another thread's flags.

## What you learned

Three ideas, building on tutorial 1's `(value, Status)` pair.

- The status surface is six predicates: `inexact()`, `overflow()`,
  `underflow()`, `invalid()`, `div_by_zero()`, and `is_ok()`. An exact result
  sets none, so `is_ok()` is true; a rounded one sets `INEXACT`.
- `OVERFLOW` is reachable but rare, because the exponent is an `i64`. When it
  fires, the result saturates to a finite value rather than infinity, and the
  flag is how you know.
- Under `std`, flags also accumulate in a per thread sticky cell. `flags::clear`
  resets it and returns the old value; `flags::test` reads it. Use the
  accumulator for the IEEE "global flags" model, the per call return when a
  function must not observe its caller's flags.

## Next

- Tutorial 3 walks the five rounding modes and shows when a directed mode, not
  `NearestEven`, is the right tool, including how the flags differ under
  directed rounding.
- Tutorial 4 moves from a single rounded value with a flag to a rigorous
  enclosure with `pfloat-ball`, when you need a guaranteed bound rather than a
  flag on a point.