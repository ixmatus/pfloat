# How to: run pfloat in a `no_std` or embedded build

Type: how to (task oriented). Time: ten minutes. Runnable companion:
`examples/fixed_precision.rs` in the `pfloat` crate.

This guide swaps the heap backed `BigFloat` for `FixedFloat<PREC>`, the type
that carries its precision at compile time and stacks its mantissa in a fixed
array. You will build a binary128 width value, take a square root, and learn the
two feature flags an embedded target needs, plus the one caveat that still
reaches for `alloc`.

## Goal

Compute correctly rounded square roots at a width fixed by the type rather than
chosen at runtime, in a build that drops `std` and the dynamic precision
`BigFloat`. `FixedFloat<113>` is binary128: the IEEE 754 quad mantissa is 113
bits, and the storage size is known to the compiler, so the value lives on the
stack.

## 1. Set the cargo line for `no_std`

The default features pull in `std`, the thread local exception accumulator, and
the heap backed `BigFloat`. An embedded target turns all of that off and asks
only for `fixed`.

```toml
[dependencies]
# Not yet on crates.io; depend from git. Needs a nightly toolchain.
pfloat = { git = "https://github.com/ixmatus/pfloat", tag = "v1.0.0", default-features = false, features = ["fixed"] }
```

`default-features = false` drops `std`, `fmt`, and `big` from the build.
`fixed` brings in `FixedFloat<const PREC: u32>`, the const generic float. The
type compiles without `std`; the Status it returns is the same per call flag
record, minus the thread local accumulator that needs `std`.

## 2. Name the width as a type parameter

`FixedFloat<PREC>` puts the precision in the type, so a type alias reads like a
named numeric format. The const generic bound rides on the type, so every crate
that spells `FixedFloat` declares the same nightly feature the library declares;
the pinned toolchain provides it.

```rust
#![feature(generic_const_exprs)]
#![allow(incomplete_features)]

use pfloat::{FixedFloat, RoundingMode};

// binary128 width: the IEEE 754 quad mantissa is 113 bits.
type Quad = FixedFloat<113>;

const NE: RoundingMode = RoundingMode::NearestEven;
```

## 3. Build a value with no precision argument

This is the visible difference from `BigFloat`. `BigFloat::try_from_i64_exact`
takes a runtime bit count; `FixedFloat::try_from_i64_exact` takes none, because
the width is `113`, baked into the type. The constructor stays fallible:
not every integer fits exactly in every width.

```rust
let two = Quad::try_from_i64_exact(2).unwrap();
assert_eq!(Quad::PRECISION, 113);
```

`Quad::PRECISION` reads back the compile time width as a `u32`. There is no
runtime precision field to inspect; the type is the precision.

## 4. Compute, and read the same Status

The arithmetic and root surface mirrors `BigFloat` exactly: a rounding mode in,
a `(value, Status)` pair out. The flags carry the same meaning per call.

```rust
let (root, status) = two.sqrt(NE);
assert!(status.inexact());

// A perfect square is exact at this width: sqrt(4) = 2, no rounding.
let four = Quad::try_from_i64_exact(4).unwrap();
let (_two_again, exact_status) = four.sqrt(NE);
assert!(!exact_status.inexact());
```

`sqrt(2)` is irrational, so it rounds and sets `INEXACT` at any finite width.
`sqrt(4)` lands on `2` exactly, so the flag stays clear. `inexact()` is the per
call exactness witness whether the width is fixed or dynamic.

## 5. Widen to `BigFloat` to print or to call transcendentals

`FixedFloat` carries arithmetic, comparison, roots, and `fma` with no heap. The
decimal string formatter and the elementary functions live on `BigFloat`, so you
widen with `to_big()` when you reach for them.

```rust
let root_big = root.to_big();
assert_eq!(root_big.precision(), 113);

let digits = root_big.to_decimal_string(30, NE);
println!("sqrt(2) at binary128 width = {digits}");
assert!(digits.starts_with("1.4142135623730950488016887242"));
```

`to_big()` is the same widening as `From<FixedFloat<PREC>> for BigFloat`. The
result is a `BigFloat` at precision `113`, holding the identical value.

## The alloc caveat

`FixedFloat` is alloc free for arithmetic and roots, but the transcendentals are
not, and the reason is structural, not an oversight. Correct rounding of `exp`,
`ln`, `sin`, and the rest grows the working precision on demand: when the true
result sits too close to a rounding boundary, the kernel recomputes wider until
the rounding is decided (the Ziv strategy). That growth can pass any fixed width,
so a transcendental cannot promise to stay inside a `[u64; N]` array. The
elementary surface therefore lives on `BigFloat` behind `alloc`, and a fully
heapless target gets arithmetic, comparison, square root, cube root, and `fma`,
but not `exp` or `sin`.

Two consequences follow. First, an embedded build that needs only the algebraic
operations stays on `FixedFloat` with no allocator. Second, a build that needs
transcendentals enables `alloc` (and `exp-log` or `trig`) and accepts a heap,
even if it never enables `std`; `alloc` over `std` is the smaller step.

## What you achieved

- A `no_std` build line: `default-features = false` plus `fixed`, dropping
  `std`, `fmt`, and `big`.
- Compile time precision via `FixedFloat<113>`: the width is the type, the
  constructor takes no precision argument, and the mantissa stacks in a fixed
  array.
- The alloc boundary: arithmetic and roots are heap free; correct rounding of
  transcendentals grows precision past any fixed width, so they stay on
  `BigFloat` behind `alloc`.

## Next

- Guide 01 builds the same square root with the dynamic precision `BigFloat`,
  where the precision is a runtime argument.
- Guide 02 reads the full Status surface, including the thread local accumulator
  that `no_std` builds leave behind.
- Guide 03 walks the five rounding modes that the `mode` argument selects, here
  and across the whole family.