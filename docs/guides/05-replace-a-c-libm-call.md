# How to: replace a C libm call with a correctly rounded pure Rust one
Type: how to (task oriented). Time: ten minutes. Runnable companion: pfloat-libm/examples/correctly_rounded_exp.rs.

You have code that calls `f32::exp`, `f64::sin`, or another hardware float transcendental, and you want the last bit to be correct on every input. This guide swaps those calls for pfloat-libm, a pure Rust `no_std` plus `alloc` shell that computes each value at high working precision and rounds to the hardware width only once an enclosure proves the rounding. No C toolchain, no `libm` link step, no platform variance in the low bits.

## Goal

Turn a fast but last bit unreliable `std` float call into a correctly rounded one, and reach for the directed `_round` form when you need the enclosing bracket and the IEEE 754 status flags.

## 1. Add the dependency

pfloat-libm pulls in pfloat as its kernel. The default features cover the elementary surface; the example below needs `exp`, which lives behind the `exp-log` kernel feature the default build already enables. The crates are not yet published to crates.io, so depend on pfloat-libm from git, and note that pfloat needs a nightly toolchain for `feature(generic_const_exprs)`.

```toml
[dependencies]
pfloat-libm = { git = "https://github.com/ixmatus/pfloat" }
```

There is nothing else to install. A C `libm` binding needs a C compiler, a linked system library, and whatever rounding the platform `libm` happened to ship; pfloat-libm needs only the Rust nightly toolchain, and the rounding is the same on every target.

## 2. Swap the bare call

The bare `std` method becomes a free function in the `f32` or `f64` module. The signature is the same shape: one float in, one float out. The difference is the guarantee, not the call site.

```rust
use pfloat_libm::f32 as lm;

let x: f32 = 1.5;

// Before: std, fast, last bit not guaranteed.
let _std_value = x.exp();

// After: correctly rounded to nearest even.
let value = lm::exp(x);
println!("correctly rounded exp(1.5) = {value}");
```

`lm::exp` rounds to nearest even, the IEEE 754 default, and the returned `f32` is the grid point closest to the true `exp(1.5)`. The same module exposes `ln`, `log2`, `log10`, `sin`, `cos`, `tan`, `sqrt`, `cbrt`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, and the rest of the surface; the `f64` module mirrors it for double width.

## 3. Reach for the directed form when you need the bracket

The `_round` form takes an explicit `RoundingMode` and returns a `(value, Status)` pair. Rounding the same true result toward `+inf` and toward `-inf` gives you a two sided bracket: the true value sits between `down` and `up`. For a transcendental argument the two roundings land on adjacent grid points, one ULP apart, and the status reports the result is inexact.

```rust
use pfloat_libm::{f32 as lm, RoundingMode};

let x: f32 = 1.5;
let (up, status) = lm::exp_round(x, RoundingMode::TowardPositive);
let (down, _) = lm::exp_round(x, RoundingMode::TowardNegative);

assert!(up > down);
assert_eq!(up.to_bits() - down.to_bits(), 1); // one ULP apart
assert!(status.inexact());
```

This is the move when an error analysis needs a guaranteed enclosure of the true result rather than a single nearest value. `TowardZero`, `TowardPositive`, and `TowardNegative` give directed rounding; `NearestEven` and `NearestAway` give the two tie breaking nearest modes.

## 4. Read the status flags

The `Status` returned alongside the directed value carries the sticky IEEE 754 flags. Each predicate is a `const fn` returning `bool`: `inexact()`, `overflow()`, `underflow()`, `invalid()`, `div_by_zero()`, and `is_ok()`. An argument far past the `f32` overflow point saturates to infinity and raises both overflow and inexact through a fast path that skips argument reduction.

```rust
use pfloat_libm::{f32 as lm, RoundingMode};

let (big, st) = lm::exp_round(1000.0, RoundingMode::NearestEven);
assert!(big.is_infinite() && st.overflow() && st.inexact());
```

## When correct rounding matters

A `std` or C `libm` transcendental is usually within one ULP of the true value, which is fine for graphics, physics, and most numerics. Correct rounding earns its keep when the last bit is load bearing:

- Reproducibility across platforms. The system `libm` rounds differently on different operating systems and architectures, so the same program prints different low bits. pfloat-libm returns the same correctly rounded value everywhere, which makes float output a stable artifact you can hash, snapshot, or commit.
- A correct rounding contract a downstream consumer depends on. When the value feeds a decimal formatter, a comparison against a stored result, or a second rounding step, a one ULP wobble can flip the visible output. Correct rounding removes that source of drift.
- Directed rounding for verified error bounds. The `TowardPositive` and `TowardNegative` pair gives a true enclosure of the result, which a `std` nearest only call cannot provide.

The cost is speed: pfloat-libm widens to arbitrary precision and runs a Ziv loop, so it allocates and is far slower than a hardware `libm` call. Use it where the bit matters and keep `std` where throughput does.

## Next

- Guide 01 introduces the underlying `BigFloat` value, `Status`, and `to_decimal_string` that pfloat-libm rounds through.
- Guide 02 covers the thread local flag accumulator for the IEEE 754 global flag mental model.
- Guide 04 builds enclosures with pfloat-ball when you want a tracked error radius rather than a single directed bracket.