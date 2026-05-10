//! pfloat: pure Rust correctly-rounded arbitrary-precision floats.
//!
//! This crate is pre-1.0. The public surface is unstable and the
//! arithmetic kernels are not yet implemented. See `DESIGN.md` at
//! the repository root for the full design and `docs/decisions/` for
//! the architecture decision records that capture the load-bearing
//! choices.
//!
//! # Scope target
//!
//! v1.0 ships an MPFR-equivalent surface: arithmetic with all five
//! IEEE 754-2019 rounding modes, sticky exception flags, correctly-
//! rounded transcendentals (exp, log, trig, hyperbolic, pow), and
//! special functions (gamma family, erf family, Bessel, zeta,
//! Ei/Si/Ci, Airy, AGM).
//!
//! Two precision profiles share the same operations:
//!
//! - [`BigFloat`]: runtime-determined precision. Requires the
//!   `alloc` feature.
//! - `FixedFloat<const PREC: u32>`: compile-time precision via a
//!   const generic. Stack-allocated, works without `alloc`. Lands
//!   in slice 1g.
//!
//! # Slices 1a–1c (currently shipped)
//!
//! 1a: [`BigFloat`] type, classification predicates, comparison
//! ([`partial_cmp`](BigFloat::partial_cmp),
//! [`total_cmp`](BigFloat::total_cmp), [`min`](BigFloat::min),
//! [`max`](BigFloat::max)), and exact integer construction
//! ([`try_from_i64_exact`](BigFloat::try_from_i64_exact)).
//!
//! 1b: full [`Status`] (all five IEEE flags), [`RoundingMode`], the
//! universal rounding pipeline, std-only thread-local flag accessors
//! ([`flags`](status::flags)), and rounding-required constructors
//! ([`try_from_i64_round`](BigFloat::try_from_i64_round) and
//! [`round_to_precision`](BigFloat::round_to_precision)).
//!
//! 1c: first arithmetic kernel.
//! [`add`](BigFloat::add) and [`sub`](BigFloat::sub) (plus
//! [`add_round`](BigFloat::add_round) / [`sub_round`](BigFloat::sub_round)
//! and `_with_flags` siblings) handle NaN propagation, signed
//! infinities, signed-zero arithmetic (including the IEEE sign rule
//! for `±0 ± ±0` under `TowardNegative`), and mantissa alignment by
//! `2^Δ` for arbitrary exponent gaps.
//!
//! 1d: [`mul`](BigFloat::mul) (plus
//! [`mul_round`](BigFloat::mul_round) and `_with_flags` siblings).
//! Schoolbook + Karatsuba multiplication via the shared
//! `ops::limbs` module; FFT (Schönhage-Strassen) deferred to 1.x
//! per ADR-0010. Handles `0 × ∞ → qNaN + INVALID` per IEEE 754
//! §7.2.
//!
//! 1e: [`div`](BigFloat::div) (plus
//! [`div_round`](BigFloat::div_round) and `_with_flags` siblings).
//! Bit-by-bit long division of the mantissas with sticky-bit
//! tracking from the remainder; routes through the rounding
//! pipeline. Raises `DIV_BY_ZERO` per IEEE 754 §7.3 for
//! `finite_nonzero / 0` and `INVALID` for `0 / 0` and `∞ / ∞`.
//!
//! 1f: [`sqrt`](BigFloat::sqrt) and
//! [`fma`](BigFloat::fma), the last two arithmetic primitives in
//! Phase 1 (plus `_round` and `_with_flags` siblings each). `sqrt`
//! uses bit-by-bit integer square root with parity-adjusted shift
//! so the result exponent splits cleanly; `fma` builds an exact
//! product `BigFloat` then re-rounds via `add_round` for the IEEE
//! 754 §9.4 single-rounding guarantee. The slice also fixes a
//! latent `addsub` bug exposed by FMA (the result-exponent formula
//! used `e_s - p_s + 1` instead of the genuine
//! `min(scale_l, scale_s)`; both happened to coincide whenever
//! `scale_s` was the minimum, which all of slice 1c's tests
//! exercised, but cross-precision FMA does not).
//!
//! 1g: [`FixedFloat<PREC>`](FixedFloat) (const-generic counterpart
//! to [`BigFloat`] with stack-allocated mantissa) plus `core::ops`
//! operator overloads (`Add`, `Sub`, `Mul`, `Div`, `Neg`, and the
//! `*Assign` siblings) for both types behind the `ops` feature.
//! Phase 1 arithmetic surface complete on both types.
//!
//! # Phase 2 (in progress)
//!
//! 2a: Decimal string parsing via
//! [`BigFloat::parse_str`](BigFloat::parse_str) and
//! [`FixedFloat::parse_str`](FixedFloat::parse_str). Accepts
//! signed decimal numbers with optional decimal point and `e/E`
//! exponent, plus case-insensitive `nan` / `inf` / `infinity`.
//! Routes through the rounding pipeline via a multi-precision
//! `m × 5^exp × 2^exp` decomposition.
//!
//! 2b: Decimal formatting via
//! [`BigFloat::to_decimal_string`](BigFloat::to_decimal_string)
//! and the standard [`Display`](core::fmt::Display) trait for both
//! types. `Display` uses
//! [`BigFloat::round_trip_digit_count`](BigFloat::round_trip_digit_count)
//! digits — enough that `parse_str` at the same precision recovers
//! the exact value. Output uses fixed-point notation for moderate
//! magnitudes and scientific notation otherwise. Shortest
//! round-trip (Dragon4 / Steele-White) can be a future
//! optimization; the current output is correct, just not always
//! minimal.
//!
//! # Phase 3 (in progress)
//!
//! 3a: First transcendental — [`BigFloat::exp`](BigFloat::exp) (and
//! `FixedFloat<PREC>::exp`) behind the `exp-log` cluster feature.
//! Algorithm: range-reduce by `ln(2)` (hardcoded 1024-bit constant)
//! so the residual `|r| ≤ ln(2)/2`, evaluate the Taylor series, and
//! compose with a free exponent shift for the `2^k` factor. Fixed
//! 64-bit guard above target precision; full Ziv-strategy retry
//! deferred. Lefèvre–Muller worst-case verification wires in during
//! Phase 6.
//!
//! 3b: [`BigFloat::ln`](BigFloat::ln) (and `FixedFloat<PREC>::ln`).
//! Range-reduce by the binary exponent
//! (`x = m × 2^e` with `m ∈ [1, 2)`), then `ln(x) = ln(m) + e · ln(2)`.
//! The `ln(m)` part uses the atanh series
//! `ln(m) = 2·atanh((m-1)/(m+1))`, which converges roughly three
//! bits per term over `m ∈ [1, 2)`. `ln(2)` is shared with `exp`.
//! Same fixed 64-bit guard as 3a.
//!
//! 3c: [`BigFloat::pow`](BigFloat::pow) (and `FixedFloat<PREC>::pow`).
//! For positive finite `x` and finite `y`, evaluates
//! `exp(y · ln(x))` at working precision and rounds back to target.
//! Negative bases dispatch on integer parity of `y` (per IEEE 754-2019
//! §9.2.1): odd integer flips the sign, even integer leaves it
//! positive, non-integer raises `INVALID` and returns qNaN. Full
//! special-case table for `0`, `±∞`, `NaN`, and the `pow(±1, ±∞) = 1`
//! rule.
//!
//! 3d: thin wrappers around `exp` and `ln`.
//! [`expm1`](BigFloat::expm1) and [`log1p`](BigFloat::log1p) handle
//! the small-argument cancellation regime by boosting working
//! precision proportional to `-exponent(x)` before invoking the base
//! transcendental, then rounding back to target.
//! [`exp2`](BigFloat::exp2) / [`exp10`](BigFloat::exp10) compose as
//! `exp(x · ln(b))` for `b ∈ {2, 10}`; `ln(2)` reuses the hardcoded
//! 1024-bit constant, `ln(10)` is computed lazily per call.
//! [`log2`](BigFloat::log2) / [`log10`](BigFloat::log10) compose as
//! `ln(x) / ln(b)`. All six are mirrored on `FixedFloat`.
//!
//! 3f: forward trig family.
//! [`sin`](BigFloat::sin), [`cos`](BigFloat::cos), and
//! [`tan`](BigFloat::tan) on `BigFloat` and `FixedFloat`. Behind a
//! separate `trig` feature flag. Argument reduction multiplies `x`
//! by a hardcoded 4096-bit table of `2/π` (Payne-Hanek style),
//! yielding `q mod 4` (the quadrant index) and a reduced
//! `r ∈ [−π/4, π/4]`. The kernel then evaluates `sin(r)` or
//! `cos(r)` via Taylor and dispatches on the quadrant. The table
//! caps the supported input range at roughly `|x| < 2^3000`; past
//! that the kernel raises `INVALID` and returns qNaN. `tan` is
//! defined as `sin/cos` over the same reduction.
//!
//! 3g: inverse trig family.
//! [`atan`](BigFloat::atan) uses the identity
//! `atan(x) = π/2 − atan(1/x)` for `|x| > 1` to bring the argument
//! into `[0, 1]`, then applies the half-angle identity
//! `atan(y) = 2 · atan(y / (1 + sqrt(1 + y²)))` repeatedly until
//! `|y| < 1/16`, then sums the Taylor series. `asin(x)` and
//! `acos(x)` route through `atan` via cancellation-free
//! identities: `asin(x) = 2 · atan(x / (1 + sqrt(1 − x²)))`;
//! `acos(x) = 2 · atan(sqrt((1 − x)/(1 + x)))` for `x ≥ 0`, and
//! the reflected variant for `x < 0`. [`atan2`](BigFloat::atan2)
//! dispatches the full IEEE 754-2019 §9.2.1 special-case table for
//! `(y, x)` signs at zero, infinity, and the four quadrants, then
//! reduces to `atan(y/x)` plus a `±π` shift for the second and
//! third quadrants. All four are mirrored on `FixedFloat`.
//!
//! 3e: hyperbolic family.
//! [`sinh`](BigFloat::sinh) uses
//! `(expm1(x) − expm1(−x))/2` to avoid the cancellation of
//! `exp(x) − exp(−x)` for small `|x|`.
//! [`cosh`](BigFloat::cosh) is the direct
//! `(exp(x) + exp(−x))/2` (both summands non-negative, no
//! cancellation).
//! [`tanh`](BigFloat::tanh) computes
//! `(1 − exp(−2|x|))/(1 + exp(−2|x|))` then flips the sign for
//! negative `x`; the form is well-behaved at `±∞` (returns ±1
//! directly).
//! [`asinh`](BigFloat::asinh) uses
//! `sign(x) · log1p(|x| + |x|²/(sqrt(|x|² + 1) + 1))`, avoiding the
//! catastrophic cancellation in `ln(x + sqrt(x² + 1))` for large
//! negative `x`.
//! [`acosh`](BigFloat::acosh) uses
//! `log1p((x − 1) + sqrt((x − 1)(x + 1)))` so the near-1 case
//! routes through `log1p`'s cancellation-aware path.
//! [`atanh`](BigFloat::atanh) uses
//! `(log1p(x) − log1p(−x))/2`. Domain check rejects `|x| > 1` with
//! `INVALID` and `|x| = 1` with `DIV_BY_ZERO`. All six are mirrored
//! on `FixedFloat`.

// pfloat depends on `feature(generic_const_exprs)` for the
// `FixedFloat<const PREC: u32>` storage spelling that lands in
// slice 1g. ADR-0011 records the trade-off (nightly toolchain
// required); the feature is `incomplete` upstream, so its lint is
// allowed at the crate root.
#![feature(generic_const_exprs)]
#![allow(incomplete_features)]
#![cfg_attr(not(feature = "std"), no_std)]
#![forbid(unsafe_code)]

#[cfg(feature = "alloc")]
extern crate alloc;

#[cfg(feature = "std")]
extern crate std;

#[cfg(feature = "big")]
mod big;
#[cfg(feature = "big")]
mod class;
#[cfg(feature = "big")]
mod classify;
#[cfg(feature = "big")]
mod cmp;
#[cfg(all(feature = "big", feature = "fixed"))]
mod fixed;
#[cfg(feature = "big")]
mod fmt;
mod mantissa;
#[cfg(feature = "exp-log")]
mod math;
#[cfg(feature = "big")]
mod ops;
#[cfg(all(feature = "big", feature = "ops"))]
mod ops_traits;
#[cfg(feature = "big")]
mod parse;
#[cfg(feature = "big")]
mod rounding;
mod sign;
mod status;

pub use sign::Sign;
pub use status::Status;

#[cfg(feature = "std")]
pub use status::flags;

#[cfg(feature = "big")]
pub use big::{BigFloat, BuildError};
#[cfg(feature = "big")]
pub use classify::IeeeClass;
#[cfg(all(feature = "big", feature = "fixed"))]
pub use fixed::{ClassFixed, FixedFloat};
#[cfg(feature = "big")]
pub use parse::ParseError;
#[cfg(feature = "big")]
pub use rounding::RoundingMode;
