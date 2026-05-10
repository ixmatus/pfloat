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
//! product BigFloat then re-rounds via `add_round` for the IEEE
//! 754 §9.4 single-rounding guarantee. The slice also fixes a
//! latent `addsub` bug exposed by FMA (the result-exponent formula
//! used `e_s - p_s + 1` instead of the genuine
//! `min(scale_l, scale_s)`; both happened to coincide whenever
//! `scale_s` was the minimum, which all of slice 1c's tests
//! exercised, but cross-precision FMA does not).

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
mod mantissa;
#[cfg(feature = "big")]
mod ops;
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
#[cfg(feature = "big")]
pub use rounding::RoundingMode;
