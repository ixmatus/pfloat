//! pfloat: pure Rust correctly-rounded arbitrary-precision floats.
//!
//! This crate is pre-1.0. The public surface is unstable and the
//! algorithmic kernels are not yet implemented. See `DESIGN.md` at
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
//! - [`FixedFloat`]: compile-time precision via a const generic.
//!   Stack-allocated, works without `alloc`.
//!
//! Both types are stubs at this stage of the project.

#![cfg_attr(not(feature = "std"), no_std)]
#![forbid(unsafe_code)]

#[cfg(feature = "alloc")]
extern crate alloc;

#[cfg(feature = "std")]
extern crate std;

// `BigFloat` and `FixedFloat` will land as separate modules in Phase 1.
// The names are reserved here so the crate surface is visible from
// the start; the types do not exist yet.
#[cfg(feature = "big")]
#[doc(hidden)]
pub struct BigFloat;

#[cfg(feature = "fixed")]
#[doc(hidden)]
pub struct FixedFloat<const PREC: u32>;
