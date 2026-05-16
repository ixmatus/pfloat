//! Kani verification harnesses for pfloat.
//!
//! Compiled only under `cfg(kani)`. The parent declaration in
//! `src/lib.rs` gates the entire module; individual harnesses use
//! `#[kani::proof]` without an inner `cfg` gate.
//!
//! Slice 6a lands the scaffold and a canonical example of each
//! harness shape for the `add` operation. Slices 6b–6f add coverage
//! across the rest of the surface. ADR-0012 records the architecture
//! and the load-bearing decision: the CI Kani lane is advisory
//! (`continue-on-error: true`), not blocking, per the
//! `feedback_kani_ci_timeout_ok.md` engineering memory.
//!
//! ## Running locally
//!
//! ```sh
//! cargo kani
//! ```
//!
//! Kani sets `cfg(kani)` automatically when it invokes the compiler;
//! the `kani = []` feature in `Cargo.toml` is a placeholder that
//! exists so external tooling can opt in symmetrically without
//! tripping `unexpected_cfgs`.
//!
//! ## Operand bounding
//!
//! pfloat does not have ±MAX / ±MIN_POSITIVE constants — its
//! exponent is `i64` and precision is arbitrary. The
//! [`helpers::nondet_constant_at`] selector returns one of eight
//! canonical values (qNaN, sNaN, ±∞, ±0, ±1) at a fixed precision.
//! Bounded-normal generators with parameterized exponent ranges
//! land in slice 6b alongside the rounding-direction harnesses.

pub(super) mod helpers;

mod acos;
mod acosh;
mod add;
#[cfg(feature = "agm")]
mod agm;
mod asin;
mod asinh;
mod atan;
mod atan2;
mod atanh;
mod beta;
mod classify;
mod cmp;
mod cos;
mod cosh;
mod digamma;
mod div;
mod erf;
mod erfc;
mod exp;
mod exp10;
mod exp2;
mod expm1;
mod fma;
mod fmt;
mod gamma;
mod lgamma;
mod ln;
mod log10;
mod log1p;
mod log2;
mod mul;
mod parse;
mod pow;
mod sin;
mod sinh;
mod sqrt;
mod sub;
mod tan;
mod tanh;

#[cfg(feature = "integrals")]
mod ci;
#[cfg(feature = "integrals")]
mod ei;
#[cfg(feature = "integrals")]
mod li;
#[cfg(feature = "integrals")]
mod si;
