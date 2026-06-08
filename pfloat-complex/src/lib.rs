//! pfloat-complex: componentwise correctly-rounded complex arithmetic.
//!
//! A [`Complex<T>`] is a pair `(re, im)` of pfloat scalars. Each operation
//! is correctly rounded *componentwise*: the real and imaginary parts are
//! each rounded under their own real rounding mode, which is the only
//! coherent strong rounding claim for complex numbers (they carry no total
//! order, so a single complex directed rounding has no meaning). This is
//! the model MPC uses, and the honest lift of pfloat's scalar contract.
//!
//! This crate is the MPC analog in the pfloat family: arbitrary-precision
//! complex arithmetic in pure Rust, where `num-complex` is a bare container
//! and `rug`/MPC require a C toolchain. It depends only on `pfloat`.
//!
//! # The sealed scalar trait
//!
//! `Complex<T>` is generic over `T:` [`RealScalar`], a *sealed* trait
//! implemented only for [`pfloat::BigFloat`] and
//! [`pfloat::FixedFloat<PREC>`]. Sealing makes "each component is a
//! correctly-rounded pfloat scalar" a fact this crate's own surface cannot
//! be made to break: a `Complex` here can never be instantiated over an
//! unverified or wrongly-rounded scalar type. The seal is scoped, not
//! universal: because Phase 3 shipped `num_traits::Num` for
//! `FixedFloat<PREC>` (pfloat ADR-0070), a third party can still build a
//! `num_complex::Complex<FixedFloat<P>>` outside this crate. `RealScalar`
//! closes *pfloat-complex's* inhabitant set, not the universe of generic
//! numeric code.
//!
//! The trait is defined here, independently of `pfloat-ball`'s `RealScalar`:
//! complex does not depend on the ball (the two are sibling Phase 4 stars
//! with no edge between them). A shared trait, if one is ever extracted,
//! comes after the concrete crates exist, never before (the roadmap rule).
//!
//! # Status
//!
//! Version `1.0`: the public API is settled, and the aim is to keep `1.x`
//! changes additive (semver); as a personal project this is an intent rather
//! than a contractual guarantee (ADR-0093). This cut ships the `Complex<T>`
//! type, componentwise correctly-rounded arithmetic (`add`, `sub`, `neg`,
//! `conj`, `norm_sqr`, `mul`, `div` with the C99 Annex G §G.5.1 infinity
//! recovery), magnitude and phase (`abs`, `arg`, `to_polar`), and the
//! elementary core (`sqrt`, `exp`, `log`) with their C99 Annex G branch cuts.
//! The trigonometric, hyperbolic, and inverse functions, `pow`/`cis`/
//! `from_polar`, and the `ComplexBall` join are additive later work (ADR-0093).
//! Part of the pfloat workspace; built and verified alongside it.

#![feature(generic_const_exprs)]
#![allow(incomplete_features)]
#![cfg_attr(not(feature = "std"), no_std)]
#![forbid(unsafe_code)]

#[cfg(feature = "alloc")]
extern crate alloc;

// Advisory Kani proof harnesses over the componentwise `Status` merge (the one
// Vec-free invariant; the BigFloat kernels are CBMC-hostile, ADR-0062/0092).
// Compiled only under `cfg(kani)`, which `cargo kani` sets.
#[cfg(kani)]
mod kani_harness;

#[cfg(feature = "trig")]
mod cexp;
#[cfg(feature = "trig")]
mod clog;
#[cfg(feature = "big")]
pub mod complex;
#[cfg(feature = "exp-log")]
mod csqrt;
#[cfg(feature = "big")]
mod div;
#[cfg(feature = "exp-log")]
mod enclosure;
#[cfg(feature = "big")]
pub mod scalar;
#[cfg(feature = "big")]
mod specials;

#[cfg(feature = "big")]
pub use complex::Complex;
#[cfg(feature = "big")]
pub use scalar::RealScalar;
