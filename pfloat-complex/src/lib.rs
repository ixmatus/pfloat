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
//! # Status (pre-1.0)
//!
//! Version `0.x`; the API will change until `1.0`. This cut ships the
//! `Complex<T>` type and componentwise additive arithmetic (`add`, `sub`,
//! `neg`, `conj`). Multiplication and division (the cancellation-safe
//! `mul_sub_mul` / `mul_add_mul` forms), magnitude and phase, and the
//! elementary functions with their C99 Annex G branch cuts are later
//! slices. Part of the pfloat workspace; built and verified alongside it.

#![feature(generic_const_exprs)]
#![allow(incomplete_features)]
#![cfg_attr(not(feature = "std"), no_std)]
#![forbid(unsafe_code)]

#[cfg(feature = "alloc")]
extern crate alloc;

#[cfg(feature = "big")]
pub mod complex;
#[cfg(feature = "big")]
mod div;
#[cfg(feature = "big")]
pub mod scalar;
#[cfg(feature = "big")]
mod specials;

#[cfg(feature = "big")]
pub use complex::Complex;
#[cfg(feature = "big")]
pub use scalar::RealScalar;
