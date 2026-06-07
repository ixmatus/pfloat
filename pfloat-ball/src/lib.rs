//! pfloat-ball: rigorous arbitrary-precision real ball arithmetic.
//!
//! A *ball* is a midpoint-radius enclosure `[m ± r]` denoting the closed
//! real interval `[m − r, m + r]`. pfloat-ball carries the midpoint as a
//! full-precision pfloat scalar ([`pfloat::BigFloat`] or
//! [`pfloat::FixedFloat`]) and the radius as a [`Mag`], a small
//! upward-rounded unsigned magnitude. Every ball operation computes the
//! midpoint with pfloat's correctly-rounded kernel and bounds the radius
//! by the rounding error those kernels already compute, so the result is
//! a *sound* enclosure: the true mathematical result of the operation,
//! applied to any point in the input ball, lies inside the output ball.
//!
//! This crate is the first shippable cut of pfloat's Phase 4 rigorous
//! enclosure tower. The v1.0 surface is real-ball arithmetic and the
//! elementary functions over [`pfloat::BigFloat`]; special functions and
//! the complex / IEEE-1788 faces are separate later work.
//!
//! # Soundness is a type fact where it can be
//!
//! [`Mag`] makes an unsound radius unrepresentable: it has no sign (a
//! negative radius cannot be written) and no NaN, and every `Mag`
//! operation rounds toward `+∞` by the type's contract (an
//! inward-rounded radius cannot be written). The remaining soundness
//! obligations live in the in-tree enclosure spec and the verification
//! lanes, not in reviewer memory.
//!
//! # Features
//!
//! - `big` (default): `Ball<BigFloat>`, the headline dynamic-precision
//!   type. Pulls `alloc`.
//! - `fixed`: `Ball<FixedFloat<PREC>>`, compile-time precision.
//! - `std` (default): pfloat's thread-local sticky flags and std error
//!   impls.
//! - `exp-log`, `trig`: enable the matching ball elementary functions.
//! - `serde`: `Serialize`/`Deserialize` for [`Mag`] and `Ball`.
//!
//! A bare `--no-default-features` build exposes only [`Mag`] (alloc-free,
//! the minimal embedded surface and the Kani target).

#![feature(generic_const_exprs)]
#![allow(incomplete_features)]
#![cfg_attr(not(feature = "std"), no_std)]
#![forbid(unsafe_code)]

#[cfg(feature = "alloc")]
extern crate alloc;

#[cfg(feature = "big")]
pub mod accuracy;
#[cfg(feature = "big")]
pub mod arith;
#[cfg(feature = "big")]
pub mod ball;
#[cfg(feature = "exp-log")]
pub mod elem;
#[cfg(feature = "big")]
pub mod io;
#[cfg(kani)]
mod kani_harness;
pub mod mag;
#[cfg(feature = "big")]
pub mod scalar;
#[cfg(feature = "big")]
pub mod spec;

#[cfg(feature = "big")]
pub use accuracy::refine_to_accuracy;
#[cfg(feature = "big")]
pub use ball::{Ball, BallError};
#[cfg(feature = "big")]
pub use io::BallParseError;
pub use mag::Mag;
#[cfg(feature = "big")]
pub use scalar::RealScalar;
