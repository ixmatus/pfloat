//! pfloat-libm: a pure-Rust, correctly-rounded `libm`.
//!
//! pfloat-libm is a thin shell over [pfloat](https://github.com/ixmatus/pfloat).
//! A call widens the hardware float to an arbitrary-precision `BigFloat`,
//! evaluates the function correctly-rounded through pfloat's kernel, and
//! rounds the result back to the hardware width under an outer Ziv loop that
//! commits a value only once an enclosure proves it. The outer loop is what
//! makes the second rounding safe: pfloat-libm never rounds an intermediate
//! blindly, so no double-rounding error survives.
//!
//! The unary `f32` surface is designed to be exhaustively verified: every one
//! of the 2^32 `binary32` inputs is checked against an independent oracle. The
//! `f64` surface rests on differential testing plus worst-case hard-to-round
//! vectors, since the 2^64 input space cannot be enumerated.
//!
//! pfloat-libm is `no_std` + `alloc`. There is no alloc-free profile: correct
//! rounding grows the working precision at runtime past any compile-time
//! width, so the computation allocates.
//!
//! This slice wires the outer Ziv loop ([`round`]) and the full v0.1
//! surface ([`f32`], [`f64`]); the exhaustive `f32` sweep and `f64`
//! differential that certify it land in a following slice. See
//! `docs/kernel-list.md` for the surface and its verification tiers, and
//! the architecture decision records for the outer-loop design.

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub use pfloat::{RoundingMode, Status};

pub mod f32;
pub mod f64;
mod round;
mod saturate;
