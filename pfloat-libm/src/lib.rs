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
//! This crate is a scaffold under active construction; see
//! `docs/kernel-list.md` for the planned surface and its verification tiers.

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;
