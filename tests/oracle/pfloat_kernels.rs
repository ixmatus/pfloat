//! Pfloat-side kernel dispatch for the Phase 1 Oracle harness.
//!
//! For each `FnId` the harness verifies, this module routes to the
//! corresponding pfloat `BigFloat` method at `p = 24` (f32
//! precision) under the caller's rounding mode, returning the
//! result's binary32 bit pattern. The shape mirrors the MPFR
//! backend's `enclose` dispatch in `tests/oracle/mpfr.rs`.
//!
//! All 47 frozen v1.0 surface entries are wired (the surface that
//! `docs/v1.0-surface.md` enumerates). The Arb-primary entries (no
//! MPFR primitive) are wired here even though slice p1.3's MPFR
//! backend does not verify them; the Arb backend (next slice) will
//! route those same `FnId` variants without changing the pfloat
//! dispatch.

#![cfg(all(unix, feature = "differential-mpfr"))]

use pfloat::RoundingMode;

use super::convert::{bf24_of_bits, bf_to_f32_bits};
use super::types::FnId;

/// Evaluate pfloat's kernel for `f` at the binary32 bit pattern
/// `input` under `mode`. Returns the result's binary32 bit
/// pattern. The `f32 → BigFloat → kernel → BigFloat → f32` round
/// trip uses the bit-exact bridges from `super::convert` so
/// subnormal inputs and subnormal outputs are preserved.
pub fn pfloat_kernel(f: FnId, input: u32, mode: RoundingMode) -> u32 {
    let x = bf24_of_bits(input);
    let result = match f {
        // Elementary.
        FnId::Sqrt => x.sqrt(mode).0,
        FnId::Exp => x.exp(mode).0,
        FnId::Exp2 => x.exp2(mode).0,
        FnId::Exp10 => x.exp10(mode).0,
        FnId::Expm1 => x.expm1(mode).0,
        FnId::Ln => x.ln(mode).0,
        FnId::Log1p => x.log1p(mode).0,
        FnId::Log2 => x.log2(mode).0,
        FnId::Log10 => x.log10(mode).0,
        FnId::Sin => x.sin(mode).0,
        FnId::Cos => x.cos(mode).0,
        FnId::Tan => x.tan(mode).0,
        FnId::Asin => x.asin(mode).0,
        FnId::Acos => x.acos(mode).0,
        FnId::Atan => x.atan(mode).0,
        FnId::Sinh => x.sinh(mode).0,
        FnId::Cosh => x.cosh(mode).0,
        FnId::Tanh => x.tanh(mode).0,
        FnId::Asinh => x.asinh(mode).0,
        FnId::Acosh => x.acosh(mode).0,
        FnId::Atanh => x.atanh(mode).0,
        // Specials.
        FnId::Erf => x.erf(mode).0,
        FnId::Erfc => x.erfc(mode).0,
        FnId::Gamma => x.gamma(mode).0,
        FnId::Lgamma => x.lgamma(mode).0,
        FnId::Digamma => x.digamma(mode).0,
        FnId::Zeta => x.zeta(mode).0,
        FnId::Ei => x.ei(mode).0,
        FnId::Si => x.si(mode).0,
        FnId::Ci => x.ci(mode).0,
        FnId::Li => x.li(mode).0,
        // Airy.
        FnId::Ai => x.ai(mode).0,
        FnId::Bi => x.bi(mode).0,
        FnId::AiPrime => x.ai_prime(mode).0,
        FnId::BiPrime => x.bi_prime(mode).0,
        // Bessel J / Y.
        FnId::BesselJ0 => x.j0(mode).0,
        FnId::BesselJ1 => x.j1(mode).0,
        FnId::BesselJn(n) => x.jn(n, mode).0,
        FnId::BesselY0 => x.y0(mode).0,
        FnId::BesselY1 => x.y1(mode).0,
        FnId::BesselYn(n) => x.yn(n, mode).0,
        // Bessel I / K.
        FnId::BesselI0 => x.i0(mode).0,
        FnId::BesselI1 => x.i1(mode).0,
        FnId::BesselIn(n) => x.in_(n, mode).0,
        FnId::BesselK0 => x.k0(mode).0,
        FnId::BesselK1 => x.k1(mode).0,
        FnId::BesselKn(n) => x.kn(n, mode).0,
    };
    bf_to_f32_bits(&result)
}
