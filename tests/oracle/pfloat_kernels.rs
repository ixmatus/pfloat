//! Pfloat-side kernel dispatch for the Phase 1 Oracle harness.
//!
//! For each `FnId` the harness verifies, this module routes to the
//! corresponding pfloat `BigFloat` method at `p = VERIFICATION_PRECISION`
//! (= 53, f64 precision) under the caller's rounding mode, returning
//! the result's binary32 bit pattern. The shape mirrors the MPFR
//! backend's `enclose` dispatch in `tests/oracle/mpfr.rs`.
//!
//! Slice p1.4 raised the kernel-call precision from `p = 24` (f32) to
//! `p = 53` (f64). The input bit pattern still maps to its exact f32
//! value via [`bf24_of_bits`]; converting that value to `p = 53` is
//! lossless (no f32 has more than 24 bits of significand). The
//! kernel then runs at `p = 53` and the result rounds to f32 via the
//! [`bf_to_f32_bits`] decimal bridge.
//!
//! Why `p = 53` rather than the original `p = 24`: for `f32`
//! *subnormal* outputs, the `BigFloat` at `p = 24` carries fewer
//! bits per ULP than the `f32` subnormal grid spacing, so the round
//! to `p = 24` can land the value on an f32-subnormal-grid midpoint
//! even when the true value sits just off the midpoint. The
//! subsequent `f32` conversion ties to even and may pick the wrong
//! neighbor. At `p = 53` the `BigFloat` ULP at any f32-subnormal
//! exponent is far finer than the `f32` grid, so the midpoint trap
//! never closes and the conversion lands on the `f32` the true
//! value is closer to (slice p1.4 closes pf-z0f for `erf`). Functions
//! whose tiny-input correction lives below `2^-(53 + 64) = 2^-117`
//! relative magnitude (notably `J1`'s cubic correction at
//! `~2^-298`) still need a kernel-side tiny-x precision boost; that
//! is wired separately in `bessel_j_tiny`.
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

/// Kernel-call precision for the f32 oracle pipeline. See the module
/// doc for the rationale (the BigFloat-at-p=24 → f32 conversion has
/// a midpoint-trap on f32 subnormals; p=53 sidesteps it).
const VERIFICATION_PRECISION: u32 = 53;

/// Evaluate pfloat's kernel for `f` at the binary32 bit pattern
/// `input` under `mode`. Returns the result's binary32 bit
/// pattern. The `f32 → BigFloat(p=53) → kernel → BigFloat(p=53)
/// → f32` round trip uses the bit-exact bridges from
/// `super::convert` so subnormal inputs and subnormal outputs are
/// preserved. The f32 → BigFloat(p=53) step lifts the exact f32
/// value into a 53-bit container (lossless because f32 carries at
/// most 24 bits of significand).
pub fn pfloat_kernel(f: FnId, input: u32, mode: RoundingMode) -> u32 {
    let x24 = bf24_of_bits(input);
    let x = x24
        .round_to_precision(VERIFICATION_PRECISION, RoundingMode::NearestEven)
        .expect("VERIFICATION_PRECISION >= 1; lift from p=24 is lossless")
        .0;
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
