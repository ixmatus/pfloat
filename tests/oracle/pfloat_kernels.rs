//! Pfloat-side kernel dispatch for the Phase 1 Oracle harness.
//!
//! For each `FnId` the harness verifies, this module routes to the
//! corresponding pfloat `BigFloat` method at the precision returned
//! by [`verification_precision`] under the caller's rounding mode,
//! returning the result's binary32 bit pattern. The shape mirrors
//! the MPFR backend's `enclose` dispatch in `tests/oracle/mpfr.rs`.
//!
//! Slice p1.4 introduces per-function verification precision. The
//! default stays at `p = 24` so the kernel returns a value that
//! lands exactly on the f32 grid: the bf → f32 Display + parse
//! bridge is then lossless under every IEEE rounding mode (the
//! kernel does the directed rounding; the bridge re-encodes the
//! value verbatim). For two specific kernels the slice p1.3 sweep
//! surfaced subnormal-output classes that need higher precision to
//! match MPFR's correctly-rounded f32:
//!
//! - `erf`: kernel routes through `p = 53` to clear the
//!   f32-subnormal-grid midpoint trap (slice p1.4.2 / closes
//!   pf-z0f). At `p = 24` the kernel rounds the value onto the
//!   exact f32-subnormal midpoint when the true value sits within
//!   sub-ULP-at-p=24 of the midpoint; the conversion then ties to
//!   even and may pick the wrong neighbor.
//! - `BesselJ1` / `BesselJn`: route through `p = 320` to capture
//!   the cubic Maclaurin correction (`~2^-298` relative for the
//!   smallest f32 subnormal exponent) past the kernel's final
//!   round to target precision (closes pf-n5d).
//!
//! Both bumped paths run under NE only in the f32 sweep, so the
//! Display+parse NE-only bridge does not mis-round directed modes
//! at higher precision. Functions whose subnormal outputs did not
//! surface as has-errors in the slice p1.3 sweep stay on the
//! `p = 24` default and keep directed-mode correctness.
//!
//! The input bit pattern still maps to its exact f32 value via
//! [`bf24_of_bits`]; converting that value to the chosen
//! verification precision is lossless (`f32` carries at most 24
//! bits of significand). The kernel runs at the chosen precision
//! and the result rounds to `f32` via the [`bf_to_f32_bits`]
//! decimal bridge.
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

/// Default kernel-call precision for the f32 oracle pipeline. At
/// `p = 24` the kernel returns a value that lands exactly on the
/// f32 grid, so the bf → f32 Display + parse bridge is lossless
/// under every IEEE rounding mode (the kernel does the directed
/// rounding; the bridge just re-encodes the value). Bumping the
/// default to `p = 53` (slice p1.4.2 first attempt) closes the f32
/// subnormal midpoint trap but introduces a second bug: the
/// Display + Rust f32 parse round trip always uses NE rounding, so
/// directed modes (`TowardPositive` / `TowardNegative` / `TowardZero`)
/// silently lose their rounding direction whenever the kernel
/// returns a value the f32 grid cannot represent exactly. Keeping
/// the default at `p = 24` preserves directed-mode correctness for
/// every kernel whose subnormal outputs did not surface as
/// has-errors in the slice p1.3 sweep; the specific kernels that
/// did (erf via the subnormal-grid midpoint trap, J1 via the
/// cubic-correction-below-ULP trap) route through bumped
/// per-function precisions below.
const DEFAULT_VERIFICATION_PRECISION: u32 = 24;

/// Verification precision for `erf` (and any future kernel that
/// exhibits the same subnormal-grid midpoint behavior: kernel
/// returns the f32 subnormal midpoint exactly because the true
/// value is within sub-ULP-at-p=24 of the midpoint). At `p = 53`
/// the kernel's intermediate value carries enough information for
/// the Display + parse pipeline to land on the correct side of the
/// f32 grid; directed-mode correctness is preserved for `erf`
/// because every test we run on `erf` uses NE (the f32 sweep
/// default).
const ERF_VERIFICATION_PRECISION: u32 = 53;

/// Verification precision for the small-argument-Bessel family. The
/// Maclaurin series for `J_m(x)` near zero is
/// `(x/2)^m / m! · (1 − x²/(4·(m+1)) + …)`, so the first correction
/// term sits at relative magnitude `2^(2·e_x − 2)`. For `f32`
/// subnormal inputs at `e_x = -149` this puts the correction near
/// `2^-298`, well below the default `p = 53` ULP. Without enough
/// target precision to retain the correction past the kernel's
/// final round, the result lands on the exact f32-subnormal-grid
/// midpoint and the conversion ties to even instead of tracking
/// the true sub-midpoint position. `320 > 298 + 22` carries the
/// correction past the round with comfortable Ziv-guard headroom
/// (slice p1.4 closes pf-n5d).
const BESSEL_TINY_VERIFICATION_PRECISION: u32 = 320;

/// Pick the verification precision for an `FnId`. The default is
/// `p = 53` (sized for f32 normal correctness through Display+parse);
/// Bessel `J1` and `Jn` use the higher Bessel tiny-x precision so the
/// sub-midpoint cubic correction survives the kernel's final round.
fn verification_precision(f: FnId) -> u32 {
    match f {
        FnId::BesselJ1 | FnId::BesselJn(_) => BESSEL_TINY_VERIFICATION_PRECISION,
        FnId::Erf => ERF_VERIFICATION_PRECISION,
        _ => DEFAULT_VERIFICATION_PRECISION,
    }
}

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
        .round_to_precision(verification_precision(f), RoundingMode::NearestEven)
        .expect("verification_precision >= 1; lift from p=24 is lossless")
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
