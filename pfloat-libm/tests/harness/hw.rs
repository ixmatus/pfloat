//! The hardware-width abstraction.
//!
//! [`Hw`] is implemented by `f32` and `f64`. It carries everything the
//! width-generic verifier and driver need: lifting an input to an exact
//! `rug::Float`, rounding a `rug::Float` back to the width, and calling
//! the shell entry point under test. One implementation of the Ziv-at-
//! oracle loop then serves both widths, so the f32 and f64 lanes cannot
//! drift apart. This is the harness's analogue of the shell's own
//! `Shell` trait (`pfloat-libm/src/round.rs`).

#![cfg(all(unix, feature = "differential-mpfr"))]

use pfloat_libm::{f32 as lm32, f64 as lm64, RoundingMode, Status};
use rug::Float;

use super::convert::{round_f32, round_f64};
use super::types::{LibmArg, LibmFnId, Width};

/// A hardware float the harness can verify against the MPFR oracle.
pub trait Hw: Copy + 'static {
    /// The width's bit pattern (`u32` for `f32`, `u64` for `f64`).
    type Bits: Copy + Eq + core::fmt::Debug;
    const WIDTH: Width;

    fn from_bits(b: Self::Bits) -> Self;
    fn to_bits(self) -> Self::Bits;
    fn is_nan(self) -> bool;
    fn is_finite(self) -> bool;
    /// Widen exactly to `f64` for the pole classifier (every `f32` is an
    /// exact `f64`; `f64` is the identity).
    fn as_f64(self) -> f64;
    /// Widen the bit pattern to `u64` (f32 in the low 32 bits) for the
    /// width-tagged `Verdict`/corpus.
    fn bits_to_u64(b: Self::Bits) -> u64;
    /// Reconstruct this width's bit pattern from a `u64` sweep counter
    /// (f32 takes the low 32 bits, which is exact within `[0, 2^32)`).
    fn bits_from_u64(u: u64) -> Self::Bits;
    /// Cast a binary64 Lefevre-Muller seed (given as `f64` bits) to this
    /// width's nearest value's bits, for folding the corpus into a sweep.
    fn seed_from_f64_bits(seed: u64) -> Self::Bits;
    /// Extract this width's partner bits from a `HypotY` payload.
    fn partner_bits(yb: u64) -> Self::Bits;
    /// Lift the input exactly to a `rug::Float` at `prec` bits
    /// (`prec >= 24/53`, so the lift is lossless).
    fn lift(b: Self::Bits, prec: u32) -> Float;
    /// Round a `rug::Float` to this width under `mode`.
    fn round(v: &Float, mode: RoundingMode) -> Option<Self>;
    fn nan() -> Self;
    /// Call the shell entry point for `f` under `mode`, returning the
    /// result's bits and the IEEE status it raised.
    fn shell(f: LibmFnId, b: Self::Bits, arg: LibmArg, mode: RoundingMode) -> (Self::Bits, Status);
}

impl Hw for f32 {
    type Bits = u32;
    const WIDTH: Width = Width::F32;

    fn from_bits(b: u32) -> f32 {
        f32::from_bits(b)
    }
    fn to_bits(self) -> u32 {
        f32::to_bits(self)
    }
    fn is_nan(self) -> bool {
        f32::is_nan(self)
    }
    fn is_finite(self) -> bool {
        f32::is_finite(self)
    }
    fn as_f64(self) -> f64 {
        self as f64
    }
    fn bits_to_u64(b: u32) -> u64 {
        u64::from(b)
    }
    fn bits_from_u64(u: u64) -> u32 {
        u as u32
    }
    fn seed_from_f64_bits(seed: u64) -> u32 {
        (f64::from_bits(seed) as f32).to_bits()
    }
    fn partner_bits(yb: u64) -> u32 {
        yb as u32
    }
    fn lift(b: u32, prec: u32) -> Float {
        Float::with_val(prec, f32::from_bits(b))
    }
    fn round(v: &Float, mode: RoundingMode) -> Option<f32> {
        round_f32(v, mode)
    }
    fn nan() -> f32 {
        f32::NAN
    }
    fn shell(f: LibmFnId, b: u32, arg: LibmArg, mode: RoundingMode) -> (u32, Status) {
        let x = f32::from_bits(b);
        let (v, s): (f32, Status) = match f {
            LibmFnId::Exp => lm32::exp_round(x, mode),
            LibmFnId::Exp2 => lm32::exp2_round(x, mode),
            LibmFnId::Exp10 => lm32::exp10_round(x, mode),
            LibmFnId::Expm1 => lm32::expm1_round(x, mode),
            LibmFnId::Ln => lm32::ln_round(x, mode),
            LibmFnId::Log2 => lm32::log2_round(x, mode),
            LibmFnId::Log10 => lm32::log10_round(x, mode),
            LibmFnId::Log1p => lm32::log1p_round(x, mode),
            LibmFnId::Sqrt => lm32::sqrt_round(x, mode),
            LibmFnId::Cbrt => lm32::cbrt_round(x, mode),
            LibmFnId::Sin => lm32::sin_round(x, mode),
            LibmFnId::Cos => lm32::cos_round(x, mode),
            LibmFnId::Tan => lm32::tan_round(x, mode),
            LibmFnId::Cot => lm32::cot_round(x, mode),
            LibmFnId::Sec => lm32::sec_round(x, mode),
            LibmFnId::Csc => lm32::csc_round(x, mode),
            LibmFnId::Asin => lm32::asin_round(x, mode),
            LibmFnId::Acos => lm32::acos_round(x, mode),
            LibmFnId::Atan => lm32::atan_round(x, mode),
            LibmFnId::Sinh => lm32::sinh_round(x, mode),
            LibmFnId::Cosh => lm32::cosh_round(x, mode),
            LibmFnId::Tanh => lm32::tanh_round(x, mode),
            LibmFnId::Asinh => lm32::asinh_round(x, mode),
            LibmFnId::Acosh => lm32::acosh_round(x, mode),
            LibmFnId::Atanh => lm32::atanh_round(x, mode),
            LibmFnId::Hypot => {
                let LibmArg::HypotY(yb) = arg else {
                    panic!("hypot requires a HypotY arg");
                };
                lm32::hypot_round(x, f32::from_bits(yb as u32), mode)
            }
            LibmFnId::Rootn(n) => lm32::rootn_round(x, n, mode),
        };
        (v.to_bits(), s)
    }
}

impl Hw for f64 {
    type Bits = u64;
    const WIDTH: Width = Width::F64;

    fn from_bits(b: u64) -> f64 {
        f64::from_bits(b)
    }
    fn to_bits(self) -> u64 {
        f64::to_bits(self)
    }
    fn is_nan(self) -> bool {
        f64::is_nan(self)
    }
    fn is_finite(self) -> bool {
        f64::is_finite(self)
    }
    fn as_f64(self) -> f64 {
        self
    }
    fn bits_to_u64(b: u64) -> u64 {
        b
    }
    fn bits_from_u64(u: u64) -> u64 {
        u
    }
    fn seed_from_f64_bits(seed: u64) -> u64 {
        seed
    }
    fn partner_bits(yb: u64) -> u64 {
        yb
    }
    fn lift(b: u64, prec: u32) -> Float {
        Float::with_val(prec, f64::from_bits(b))
    }
    fn round(v: &Float, mode: RoundingMode) -> Option<f64> {
        round_f64(v, mode)
    }
    fn nan() -> f64 {
        f64::NAN
    }
    fn shell(f: LibmFnId, b: u64, arg: LibmArg, mode: RoundingMode) -> (u64, Status) {
        let x = f64::from_bits(b);
        let (v, s): (f64, Status) = match f {
            LibmFnId::Exp => lm64::exp_round(x, mode),
            LibmFnId::Exp2 => lm64::exp2_round(x, mode),
            LibmFnId::Exp10 => lm64::exp10_round(x, mode),
            LibmFnId::Expm1 => lm64::expm1_round(x, mode),
            LibmFnId::Ln => lm64::ln_round(x, mode),
            LibmFnId::Log2 => lm64::log2_round(x, mode),
            LibmFnId::Log10 => lm64::log10_round(x, mode),
            LibmFnId::Log1p => lm64::log1p_round(x, mode),
            LibmFnId::Sqrt => lm64::sqrt_round(x, mode),
            LibmFnId::Cbrt => lm64::cbrt_round(x, mode),
            LibmFnId::Sin => lm64::sin_round(x, mode),
            LibmFnId::Cos => lm64::cos_round(x, mode),
            LibmFnId::Tan => lm64::tan_round(x, mode),
            LibmFnId::Cot => lm64::cot_round(x, mode),
            LibmFnId::Sec => lm64::sec_round(x, mode),
            LibmFnId::Csc => lm64::csc_round(x, mode),
            LibmFnId::Asin => lm64::asin_round(x, mode),
            LibmFnId::Acos => lm64::acos_round(x, mode),
            LibmFnId::Atan => lm64::atan_round(x, mode),
            LibmFnId::Sinh => lm64::sinh_round(x, mode),
            LibmFnId::Cosh => lm64::cosh_round(x, mode),
            LibmFnId::Tanh => lm64::tanh_round(x, mode),
            LibmFnId::Asinh => lm64::asinh_round(x, mode),
            LibmFnId::Acosh => lm64::acosh_round(x, mode),
            LibmFnId::Atanh => lm64::atanh_round(x, mode),
            LibmFnId::Hypot => {
                let LibmArg::HypotY(yb) = arg else {
                    panic!("hypot requires a HypotY arg");
                };
                lm64::hypot_round(x, f64::from_bits(yb), mode)
            }
            LibmFnId::Rootn(n) => lm64::rootn_round(x, n, mode),
        };
        (v.to_bits(), s)
    }
}
