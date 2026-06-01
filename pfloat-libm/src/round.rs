//! The outer correctly-rounded Ziv loop: `BigFloat -> f32/f64` with no
//! double rounding.
//!
//! A libm call widens its hardware float to an exact [`BigFloat`],
//! evaluates pfloat's correctly-rounded kernel at a working precision
//! `w`, and then rounds the result back to the hardware width. That
//! second rounding is the double-rounding hazard: rounding an
//! approximation to one width and then to the format can land on the
//! wrong neighbour. This module removes it by *enclosing* the true value
//! and committing a hardware float only once both ends of the enclosure
//! round there.
//!
//! # The directed enclosure
//!
//! At a working precision `w` ask the kernel for two values:
//!
//! - `lo`, the result rounded **toward negative infinity** to `w` bits:
//!   the largest `w`-bit value not exceeding `f(x)`.
//! - `hi`, the result rounded **toward positive infinity** to `w` bits:
//!   the smallest `w`-bit value not below `f(x)`.
//!
//! `[lo, hi]` is a true enclosure of `f(x)`, at most one `ulp_w` wide.
//! Round both ends to the format under the requested mode. Every IEEE 754
//! rounding mode is monotone non-decreasing, so when the two ends land on
//! the same hardware float the true value between them lands there too:
//! commit it. Otherwise grow `w`; the directed kernel calls make the
//! enclosure tighten around `f(x)` from both sides.
//!
//! The directed pair, rather than a single nearest-even call plus a
//! half-ULP bracket, is what makes directed rounding correct. A single
//! nearest-even call discards the sign of the residual `f(x) - r`: when
//! `f(x)` lies a hair off a grid point the kernel can round to that grid
//! point and even report the result exact (the residual falls below its
//! working precision), and a symmetric bracket around that point can
//! never decide which side `f(x)` is on. The directed calls force the
//! kernel's own Ziv loop to resolve the residual in each direction, so
//! the enclosure carries the information directed rounding needs. See
//! ADR-0057.
//!
//! # Schedule and the hard-to-round fallback
//!
//! The working precision is `prec + guard` for `guard` drawn in turn from
//! [`GUARD_SCHEDULE`], five iterations that mirror pfloat's own Ziv
//! schedule (`ZIV_BASE_GUARD = 64`, doubling, `ZIV_GUARD_CAP = 1024`,
//! `ZIV_MAX_ITERS = 5`; those constants are private to pfloat). On the
//! measure-zero hard-to-round input that exhausts the schedule the loop
//! returns the nearest-rounded best effort with `INEXACT` set, the caveat
//! pfloat and MPFR document. The exhaustive `f32` sweep and `f64`
//! differential against an independent oracle (the following slice) are
//! what certify the grid.
//!
//! # Status
//!
//! The returned [`Status`] is the union of the two kernel statuses and
//! the statuses of converting the two ends. The kernel contributes
//! `INVALID` and `DIV_BY_ZERO` from its domain handling; the conversions
//! contribute `INEXACT`, `OVERFLOW`, and `UNDERFLOW`. The union reports
//! `INEXACT` exactly when `f(x)` differs from the committed float: it is
//! clear only when both ends equal the committed value exactly (`f(x)` is
//! exactly that format float). The returned per-call status is
//! authoritative; under `std`, the conversions also raise into pfloat's
//! thread-local flag set, so that thread-local is not meaningful across a
//! shell call.

use pfloat::{BigFloat, Parts, RoundingMode, Status};

/// Guard bits tried above the format precision, in order. Mirrors
/// pfloat's private Ziv schedule: `ZIV_BASE_GUARD = 64`, doubling, capped
/// at `ZIV_GUARD_CAP = 1024`, for `ZIV_MAX_ITERS = 5` iterations
/// (`src/math/ziv.rs`).
const GUARD_SCHEDULE: [u32; 5] = [64, 128, 256, 512, 1024];

/// Format-specific glue the generic driver needs. Implemented by the
/// zero-sized [`F32Shell`] / [`F64Shell`] tags.
pub(crate) trait Shell {
    /// The hardware float type (`f32` or `f64`).
    type Hw: Copy;
    /// Significand precision including the implicit bit (24 / 53).
    const PREC: u32;

    /// Widen exactly to a `BigFloat` (lossless; precision [`PREC`]).
    fn widen(x: Self::Hw) -> BigFloat;
    /// Round to the format under `mode`, straight to the grid.
    fn convert(v: &BigFloat, mode: RoundingMode) -> (Self::Hw, Status);
    /// Bit-identity of two hardware floats (distinguishes `±0`).
    fn bits_eq(a: Self::Hw, b: Self::Hw) -> bool;
}

/// `f32` format tag.
pub(crate) struct F32Shell;
/// `f64` format tag.
pub(crate) struct F64Shell;

impl Shell for F32Shell {
    type Hw = f32;
    const PREC: u32 = 24;

    #[inline]
    fn widen(x: f32) -> BigFloat {
        BigFloat::from_f32(x)
    }
    #[inline]
    fn convert(v: &BigFloat, mode: RoundingMode) -> (f32, Status) {
        v.to_f32_round(mode)
    }
    #[inline]
    fn bits_eq(a: f32, b: f32) -> bool {
        a.to_bits() == b.to_bits()
    }
}

impl Shell for F64Shell {
    type Hw = f64;
    const PREC: u32 = 53;

    #[inline]
    fn widen(x: f64) -> BigFloat {
        BigFloat::from_f64(x)
    }
    #[inline]
    fn convert(v: &BigFloat, mode: RoundingMode) -> (f64, Status) {
        v.to_f64_round(mode)
    }
    #[inline]
    fn bits_eq(a: f64, b: f64) -> bool {
        a.to_bits() == b.to_bits()
    }
}

#[inline]
fn is_normal(v: &BigFloat) -> bool {
    matches!(v.parts(), Parts::Normal { .. })
}

/// Drive the outer Ziv loop for one hardware input.
///
/// `kernel(xb, w, dir)` returns `f(x)` correctly rounded to `w` bits
/// under the rounding mode `dir`, paired with the kernel's status. It is
/// supplied by each public entry point (it absorbs the `Result` /
/// non-`Result` split and captures any extra arguments).
pub(crate) fn drive<S: Shell>(
    x: S::Hw,
    mode: RoundingMode,
    kernel: impl Fn(&BigFloat, u32, RoundingMode) -> (BigFloat, Status),
) -> (S::Hw, Status) {
    let xb = S::widen(x);

    for &guard in &GUARD_SCHEDULE {
        let w = S::PREC + guard;
        let (lo, lo_s) = kernel(&xb, w, RoundingMode::TowardNegative);
        let (hi, hi_s) = kernel(&xb, w, RoundingMode::TowardPositive);

        if !is_normal(&lo) || !is_normal(&hi) {
            // NaN / ±inf / ±0: the kernel already applied IEEE semantics
            // and both directed roundings agree. Convert and merge.
            let (val, cs) = S::convert(&lo, mode);
            return (val, lo_s | hi_s | cs);
        }

        let (lo_f, lo_cs) = S::convert(&lo, mode);
        let (hi_f, hi_cs) = S::convert(&hi, mode);
        if S::bits_eq(lo_f, hi_f) {
            // The whole enclosure rounds to one float: settled.
            return (lo_f, lo_s | hi_s | lo_cs | hi_cs);
        }
    }

    // Schedule exhausted on a measure-zero hard-to-round input: best
    // effort is the nearest-rounded value at the finest precision, with
    // INEXACT set, mirroring pfloat's own ZIV_MAX_ITERS fallback.
    let w = S::PREC + GUARD_SCHEDULE[GUARD_SCHEDULE.len() - 1];
    let (ne, ne_s) = kernel(&xb, w, RoundingMode::NearestEven);
    let (val, cs) = S::convert(&ne, mode);
    (val, ne_s | cs | Status::INEXACT)
}

/// Emit a hardware-float libm entry pair for one unary kernel: a
/// nearest-even convenience `fn $name(x) -> Hw` and a mode-aware
/// `fn $round(x, mode) -> (Hw, Status)`. The `$round` metavariable names
/// both the public mode-aware entry and the `BigFloat` kernel method
/// (they share the spelling, e.g. `exp_round`); the `$round(...)` inside
/// the convenience fn is the sibling free fn, while `xb.$round(...)` is
/// the inherent method. The `result` / `direct` discriminator selects
/// whether the kernel method returns a `Result` (every kernel but `exp`).
macro_rules! unary {
    ($hw:ty, $shell:ident, $name:ident, $round:ident, result, $disp:literal) => {
        #[doc = concat!($disp, " correctly rounded to the format under round-to-nearest-even.")]
        #[must_use]
        pub fn $name(x: $hw) -> $hw {
            $round(x, $crate::RoundingMode::NearestEven).0
        }

        #[doc = concat!($disp, " correctly rounded under `mode`, with IEEE 754 status flags.")]
        #[must_use]
        pub fn $round(x: $hw, mode: $crate::RoundingMode) -> ($hw, $crate::Status) {
            $crate::round::drive::<$crate::round::$shell>(x, mode, |xb, w, dir| {
                xb.$round(w, dir)
                    .expect("w = PREC + guard >= 1: BuildError only on precision 0")
            })
        }
    };
    ($hw:ty, $shell:ident, $name:ident, $round:ident, direct, $disp:literal) => {
        #[doc = concat!($disp, " correctly rounded to the format under round-to-nearest-even.")]
        #[must_use]
        pub fn $name(x: $hw) -> $hw {
            $round(x, $crate::RoundingMode::NearestEven).0
        }

        #[doc = concat!($disp, " correctly rounded under `mode`, with IEEE 754 status flags.")]
        #[must_use]
        pub fn $round(x: $hw, mode: $crate::RoundingMode) -> ($hw, $crate::Status) {
            $crate::round::drive::<$crate::round::$shell>(x, mode, |xb, w, dir| xb.$round(w, dir))
        }
    };
}

pub(crate) use unary;

#[cfg(test)]
mod tests {
    use super::*;

    // The non-circular gold standard: widen exactly, evaluate the kernel
    // at 2048 bits (far more than any hardware float needs), and round
    // once to the format. `2048` dwarfs the shell's working precision, so
    // for all but a measure-zero set of hard-to-round inputs this single
    // rounding is the correct result. It is not an independent oracle
    // (both sides use pfloat's kernel), so it certifies the shell's
    // *rounding logic*, not the kernel; the exhaustive sweep against
    // MPFR/Arb (next slice) certifies the kernel.
    const ORACLE_PREC: u32 = 2048;

    /// Samples per strided sweep. A smoke over the range that catches
    /// rounding-logic bugs in every mode; the exhaustive 2^32 density is
    /// the next slice's (pf-lm3) job, not this routine unit gate.
    const SWEEP_SAMPLES: u32 = 300;

    const MODES: [RoundingMode; 5] = [
        RoundingMode::NearestEven,
        RoundingMode::NearestAway,
        RoundingMode::TowardZero,
        RoundingMode::TowardPositive,
        RoundingMode::TowardNegative,
    ];

    fn oracle_f32(x: f32, mode: RoundingMode, kr: fn(&BigFloat) -> BigFloat) -> (f32, Status) {
        kr(&BigFloat::from_f32(x)).to_f32_round(mode)
    }

    fn oracle_f64(x: f64, mode: RoundingMode, kr: fn(&BigFloat) -> BigFloat) -> (f64, Status) {
        kr(&BigFloat::from_f64(x)).to_f64_round(mode)
    }

    fn check_f32(
        shell: fn(f32, RoundingMode) -> (f32, Status),
        kr: fn(&BigFloat) -> BigFloat,
        extra: &[f32],
    ) {
        let battery = [
            0.0f32,
            -0.0,
            f32::from_bits(0x0000_0001), // smallest subnormal
            f32::from_bits(0x007F_FFFF), // largest subnormal
            f32::from_bits(0x0080_0000), // smallest normal
            1.0,
            -1.0,
            2.0,
            0.5,
            f32::MAX,
            -f32::MAX,
        ];
        for &x in battery.iter().chain(extra) {
            for &mode in &MODES {
                let (got, _) = shell(x, mode);
                let (want, _) = oracle_f32(x, mode, kr);
                assert!(
                    got.to_bits() == want.to_bits() || (got.is_nan() && want.is_nan()),
                    "f32 x={x} ({:#010x}) mode={mode:?}: got {:#010x}, want {:#010x}",
                    x.to_bits(),
                    got.to_bits(),
                    want.to_bits()
                );
            }
        }
        let mut bits: u32 = 0x1357_9BDF;
        for i in 0..SWEEP_SAMPLES {
            let x = f32::from_bits(bits);
            let mode = MODES[(i % 5) as usize];
            let (got, _) = shell(x, mode);
            let (want, _) = oracle_f32(x, mode, kr);
            assert!(
                got.to_bits() == want.to_bits() || (got.is_nan() && want.is_nan()),
                "f32 sweep x={x} ({bits:#010x}) mode={mode:?}: got {:#010x}, want {:#010x}",
                got.to_bits(),
                want.to_bits()
            );
            bits = bits.wrapping_add(0x9E37_79B1);
        }
    }

    fn check_f64(
        shell: fn(f64, RoundingMode) -> (f64, Status),
        kr: fn(&BigFloat) -> BigFloat,
        extra: &[f64],
    ) {
        let battery = [
            0.0f64,
            -0.0,
            f64::from_bits(0x0000_0000_0000_0001),
            f64::from_bits(0x000F_FFFF_FFFF_FFFF),
            f64::from_bits(0x0010_0000_0000_0000),
            1.0,
            -1.0,
            2.0,
            0.5,
            f64::MAX,
            -f64::MAX,
        ];
        for &x in battery.iter().chain(extra) {
            for &mode in &MODES {
                let (got, _) = shell(x, mode);
                let (want, _) = oracle_f64(x, mode, kr);
                assert!(
                    got.to_bits() == want.to_bits() || (got.is_nan() && want.is_nan()),
                    "f64 x={x} mode={mode:?}: got {:#018x}, want {:#018x}",
                    got.to_bits(),
                    want.to_bits()
                );
            }
        }
        let mut bits: u64 = 0x1357_9BDF_2468_ACE0;
        for i in 0..SWEEP_SAMPLES {
            let x = f64::from_bits(bits);
            let mode = MODES[(i % 5) as usize];
            let (got, _) = shell(x, mode);
            let (want, _) = oracle_f64(x, mode, kr);
            assert!(
                got.to_bits() == want.to_bits() || (got.is_nan() && want.is_nan()),
                "f64 sweep x={x} ({bits:#018x}) mode={mode:?}: got {:#018x}, want {:#018x}",
                got.to_bits(),
                want.to_bits()
            );
            bits = bits.wrapping_add(0x9E37_79B9_7F4A_7C15);
        }
    }

    #[test]
    fn exp_matches_single_rounding() {
        check_f32(
            crate::f32::exp_round,
            |xb| xb.exp_round(ORACLE_PREC, RoundingMode::NearestEven).0,
            &[0.5, 1.0, 10.0, -10.0, 88.0, -88.0],
        );
        check_f64(
            crate::f64::exp_round,
            |xb| xb.exp_round(ORACLE_PREC, RoundingMode::NearestEven).0,
            &[1.0, 100.0, 700.0, -700.0],
        );
    }

    #[test]
    fn ln_matches_single_rounding() {
        check_f32(
            crate::f32::ln_round,
            |xb| {
                xb.ln_round(ORACLE_PREC, RoundingMode::NearestEven)
                    .unwrap()
                    .0
            },
            &[2.0, 10.0, 0.1, 1.5],
        );
        check_f64(
            crate::f64::ln_round,
            |xb| {
                xb.ln_round(ORACLE_PREC, RoundingMode::NearestEven)
                    .unwrap()
                    .0
            },
            &[2.0, 10.0, 0.5],
        );
    }

    #[test]
    fn sqrt_matches_single_rounding() {
        check_f32(
            crate::f32::sqrt_round,
            |xb| {
                xb.sqrt_round(ORACLE_PREC, RoundingMode::NearestEven)
                    .unwrap()
                    .0
            },
            &[2.0, 3.0, 0.5],
        );
        check_f64(
            crate::f64::sqrt_round,
            |xb| {
                xb.sqrt_round(ORACLE_PREC, RoundingMode::NearestEven)
                    .unwrap()
                    .0
            },
            &[2.0, 3.0],
        );
    }

    #[test]
    fn sin_matches_single_rounding() {
        check_f32(
            crate::f32::sin_round,
            |xb| {
                xb.sin_round(ORACLE_PREC, RoundingMode::NearestEven)
                    .unwrap()
                    .0
            },
            &[0.5, 1.0, 3.0, -3.0],
        );
    }

    #[test]
    fn cbrt_matches_single_rounding() {
        check_f32(
            crate::f32::cbrt_round,
            |xb| {
                xb.cbrt_round(ORACLE_PREC, RoundingMode::NearestEven)
                    .unwrap()
                    .0
            },
            &[8.0, -8.0, 2.0, -2.0],
        );
    }

    #[test]
    fn cot_matches_single_rounding() {
        // cot/sec/csc have no MPFR primitive; the in-crate single-rounding
        // oracle still exercises the shell's rounding logic for them.
        check_f32(
            crate::f32::cot_round,
            |xb| {
                xb.cot_round(ORACLE_PREC, RoundingMode::NearestEven)
                    .unwrap()
                    .0
            },
            &[0.5, 1.0, 2.0],
        );
    }
}
