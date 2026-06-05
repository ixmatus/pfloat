//! Saturation fast-paths for the over/underflowing transcendentals.
//!
//! `exp`, `exp2`, `exp10`, `expm1`, `sinh`, `cosh`, and `tanh` are the
//! functions whose result saturates (to `±∞`/`±MAX`, `±0`, `±1`) once the
//! argument is large enough. Evaluating the kernel there is both wasteful
//! and very slow: pfloat's argument reduction cost grows with the
//! argument's exponent, so `exp(2^113)` costs ~700 µs and `tanh(2^113)`
//! ~34 ms (it forms `e^(2·2^113)` internally), where the answer is a fixed
//! saturated value. Brute-forcing those inputs makes the exhaustive `f32`
//! sweep infeasible (a single huge-magnitude shard runs for days).
//!
//! This module returns the saturated result directly, without the kernel,
//! for inputs beyond a per-width threshold set above every break point of
//! the family. The kernel still handles the boundary region (normal
//! magnitude, so it is fast there). The saturated values and status flags
//! reproduce exactly what the shell already returns for these inputs:
//! `over_pos`/`over_neg`/`under_pos` mirror pfloat's conversion
//! (`src/convert.rs` `overflow` / `tiny`), and `pos_one`/`neg_one` mirror
//! the kernel's `±1` saturation (past the threshold the residual has
//! fallen below the kernel's working precision, so it rounds `tanh` /
//! `expm1` to `±1`, with `INEXACT` because the true value is strictly
//! inside the saturation limit; ADR-0063). The fast-path is therefore a
//! transparent, behavior-preserving speedup: bit-for-bit identical to the
//! kernel path for every input it handles. That identity is the
//! load-bearing claim and is checked against the kernel in the crate tests.
//!
//! Non-finite inputs are never fast-pathed (the guard requires a finite
//! argument): `exp(∞) = ∞` with no `OVERFLOW`, `tanh(∞) = 1` exactly with
//! no `INEXACT`, NaN propagates — all of which the kernel handles
//! correctly, and none of which the saturation values would match.

use crate::{RoundingMode, Status};

/// Hardware float that can name its saturated values. `over_pos`/
/// `over_neg` mirror `src/convert.rs` `overflow()`, `under_pos` mirrors
/// `tiny()` for a vanishing positive value, and `pos_one`/`neg_one` are
/// the kernel's exact `±1` saturation.
pub(crate) trait SatHw: Copy {
    /// Magnitude at or beyond which the kernel's result is both saturated
    /// AND uniform for every saturating function of this width, so the
    /// fast-path can name it without consulting the kernel. It is set
    /// above every break point of the family: the exp-family
    /// overflow/underflow points (widest is `exp2`, 128/1024) and the
    /// argument past which `tanh`/`expm1` saturate to exactly `±1` because
    /// their residual falls below the kernel's largest working precision
    /// (widest is `expm1` near ~726 for `f32`, ~746 for `f64`). Below the
    /// threshold the kernel still owns the band (and is fast there); the
    /// slow huge-magnitude inputs (|x| >~ 2^49) are far beyond it. Keeping
    /// the threshold above every break point is what makes the fast-path a
    /// transparent, behavior-preserving speedup: it reproduces exactly
    /// what the shell already returns for these inputs.
    const THRESHOLD: f64;

    fn is_finite(self) -> bool;
    /// Exact widening to `f64` for the threshold comparison (`f32 -> f64`
    /// is lossless; `f64` is the identity).
    fn as_f64(self) -> f64;

    /// `+∞` (NE/NA/TP) or `+MAX` (TZ/TN), `OVERFLOW | INEXACT`.
    fn over_pos(mode: RoundingMode) -> (Self, Status);
    /// `-∞` (NE/NA/TN) or `-MAX` (TZ/TP), `OVERFLOW | INEXACT`.
    fn over_neg(mode: RoundingMode) -> (Self, Status);
    /// `+0` (NE/NA/TZ/TN) or the smallest positive subnormal (TP),
    /// `UNDERFLOW | INEXACT`.
    fn under_pos(mode: RoundingMode) -> (Self, Status);
    /// Exactly `+1` with `INEXACT`, for every mode. Past the threshold the
    /// kernel rounds `tanh` to `1` (its residual `2·e^{-2x}` has fallen
    /// below the working precision), but the true value is strictly inside
    /// `(−1, 1)`, so the result is inexact (ADR-0063, resolving the
    /// pf-njs5 question the fast-path previously deferred). The kernel now
    /// forces `INEXACT` there, and the fast-path reproduces it.
    fn pos_one() -> (Self, Status);
    /// Exactly `-1` with `INEXACT`, for every mode (the negative-saturation
    /// analogue of [`SatHw::pos_one`], shared by `tanh` and `expm1`; the
    /// true value is strictly above `-1`).
    fn neg_one() -> (Self, Status);
}

/// Generate the `SatHw` impl for one width. The `over_*`/`under_pos` arms
/// encode exactly `convert::overflow` / `convert::tiny`; `pos_one`/
/// `neg_one` encode the kernel's exact `±1` saturation.
macro_rules! impl_sat_hw {
    ($hw:ty, $threshold:expr, $min_subnormal:expr) => {
        impl SatHw for $hw {
            const THRESHOLD: f64 = $threshold;
            #[inline]
            fn is_finite(self) -> bool {
                <$hw>::is_finite(self)
            }
            #[inline]
            fn as_f64(self) -> f64 {
                self as f64
            }
            #[inline]
            fn over_pos(mode: RoundingMode) -> ($hw, Status) {
                let v = match mode {
                    RoundingMode::NearestEven
                    | RoundingMode::NearestAway
                    | RoundingMode::TowardPositive => <$hw>::INFINITY,
                    RoundingMode::TowardZero | RoundingMode::TowardNegative => <$hw>::MAX,
                };
                (v, Status::OVERFLOW | Status::INEXACT)
            }
            #[inline]
            fn over_neg(mode: RoundingMode) -> ($hw, Status) {
                let v = match mode {
                    RoundingMode::NearestEven
                    | RoundingMode::NearestAway
                    | RoundingMode::TowardNegative => <$hw>::NEG_INFINITY,
                    RoundingMode::TowardZero | RoundingMode::TowardPositive => -<$hw>::MAX,
                };
                (v, Status::OVERFLOW | Status::INEXACT)
            }
            #[inline]
            fn under_pos(mode: RoundingMode) -> ($hw, Status) {
                let v = match mode {
                    RoundingMode::TowardPositive => $min_subnormal,
                    _ => 0.0,
                };
                (v, Status::UNDERFLOW | Status::INEXACT)
            }
            #[inline]
            fn pos_one() -> ($hw, Status) {
                (1.0, Status::INEXACT)
            }
            #[inline]
            fn neg_one() -> ($hw, Status) {
                (-1.0, Status::INEXACT)
            }
        }
    };
}

// f32 threshold 1024: past every f32 break point — exp-family overflow
// (exp2 at 128), underflow (exp2 at -150), and the `tanh`/`expm1` exact-
// `±1` saturation (~726).
impl_sat_hw!(
    f32,
    1024.0,
    f32::from_bits(0x0000_0001) // 2^-149, smallest positive subnormal
);
// f64 threshold 2048: past every f64 break point — exp-family overflow
// (exp2 at 1024), underflow (exp2 at -1075), and `tanh`/`expm1` saturation
// (~746). exp f64 overflows only at ~709.78, so e.g. exp(700) stays on the
// kernel path (it is finite) — the bug a too-low uniform threshold caused.
impl_sat_hw!(
    f64,
    2048.0,
    f64::from_bits(0x0000_0000_0000_0001) // 5e-324, smallest positive subnormal
);

/// `exp(x)`: overflow above the threshold, underflow to `+0` below it.
#[inline]
pub(crate) fn sat_exp<F: SatHw>(x: F, mode: RoundingMode) -> Option<(F, Status)> {
    if !x.is_finite() {
        return None;
    }
    let v = x.as_f64();
    if v >= F::THRESHOLD {
        Some(F::over_pos(mode))
    } else if v <= -F::THRESHOLD {
        Some(F::under_pos(mode))
    } else {
        None
    }
}

/// `exp2(x)` and `exp10(x)` saturate exactly like `exp`.
#[inline]
pub(crate) fn sat_exp2<F: SatHw>(x: F, mode: RoundingMode) -> Option<(F, Status)> {
    sat_exp(x, mode)
}

#[inline]
pub(crate) fn sat_exp10<F: SatHw>(x: F, mode: RoundingMode) -> Option<(F, Status)> {
    sat_exp(x, mode)
}

/// `expm1(x) = e^x - 1`: overflow above, but saturates to `-1` below
/// (`e^x - 1 -> -1` from above, never tiny).
#[inline]
pub(crate) fn sat_expm1<F: SatHw>(x: F, mode: RoundingMode) -> Option<(F, Status)> {
    if !x.is_finite() {
        return None;
    }
    let v = x.as_f64();
    if v >= F::THRESHOLD {
        Some(F::over_pos(mode))
    } else if v <= -F::THRESHOLD {
        Some(F::neg_one())
    } else {
        None
    }
}

/// `sinh(x)`: overflows to `+∞` above, `-∞` below (odd).
#[inline]
pub(crate) fn sat_sinh<F: SatHw>(x: F, mode: RoundingMode) -> Option<(F, Status)> {
    if !x.is_finite() {
        return None;
    }
    let v = x.as_f64();
    if v >= F::THRESHOLD {
        Some(F::over_pos(mode))
    } else if v <= -F::THRESHOLD {
        Some(F::over_neg(mode))
    } else {
        None
    }
}

/// `cosh(x)`: overflows to `+∞` for large `|x|` (even, always positive).
#[inline]
pub(crate) fn sat_cosh<F: SatHw>(x: F, mode: RoundingMode) -> Option<(F, Status)> {
    if !x.is_finite() {
        return None;
    }
    if x.as_f64().abs() >= F::THRESHOLD {
        Some(F::over_pos(mode))
    } else {
        None
    }
}

/// `tanh(x)`: saturates to `+1` above, `-1` below.
#[inline]
pub(crate) fn sat_tanh<F: SatHw>(x: F, _mode: RoundingMode) -> Option<(F, Status)> {
    if !x.is_finite() {
        return None;
    }
    let v = x.as_f64();
    if v >= F::THRESHOLD {
        Some(F::pos_one())
    } else if v <= -F::THRESHOLD {
        Some(F::neg_one())
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::round::{drive, F32Shell, F64Shell};
    use pfloat::BigFloat;

    /// A saturation fast-path function pointer for one width.
    type SatFn<H> = fn(H, RoundingMode) -> Option<(H, Status)>;
    /// The kernel closure `drive` consumes for one function.
    type KernelFn<'a> = dyn Fn(&BigFloat, u32, RoundingMode) -> (BigFloat, Status) + 'a;

    const MODES: [RoundingMode; 5] = [
        RoundingMode::NearestEven,
        RoundingMode::NearestAway,
        RoundingMode::TowardZero,
        RoundingMode::TowardPositive,
        RoundingMode::TowardNegative,
    ];

    // Kernel-safe huge `f32` inputs (>= the 1024 threshold), spanning the
    // boundary to 2^113. Capped at 2^113 so pfloat's i128 internal
    // exponent does not overflow when the kernel evaluates them for the
    // comparison (`tanh(2^113)` forms `e^(2^114)`, exponent ~3e34, still
    // inside i128); the fast-path itself handles up to `f32::MAX`.
    const HUGE_F32: &[f32] = &[
        1024.0,
        1025.0,
        2048.0,
        1.0e6,
        1.0e15,
        f32::from_bits(0x5d80_0000), // 2^60
        f32::from_bits(0x6c80_0000), // 2^90
        f32::from_bits(0x7800_0000), // 2^113
    ];
    // Kernel-safe huge `f64` inputs (>= the 2048 threshold), capped at
    // 2^120 (`exp(2^120)` exponent ~1.9e36, inside i128).
    const HUGE_F64: &[f64] = &[
        2048.0,
        2049.0,
        4096.0,
        1.0e6,
        1.0e30,
        f64::from_bits(0x4770_0000_0000_0000), // 2^120
    ];

    /// For every `x >= threshold` (and its negation) under every mode, the
    /// saturation fast-path must produce the bit-for-bit identical
    /// `(value, Status)` the raw kernel produces. The kernel is the
    /// oracle; matching it for the inputs the fast-path claims is the
    /// proof the saturation values and flags are correct.
    fn check_f32(sat: SatFn<f32>, kernel: &KernelFn, label: &str) {
        for &mag in HUGE_F32 {
            for &x in &[mag, -mag] {
                for &mode in &MODES {
                    let (sv, ss) =
                        sat(x, mode).unwrap_or_else(|| panic!("{label}: no fast-path at x={x}"));
                    let (kv, ks) = drive::<F32Shell>(x, mode, kernel);
                    assert_eq!(
                        sv.to_bits(),
                        kv.to_bits(),
                        "{label} VALUE x={x} {mode:?}: fast={:#010x} kernel={:#010x}",
                        sv.to_bits(),
                        kv.to_bits()
                    );
                    assert_eq!(ss, ks, "{label} STATUS x={x} {mode:?}");
                }
            }
        }
    }

    fn check_f64(sat: SatFn<f64>, kernel: &KernelFn, label: &str) {
        for &mag in HUGE_F64 {
            for &x in &[mag, -mag] {
                for &mode in &MODES {
                    let (sv, ss) =
                        sat(x, mode).unwrap_or_else(|| panic!("{label}: no fast-path at x={x}"));
                    let (kv, ks) = drive::<F64Shell>(x, mode, kernel);
                    assert_eq!(
                        sv.to_bits(),
                        kv.to_bits(),
                        "{label} VALUE x={x} {mode:?}: fast={:#018x} kernel={:#018x}",
                        sv.to_bits(),
                        kv.to_bits()
                    );
                    assert_eq!(ss, ks, "{label} STATUS x={x} {mode:?}");
                }
            }
        }
    }

    #[test]
    fn fast_path_matches_kernel_f32() {
        check_f32(sat_exp::<f32>, &|xb, w, d| xb.exp_round(w, d), "exp");
        check_f32(
            sat_exp2::<f32>,
            &|xb, w, d| xb.exp2_round(w, d).expect("w >= 1"),
            "exp2",
        );
        check_f32(
            sat_exp10::<f32>,
            &|xb, w, d| xb.exp10_round(w, d).expect("w >= 1"),
            "exp10",
        );
        check_f32(
            sat_expm1::<f32>,
            &|xb, w, d| xb.expm1_round(w, d).expect("w >= 1"),
            "expm1",
        );
        check_f32(
            sat_sinh::<f32>,
            &|xb, w, d| xb.sinh_round(w, d).expect("w >= 1"),
            "sinh",
        );
        check_f32(
            sat_cosh::<f32>,
            &|xb, w, d| xb.cosh_round(w, d).expect("w >= 1"),
            "cosh",
        );
        check_f32(
            sat_tanh::<f32>,
            &|xb, w, d| xb.tanh_round(w, d).expect("w >= 1"),
            "tanh",
        );
    }

    #[test]
    fn fast_path_matches_kernel_f64() {
        check_f64(sat_exp::<f64>, &|xb, w, d| xb.exp_round(w, d), "exp");
        check_f64(
            sat_exp2::<f64>,
            &|xb, w, d| xb.exp2_round(w, d).expect("w >= 1"),
            "exp2",
        );
        check_f64(
            sat_exp10::<f64>,
            &|xb, w, d| xb.exp10_round(w, d).expect("w >= 1"),
            "exp10",
        );
        check_f64(
            sat_expm1::<f64>,
            &|xb, w, d| xb.expm1_round(w, d).expect("w >= 1"),
            "expm1",
        );
        check_f64(
            sat_sinh::<f64>,
            &|xb, w, d| xb.sinh_round(w, d).expect("w >= 1"),
            "sinh",
        );
        check_f64(
            sat_cosh::<f64>,
            &|xb, w, d| xb.cosh_round(w, d).expect("w >= 1"),
            "cosh",
        );
        check_f64(
            sat_tanh::<f64>,
            &|xb, w, d| xb.tanh_round(w, d).expect("w >= 1"),
            "tanh",
        );
    }

    /// The fast-path must stay inactive below the threshold and on
    /// non-finite inputs, so the kernel path (with its exact `exp(∞)=∞`,
    /// `tanh(∞)=1`, NaN propagation) is preserved there.
    #[test]
    fn fast_path_inactive_below_threshold_and_nonfinite() {
        let fns32: &[SatFn<f32>] = &[
            sat_exp::<f32>,
            sat_exp2::<f32>,
            sat_exp10::<f32>,
            sat_expm1::<f32>,
            sat_sinh::<f32>,
            sat_cosh::<f32>,
            sat_tanh::<f32>,
        ];
        for f in fns32 {
            for &x in &[
                0.0f32,
                1.0,
                500.0,
                1023.0,
                -1023.0,
                f32::NAN,
                f32::INFINITY,
                f32::NEG_INFINITY,
            ] {
                assert!(
                    f(x, RoundingMode::NearestEven).is_none(),
                    "f32 fired at x={x}"
                );
            }
        }
        let fns64: &[SatFn<f64>] = &[sat_exp::<f64>, sat_tanh::<f64>, sat_cosh::<f64>];
        for f in fns64 {
            for &x in &[
                0.0f64,
                1.0,
                700.0,
                2047.0,
                -2047.0,
                f64::NAN,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ] {
                assert!(
                    f(x, RoundingMode::NearestEven).is_none(),
                    "f64 fired at x={x}"
                );
            }
        }
    }
}
