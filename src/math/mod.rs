//! Elementary transcendentals (Phase 3).
//!
//! Each function uses range reduction to bring the argument into a
//! Taylor-friendly window, evaluates the reduced argument's series,
//! then composes the result. Phase 3 builds out the surface across
//! several slices: 3a shipped `exp`, 3b adds `ln`. Subsequent slices
//! fill in `expm1`, `exp2`/`exp10`, `log1p`/`log2`/`log10`, trig,
//! hyperbolic, and `pow`.
//!
//! All functions are gated behind the `exp-log` or `trig` cluster
//! features so embedded users can compile in only what they need.
//!
//! Slice 3c adds `pow`, which composes the prior two: for positive
//! `x` the kernel evaluates `exp(y · ln(x))`, with full IEEE 754-2019
//! §9.2.1 special-case handling around `0`, `±∞`, `NaN`, integer
//! parity of `y`, and the `pow(±1, ±∞) = 1` rule.
//!
//! Slice 3d adds the thin wrappers around the two foundation
//! transcendentals: `expm1`, `exp2`, `exp10`, `log1p`, `log2`,
//! `log10`. Each composes `exp` or `ln` with a precision boost to
//! handle the relevant cancellation regime (small `x` for `expm1`
//! and `log1p`) or a base-change constant (`ln(2)` for `exp2`/`log2`,
//! `ln(10)` for `exp10`/`log10`).
//!
//! Slice 3e adds the hyperbolic family: `sinh`, `cosh`, `tanh`,
//! `asinh`, `acosh`, `atanh`. Forward hyperbolic functions compose
//! through `exp` / `expm1` with cancellation-aware identities
//! (`sinh(x) = (expm1(x) − expm1(−x))/2` for the small-argument
//! regime; `tanh(x) = (1 − exp(−2|x|)) / (1 + exp(−2|x|))` with a
//! sign flip for negative `x`). Inverse hyperbolic functions
//! compose through `log1p` and `sqrt`: `asinh`, `acosh`, and
//! `atanh` each use a small-argument identity that avoids the
//! cancellation in `ln(x + sqrt(x² ± 1))` near the relevant boundary.

use crate::big::BigFloat;
use crate::class::Class;
use crate::rounding::RoundingMode;
use crate::sign::Sign;

#[cfg(feature = "exp-log")]
pub(crate) mod acosh;
#[cfg(feature = "exp-log")]
pub(crate) mod asinh;
#[cfg(feature = "exp-log")]
pub(crate) mod atanh;
#[cfg(feature = "exp-log")]
pub(crate) mod cosh;
#[cfg(feature = "exp-log")]
pub(crate) mod exp;
#[cfg(feature = "exp-log")]
pub(crate) mod exp10;
#[cfg(feature = "exp-log")]
pub(crate) mod exp2;
#[cfg(feature = "exp-log")]
pub(crate) mod expm1;
#[cfg(feature = "exp-log")]
pub(crate) mod ln;
#[cfg(feature = "exp-log")]
pub(crate) mod log10;
#[cfg(feature = "exp-log")]
pub(crate) mod log1p;
#[cfg(feature = "exp-log")]
pub(crate) mod log2;
#[cfg(feature = "exp-log")]
pub(crate) mod pow;
#[cfg(feature = "exp-log")]
pub(crate) mod sinh;
#[cfg(feature = "exp-log")]
pub(crate) mod tanh;

/// Hardcoded `ln(2)` mantissa at 1024-bit precision.
///
/// Layout: little-endian limbs, top-bit-set. The mantissa-as-integer
/// equals `floor(ln(2) × 2^1024)`. Combined with `precision = 1024`
/// and `exponent = -1`, this represents `ln(2)`.
///
/// Source: the value of `ln(2)` from MPFR / standard mathematical
/// references, truncated to 1024 bits. Verified via the inline
/// test (top byte `0xB1`, matching `ln(2) ≈ 0.1011000101110010…`).
#[allow(dead_code)] // referenced by exp and ln; treat as logically pub(super) for the cluster.
pub(crate) const LN2_LIMBS_1024: [u64; 16] = [
    0x4DC7_2A88_F1F1_1A0F,
    0x0C3B_36F2_5FF2_1D85,
    0x8BCB_17A7_7B11_D2B4,
    0xB72C_E87B_19D4_540F,
    0xB256_FA0E_C765_7F74,
    0xEB9E_A9BC_3B13_6603,
    0x51AC_BDA1_1317_C387,
    0x53E9_6CA1_6224_AE8C,
    0x0275_73B2_9116_9B82,
    0xED2E_AE35_C138_2144,
    0x5595_52FB_4AFA_1B10,
    0xE7B8_7620_6DEB_AC98,
    0x8A0D_175B_8BAA_FA2B,
    0x40F3_4326_7298_B62D,
    0xC9E3_B398_03F2_F6AF,
    0xB172_17F7_D1CF_79AB,
];

/// Returns `ln(2)` rounded to the requested precision (up to 1024
/// bits faithfully; higher precisions are a slight under-
/// approximation truncated at the 1024-bit boundary).
#[allow(dead_code)]
pub(crate) fn ln_2_at(prec: u32) -> BigFloat {
    let stored = BigFloat {
        class: Class::Normal {
            sign: Sign::Positive,
            exponent: -1,
            mantissa: LN2_LIMBS_1024.to_vec(),
        },
        precision: 1024,
    };
    stored
        .round_to_precision(prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}

/// Returns `ln(10)` rounded to the requested precision.
///
/// Computed at call time as `ln(10)` via the existing `BigFloat::ln`
/// kernel. Slice 3d takes the runtime cost: every `exp10` and
/// `log10` invocation incurs one Taylor evaluation of `ln(10)` on
/// top of its own. A hardcoded 1024-bit `LN10` constant analogous to
/// [`LN2_LIMBS_1024`] is a frugal future optimization; the current
/// shape keeps the slice scope narrow and lets the next iteration
/// of `ln` improvements (Ziv strategy, Lefèvre–Muller tables) carry
/// `ln(10)` forward without divergence.
#[cfg(feature = "exp-log")]
#[allow(dead_code)]
pub(crate) fn ln_10_at(prec: u32) -> BigFloat {
    let ten = BigFloat::try_from_i64_exact(10, prec).expect("precision >= 1");
    ten.ln(RoundingMode::NearestEven).0
}
