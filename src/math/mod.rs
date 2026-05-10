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
//!
//! Slice 3f opens the trig family with the forward functions
//! `sin`, `cos`, `tan`. The argument reduction uses a Payne-Hanek-
//! style multiplication of `x` by a hardcoded 4096-bit table of
//! `2/π`: the result splits into an integer `q` (the quadrant
//! index, taken mod 4) and a normalized fractional remainder `f`
//! such that `x = (q + f) · π/2`. The reduced argument
//! `r = f · π/2 ∈ [−π/4, π/4]` then drives the Taylor series.
//! The table size caps the supported input range at roughly
//! `|x| < 2^3000`; beyond that, reduction loses precision and the
//! kernel raises `INVALID` with a quiet NaN result.
//!
//! Slice 3g closes the trig family with the inverse functions:
//! `atan`, `asin`, `acos`, and `atan2`. `atan` is the core
//! primitive — it reduces large `|x|` via the reciprocal identity
//! and small `|x|` via the half-angle identity, then sums the
//! Taylor series at ~4 bits per term. The other three route
//! through `atan` via cancellation-free identities so that the
//! `±1` boundaries of `asin`/`acos` produce exact `±π/2` and
//! `0` / `π` results. `atan2` dispatches the full IEEE 754-2019
//! §9.2.1 special-case table for `(y, x)` zero, infinity, and
//! quadrant pairs before delegating to `atan(y/x)` with a `±π`
//! shift for negative `x`.

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

#[cfg(feature = "trig")]
pub(crate) mod acos;
#[cfg(feature = "trig")]
pub(crate) mod asin;
#[cfg(feature = "trig")]
pub(crate) mod atan;
#[cfg(feature = "trig")]
pub(crate) mod atan2;
#[cfg(feature = "trig")]
pub(crate) mod cos;
#[cfg(feature = "trig")]
pub(crate) mod sin;
#[cfg(feature = "trig")]
pub(crate) mod tan;
#[cfg(feature = "trig")]
pub(crate) mod trig_reduce;

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

/// Hardcoded `π` mantissa at 1024-bit precision.
///
/// `floor(π · 2^1022)` as 16 little-endian u64 limbs. Combined with
/// `precision = 1024` and `exponent = 1`, this represents `π` (the
/// value lives in `[2, 4)`, so the top bit sits at position 1).
///
/// Source: `mpmath` (Python) at 5000-bit working precision,
/// generated at slice-3f authoring time. Verified by the inline
/// `pi_constants_self_consistent` test in
/// [`crate::math::trig_reduce`] which checks `π · (2/π) ≈ 2`.
#[cfg(feature = "trig")]
#[allow(dead_code)]
pub(crate) const PI_LIMBS_1024: [u64; 16] = [
    0x98DA48361C55D39A,
    0xC2007CB8A163BF05,
    0x49286651ECE45B3D,
    0xAE9F24117C4B1FE6,
    0xEE386BFB5A899FA5,
    0x0BFF5CB6F406B7ED,
    0xF44C42E9A637ED6B,
    0xE485B576625E7EC6,
    0x4FE1356D6D51C245,
    0x302B0A6DF25F1437,
    0xEF9519B3CD3A431B,
    0x514A08798E3404DD,
    0x020BBEA63B139B22,
    0x29024E088A67CC74,
    0xC4C6628B80DC1CD1,
    0xC90FDAA22168C234,
];

/// Hardcoded `2/π` mantissa at 4096-bit precision, the Payne-Hanek
/// reduction table.
///
/// `floor((2/π) · 2^4096)` as 64 little-endian u64 limbs. Combined
/// with `precision = 4096` and `exponent = −1`, this represents
/// `2/π` (the value lives in `[1/2, 1)`).
///
/// The table size caps the supported input range: an input `x` with
/// binary exponent `e_x` is correctly reducible while
/// `e_x + working_precision + slack < 4096`. For the default
/// working precision and reasonable target precisions (≤ 1024), the
/// kernel handles `|x| < 2^3000` cleanly; beyond that, the
/// reduction loses precision and the trig kernels report
/// `INVALID` with a quiet NaN result.
///
/// Source: `mpmath` (Python) at 5000-bit working precision,
/// generated at slice-3f authoring time.
#[cfg(feature = "trig")]
#[allow(dead_code)]
pub(crate) const TWO_OVER_PI_LIMBS_4096: [u64; 64] = [
    0x36D9CAD2A8288D61,
    0x818D67C12645CA55,
    0x6F63A62DCBBFF4EF,
    0x78738A5A8CAFBDD7,
    0x775C83C2A3883C61,
    0x0AB499D3F2A6067F,
    0x425FAECE616AA428,
    0x4A48D36710D8DDAA,
    0xF57FB0ADF2E91E43,
    0x6212830148835B8E,
    0x1DF35BE01834132E,
    0x08CB7DE050C017A7,
    0x4D58E232CAC616E3,
    0x9BDE2822D2E88628,
    0x5DD7DE16DE3B5892,
    0xCDC4EF09366CD43F,
    0x652289E83260BFE6,
    0x9947FBACD87F7EB7,
    0xFF319F6A1E666157,
    0x1F001B0AF1DFCE19,
    0x24778AD623545AB9,
    0xD9D63B3884A7CB23,
    0xB07AE715175649C0,
    0x64ABD770F87C6357,
    0x1810A3FC764D2A9D,
    0xA7B4D55537F63ED7,
    0x9B0062337CD2B497,
    0x467D862D71E39AC6,
    0xC4AD414D2C5D000C,
    0x15C614B59D19C3C2,
    0xFA6ED5772D30433B,
    0x87F121907C7C246A,
    0x9F3A1F35CAF27F1D,
    0xC33D26EF6B1E5EF8,
    0x32C2DE4F98327DBB,
    0xA5FF07053F7E33E8,
    0xDDAF44D15719053E,
    0x8359C4768B961CA6,
    0x19C367CDDCE8092A,
    0x60E27BC08C6B47C4,
    0x06061556CA73A8C9,
    0x8DFFD8804D732731,
    0x6599855F14A06840,
    0xA9E391615EE61B08,
    0xF0CFBC209AF4361D,
    0x56033046FC7B6BAB,
    0x6BFB5FB11F8D5D08,
    0x3D0739F78A5292EA,
    0x7527BAC7EBE5F17B,
    0x4F463F669E5FEA2D,
    0x6D367ECF27CB09B7,
    0xEF2F118B5A0A6D1F,
    0x1FF897FFDE05980F,
    0x9C845F8BBDF9283B,
    0x3991D639835339F4,
    0xE99C7026B45F7E41,
    0xE88235F52EBB4484,
    0xFE1DEB1CB129A73E,
    0x06492EEA09D1921C,
    0xB7246E3A424DD2E0,
    0xFE5163ABDEBBC561,
    0xDB6295993C439041,
    0xFC2757D1F534DDC0,
    0xA2F9836E4E441529,
];

/// Returns `π` rounded to the requested precision (up to 1024 bits
/// faithfully).
#[cfg(feature = "trig")]
#[allow(dead_code)]
pub(crate) fn pi_at(prec: u32) -> BigFloat {
    let stored = BigFloat {
        class: Class::Normal {
            sign: Sign::Positive,
            exponent: 1,
            mantissa: PI_LIMBS_1024.to_vec(),
        },
        precision: 1024,
    };
    stored
        .round_to_precision(prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}

/// Returns `π/2` rounded to the requested precision (up to 1024
/// bits faithfully). Constructed from [`PI_LIMBS_1024`] with a
/// single decremented exponent so no division is required.
#[cfg(feature = "trig")]
#[allow(dead_code)]
pub(crate) fn pi_over_2_at(prec: u32) -> BigFloat {
    let stored = BigFloat {
        class: Class::Normal {
            sign: Sign::Positive,
            exponent: 0,
            mantissa: PI_LIMBS_1024.to_vec(),
        },
        precision: 1024,
    };
    stored
        .round_to_precision(prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}

/// Returns `2/π` rounded to the requested precision (up to 4096
/// bits faithfully).
#[cfg(feature = "trig")]
#[allow(dead_code)]
pub(crate) fn two_over_pi_at(prec: u32) -> BigFloat {
    let stored = BigFloat {
        class: Class::Normal {
            sign: Sign::Positive,
            exponent: -1,
            mantissa: TWO_OVER_PI_LIMBS_4096.to_vec(),
        },
        precision: 4096,
    };
    stored
        .round_to_precision(prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}
