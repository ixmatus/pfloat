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
//!
//! Slice 4a opens Phase 4 (tier-1 special functions) with `erf`
//! and `erfc`. For `|x|` inside a cancellation-aware threshold
//! the kernel evaluates the Maclaurin series of `erf` at boosted
//! working precision; for larger `|x|` it uses the asymptotic
//! expansion of `erfc` and complements when erf is requested. A
//! hardcoded 1024-bit `2/sqrt(π)` constant supplies the leading
//! coefficient.

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
#[cfg(feature = "exp-log")]
pub(crate) mod ziv;

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

#[cfg(feature = "agm")]
pub(crate) mod agm;

#[cfg(feature = "exp-log")]
pub(crate) mod agm_constants;

#[cfg(feature = "specials")]
pub(crate) mod beta;
#[cfg(feature = "specials")]
pub(crate) mod digamma;
#[cfg(feature = "specials")]
pub(crate) mod erf;
#[cfg(feature = "specials")]
pub(crate) mod erfc;
#[cfg(feature = "specials")]
pub(crate) mod gamma;
#[cfg(feature = "specials")]
pub(crate) mod gamma_stirling;
#[cfg(feature = "specials")]
pub(crate) mod lgamma;

#[cfg(feature = "integrals")]
pub(crate) mod ci;
#[cfg(feature = "integrals")]
pub(crate) mod ei;
#[cfg(feature = "integrals")]
pub(crate) mod li;
#[cfg(feature = "integrals")]
pub(crate) mod si;

#[cfg(feature = "airy")]
pub(crate) mod airy;

#[cfg(feature = "bessel")]
pub(crate) mod bessel_j;

#[cfg(feature = "bessel")]
pub(crate) mod bessel_y;

#[cfg(feature = "bessel")]
pub(crate) mod bessel_i;

#[cfg(feature = "bessel")]
pub(crate) mod bessel_k;

#[cfg(feature = "zeta")]
pub(crate) mod zeta;

/// Hardcoded `ln(2)` mantissa at 1024-bit precision.
///
/// Layout: little-endian limbs, top-bit-set. The mantissa-as-integer
/// equals `floor(ln(2) × 2^1024)`. Combined with `precision = 1024`
/// and `exponent = -1`, this represents `ln(2)`.
///
/// Source: the correctly-rounded 1024-bit value of `ln(2)`. Slice
/// 7b2 regenerated these limbs after slice 7b's diagnostic found the
/// original encoding faithful only through the first ~450 bits (the
/// low nine limbs were wrong; the seven most-significant limbs here
/// are unchanged from the original, matching that boundary). The
/// limbs are the mantissa of `BigFloat::parse_str` applied to the
/// authoritative 1100-digit decimal in the `agm_constants` test
/// module's `LN2_REFERENCE`, independently cross-checked bit-for-bit
/// against the in-repo AGM atanh-series derivation rounded down from
/// 2048 bits. They are a mathematical fact derived from primary
/// sources, not adapted from any implementation's source. The
/// `ln_2_table_is_correctly_rounded_at_p1024` regression test pins
/// both derivations so any future drift fails CI.
#[allow(dead_code)] // referenced by exp and ln; treat as logically pub(super) for the cluster.
pub(crate) const LN2_LIMBS_1024: [u64; 16] = [
    0xDA2D_97C5_0F3F_D5C6,
    0x655F_A187_2F20_E3A2,
    0xF5DF_A6BD_3830_3248,
    0x72CE_87B1_9D65_48CA,
    0x256F_A0EC_7657_F74B,
    0xB9EA_9BC3_B136_603B,
    0x1ACB_DA11_317C_387E,
    0x3E96_CA16_224A_E8C5,
    0x2757_3B29_1169_B825,
    0xED2E_AE35_C138_2144,
    0x5595_52FB_4AFA_1B10,
    0xE7B8_7620_6DEB_AC98,
    0x8A0D_175B_8BAA_FA2B,
    0x40F3_4326_7298_B62D,
    0xC9E3_B398_03F2_F6AF,
    0xB172_17F7_D1CF_79AB,
];

/// Threshold below which [`ln_2_at`] returns the rounded hardcoded
/// `LN2_LIMBS_1024` constant. Slice 7b capped this at 448 because the
/// original table was faithful only through ~450 bits. Slice 7b2
/// regenerated the table to the correctly-rounded 1024-bit value, so
/// the cap is restored to the full table width, mirroring `pi_at`'s
/// `prec <= 1024` fast path. At `prec = 1024` the rounded table value
/// and the AGM atanh series are bit-identical (both equal the
/// correctly-rounded 1024-bit `ln(2)`), so the dispatch boundary is
/// seamless; the `ln_2_table_matches_agm_at_cap_precision` regression
/// test pins that equality. Precisions above 1024 route through the
/// AGM series, which is correct at any precision. ADR-0017 records
/// the design and the slice 7b2 regeneration.
#[cfg(feature = "exp-log")]
pub(crate) const LN2_TABLE_PRECISION_CAP: u32 = 1024;

/// Returns `ln(2)` rounded to the requested precision.
///
/// For `prec <= LN2_TABLE_PRECISION_CAP` (1024 post slice 7b2) the
/// rounded value comes from the correctly-rounded hardcoded table.
/// For larger precisions the value is computed on the fly via the
/// AGM-based atanh series in [`agm_constants::ln_2_via_atanh`].
#[cfg(feature = "exp-log")]
#[allow(dead_code)]
pub(crate) fn ln_2_at(prec: u32) -> BigFloat {
    if prec <= LN2_TABLE_PRECISION_CAP {
        ln_2_via_table(prec)
    } else {
        agm_constants::ln_2_via_atanh(prec)
    }
}

/// Returns the rounded hardcoded `ln(2)` constant. Internal: the
/// public dispatcher [`ln_2_at`] picks this for `prec <= 1024`.
#[cfg(feature = "exp-log")]
#[allow(dead_code)]
pub(crate) fn ln_2_via_table(prec: u32) -> BigFloat {
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

/// Hardcoded Euler–Mascheroni `γ` mantissa at 1024-bit precision.
///
/// Layout: little-endian limbs, top-bit-set. With `precision = 1024`,
/// `exponent = -1`, and a positive sign this represents
/// `γ ≈ 0.5772156649…`.
///
/// Source: the correctly-rounded 1024-bit value of γ (slice 6m0).
/// The limbs are the mantissa of `BigFloat::parse_str` applied to the
/// authoritative OEIS A001620 decimal (`EULER_GAMMA_REFERENCE` in the
/// `agm_constants` test module), independently cross-checked
/// bit-for-bit against the in-repo Brent–McMillan derivation rounded
/// down from 2048 bits. They are a mathematical fact derived from
/// primary sources, not adapted from any implementation. The
/// `euler_gamma_table_is_correctly_rounded_at_p1024` regression test
/// pins both derivations. ADR-0018 records the design.
#[cfg(feature = "exp-log")]
#[allow(dead_code)] // referenced by the integral specials cluster
pub(crate) const EULER_GAMMA_LIMBS_1024: [u64; 16] = [
    0x9615_42A3_CE3B_EA5E,
    0x5E6A_C2F0_BD61_C746,
    0x3EC7_C271_8279_7722,
    0xD2A1_EA1D_E62F_F864,
    0x0C09_D4C8_B6B7_B86F,
    0x8A96_D156_7899_AAAE,
    0xDBE7_BF38_154B_36CF,
    0x58DE_B878_CC86_D733,
    0xE43B_4673_D74B_AFEA,
    0x1056_AE91_3213_5A08,
    0xD064_9CCB_6210_57D1,
    0x8E4B_59FA_03A9_F0EE,
    0x0C03_DF34_709A_FFBD,
    0xA1CE_CC3A_F65C_C019,
    0xD1BE_3F81_0152_CB56,
    0x93C4_67E3_7DB0_C7A4,
];

/// Returns the Euler–Mascheroni constant `γ` rounded to the
/// requested precision.
///
/// For `prec <= 1024` the rounded value comes from the
/// correctly-rounded hardcoded [`EULER_GAMMA_LIMBS_1024`] table; for
/// larger precisions it is computed on the fly via the Brent–McMillan
/// algorithm in [`agm_constants::euler_gamma_via_bm`]. There is no
/// precision cap (γ is correct at any precision), mirroring `ln_2_at`
/// and `pi_at`.
#[cfg(feature = "exp-log")]
#[allow(dead_code)] // referenced by the integral specials cluster
pub(crate) fn euler_gamma_at(prec: u32) -> BigFloat {
    if prec <= 1024 {
        euler_gamma_via_table(prec)
    } else {
        agm_constants::euler_gamma_via_bm(prec)
    }
}

/// Returns the rounded hardcoded `γ` constant. Internal: the public
/// dispatcher [`euler_gamma_at`] picks this for `prec <= 1024`.
#[cfg(feature = "exp-log")]
#[allow(dead_code)]
pub(crate) fn euler_gamma_via_table(prec: u32) -> BigFloat {
    let stored = BigFloat {
        class: Class::Normal {
            sign: Sign::Positive,
            exponent: -1,
            mantissa: EULER_GAMMA_LIMBS_1024.to_vec(),
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
/// Computed on the fly via
/// [`agm_constants::ln_10_via_atanh`] (the identity
/// `ln(10) = 3·ln(2) + 2·atanh(1/9)`). The atanh(1/9) series
/// converges roughly `log₂(81) ≈ 6.34` bits per term, faster than
/// routing through the full `BigFloat::ln` kernel which would
/// duplicate the same factorization.
#[cfg(feature = "exp-log")]
#[allow(dead_code)]
pub(crate) fn ln_10_at(prec: u32) -> BigFloat {
    agm_constants::ln_10_via_atanh(prec)
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

/// Returns `π` rounded to the requested precision.
///
/// For `prec <= 1024` the rounded value comes from
/// [`PI_LIMBS_1024`]. For `prec > 1024` slice 7b computes the
/// value via the Brent–Salamin iteration in
/// [`agm_constants::pi_via_agm`].
#[cfg(feature = "trig")]
#[allow(dead_code)]
pub(crate) fn pi_at(prec: u32) -> BigFloat {
    if prec <= 1024 {
        pi_via_table(prec)
    } else {
        agm_constants::pi_via_agm(prec)
    }
}

/// Returns the rounded hardcoded `π` constant. Internal: the public
/// dispatcher [`pi_at`] picks this for `prec <= 1024`.
#[cfg(feature = "trig")]
#[allow(dead_code)]
pub(crate) fn pi_via_table(prec: u32) -> BigFloat {
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

/// Returns `π/2` rounded to the requested precision. Composes
/// `pi_at(prec)` with an exponent shift for the `/2`. The shift is
/// exact (multiplication by `2^-1`).
#[cfg(feature = "trig")]
#[allow(dead_code)]
pub(crate) fn pi_over_2_at(prec: u32) -> BigFloat {
    if prec <= 1024 {
        let stored = BigFloat {
            class: Class::Normal {
                sign: Sign::Positive,
                exponent: 0,
                mantissa: PI_LIMBS_1024.to_vec(),
            },
            precision: 1024,
        };
        return stored
            .round_to_precision(prec, RoundingMode::NearestEven)
            .expect("precision >= 1")
            .0;
    }
    let pi = agm_constants::pi_via_agm(prec);
    // Exact /2 via an exponent decrement on the Normal variant.
    match pi.class {
        Class::Normal {
            sign,
            exponent,
            mantissa,
        } => BigFloat {
            class: Class::Normal {
                sign,
                exponent: exponent - 1,
                mantissa,
            },
            precision: pi.precision,
        },
        _ => unreachable!("pi_via_agm returns a Normal"),
    }
}

/// Returns `2/π` rounded to the requested precision.
///
/// For `prec <= 4096` the rounded value comes from the 4096-bit
/// [`TWO_OVER_PI_LIMBS_4096`] table that drives the Payne-Hanek
/// trig argument reduction. For `prec > 4096` slice 7b computes
/// the value via [`agm_constants::two_over_pi_via_agm`].
#[cfg(feature = "trig")]
#[allow(dead_code)]
pub(crate) fn two_over_pi_at(prec: u32) -> BigFloat {
    if prec <= 4096 {
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
    } else {
        agm_constants::two_over_pi_via_agm(prec)
    }
}

/// Hardcoded `2/sqrt(π)` mantissa at 1024-bit precision.
///
/// The correctly-rounded (round-to-nearest-even) 1024-bit value of
/// `2/sqrt(π)` as 16 little-endian u64 limbs. Combined with
/// `precision = 1024` and `exponent = 0`, this represents
/// `2/sqrt(π) ≈ 1.1283791670955125…`, the leading coefficient of the
/// Maclaurin series for `erf`.
///
/// Slice-4a generated this by truncation, leaving the
/// least-significant limb 1 ULP low (`…6C03` for the correctly-rounded
/// `…6C04`); the value was returned 1 ULP below correctly-rounded for
/// precisions near 1024. Slice 6m-audit regenerated it from the
/// authoritative `mpmath` decimal via the bit-exact decimal parser,
/// cross-checked bit-for-bit against the independent in-repo
/// `2 / sqrt(π)` AGM derivation at 2048 bits rounded to 1024.
/// `two_over_sqrt_pi_table_is_correctly_rounded_at_p1024` pins this
/// three-way equality so the defect cannot silently recur.
#[cfg(feature = "specials")]
#[allow(dead_code)]
pub(crate) const TWO_OVER_SQRT_PI_LIMBS_1024: [u64; 16] = [
    0x18D3E91ADCFF6C04,
    0x50754B409E94D32D,
    0xAC2C88BBBA81B1C7,
    0xEB9FEB2436F2F272,
    0xD27A3282DADA7316,
    0x9522F2F93E16B2A3,
    0x9C22F47F7B7FB57C,
    0x52561DCC244DC65E,
    0x74F76F877FFEC251,
    0xBD1F4EEE48E1CA78,
    0x40C036096CC79AEB,
    0xC0759CF859270F11,
    0x39A15830CCE620B0,
    0x1409A0EBAC3E7517,
    0x71D48A7F6BFEC344,
    0x906EBA8214DB688D,
];

/// Returns `2/sqrt(π)` rounded to the requested precision.
///
/// For `prec <= 1024` the rounded value comes from
/// [`TWO_OVER_SQRT_PI_LIMBS_1024`]. For `prec > 1024` slice 7b
/// computes the value via
/// [`agm_constants::two_over_sqrt_pi_via_agm`].
#[cfg(feature = "specials")]
#[allow(dead_code)]
pub(crate) fn two_over_sqrt_pi_at(prec: u32) -> BigFloat {
    if prec <= 1024 {
        two_over_sqrt_pi_via_table(prec)
    } else {
        agm_constants::two_over_sqrt_pi_via_agm(prec)
    }
}

/// Returns the rounded hardcoded `2/√π` constant. Internal: the public
/// dispatcher [`two_over_sqrt_pi_at`] picks this for `prec <= 1024`.
#[cfg(feature = "specials")]
#[allow(dead_code)]
pub(crate) fn two_over_sqrt_pi_via_table(prec: u32) -> BigFloat {
    let stored = BigFloat {
        class: Class::Normal {
            sign: Sign::Positive,
            exponent: 0,
            mantissa: TWO_OVER_SQRT_PI_LIMBS_1024.to_vec(),
        },
        precision: 1024,
    };
    stored
        .round_to_precision(prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}

/// Hardcoded `ln(2π)` mantissa at 1024-bit precision.
///
/// `floor(ln(2π) · 2^1023)` as 16 little-endian u64 limbs. Combined
/// with `precision = 1024` and `exponent = 0`, this represents
/// `ln(2π) ≈ 1.8378770664093454836…`, the constant term in
/// Stirling's asymptotic series for `ln Γ(z)`.
///
/// Source: `mpmath` (Python) at 5000-bit working precision,
/// generated at slice-4b authoring time.
#[cfg(feature = "specials")]
#[allow(dead_code)]
pub(crate) const LN_2PI_LIMBS_1024: [u64; 16] = [
    0x3BD6E6FBA48AA194,
    0x54B6D36BEE63E04A,
    0xC525605F70BB125E,
    0xDED77FBEC954A0AF,
    0x27086C366978E17E,
    0x9254D1304A59FB7E,
    0x307D867635C11696,
    0x926770ECA54487A7,
    0xCF66ECE1772BADF2,
    0xB05CAB571B4CDA5B,
    0x93EABF905C5569BB,
    0x212F9D7FE00E86BF,
    0xDEC6A3133DAA155D,
    0xCFB08F8D13458B4D,
    0x94BC900144192023,
    0xEB3F8E4325F5A534,
];

/// Returns `ln(2π)` rounded to the requested precision.
///
/// For `prec <= 1024` the rounded value comes from
/// [`LN_2PI_LIMBS_1024`]. For `prec > 1024` slice 7b computes
/// the value via [`agm_constants::ln_2pi_via_agm`] (the identity
/// `ln(2π) = ln(2) + ln(π)`, with the atanh series carrying both
/// terms).
#[cfg(feature = "specials")]
#[allow(dead_code)]
pub(crate) fn ln_2pi_at(prec: u32) -> BigFloat {
    if prec <= 1024 {
        ln_2pi_via_table(prec)
    } else {
        agm_constants::ln_2pi_via_agm(prec)
    }
}

/// Returns the rounded hardcoded `ln(2π)` constant. Internal: the
/// public dispatcher [`ln_2pi_at`] picks this for `prec <= 1024`.
#[cfg(feature = "specials")]
#[allow(dead_code)]
pub(crate) fn ln_2pi_via_table(prec: u32) -> BigFloat {
    let stored = BigFloat {
        class: Class::Normal {
            sign: Sign::Positive,
            exponent: 0,
            mantissa: LN_2PI_LIMBS_1024.to_vec(),
        },
        precision: 1024,
    };
    stored
        .round_to_precision(prec, RoundingMode::NearestEven)
        .expect("precision >= 1")
        .0
}
