//! Lefevre-Muller corpus accessor: map a [`LibmFnId`] to its
//! hard-to-round binary64 cases.

#![cfg(all(unix, feature = "differential-mpfr"))]

use super::lefevre_muller_data as d;
use super::types::LibmFnId;

/// One corpus case: `(input_bits, NE_output_bits)` at binary64.
pub type Case = d::Case;

/// The Lefevre-Muller hard-to-round corpus for `f`, or empty for the
/// functions the corpus does not cover (`sqrt`, `cbrt`, `cot`, `sec`,
/// `csc`, `hypot`, `rootn`).
pub fn lm_seeds_for(f: LibmFnId) -> &'static [Case] {
    match f {
        LibmFnId::Exp => d::EXP_CASES,
        LibmFnId::Exp2 => d::EXP2_CASES,
        LibmFnId::Exp10 => d::EXP10_CASES,
        LibmFnId::Expm1 => d::EXPM1_CASES,
        LibmFnId::Ln => d::LN_CASES,
        LibmFnId::Log2 => d::LOG2_CASES,
        LibmFnId::Log10 => d::LOG10_CASES,
        LibmFnId::Log1p => d::LOG1P_CASES,
        LibmFnId::Sin => d::SIN_CASES,
        LibmFnId::Cos => d::COS_CASES,
        LibmFnId::Tan => d::TAN_CASES,
        LibmFnId::Asin => d::ASIN_CASES,
        LibmFnId::Acos => d::ACOS_CASES,
        LibmFnId::Atan => d::ATAN_CASES,
        LibmFnId::Sinh => d::SINH_CASES,
        LibmFnId::Cosh => d::COSH_CASES,
        LibmFnId::Tanh => d::TANH_CASES,
        LibmFnId::Asinh => d::ASINH_CASES,
        LibmFnId::Acosh => d::ACOSH_CASES,
        LibmFnId::Atanh => d::ATANH_CASES,
        _ => &[],
    }
}

/// The 20 functions the corpus covers, for sweep drivers that want to
/// iterate them.
pub const COVERED: &[LibmFnId] = &[
    LibmFnId::Exp,
    LibmFnId::Exp2,
    LibmFnId::Exp10,
    LibmFnId::Expm1,
    LibmFnId::Ln,
    LibmFnId::Log2,
    LibmFnId::Log10,
    LibmFnId::Log1p,
    LibmFnId::Sin,
    LibmFnId::Cos,
    LibmFnId::Tan,
    LibmFnId::Asin,
    LibmFnId::Acos,
    LibmFnId::Atan,
    LibmFnId::Sinh,
    LibmFnId::Cosh,
    LibmFnId::Tanh,
    LibmFnId::Asinh,
    LibmFnId::Acosh,
    LibmFnId::Atanh,
];
