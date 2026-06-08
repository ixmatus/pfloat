//! Shared helpers for the pfloat-complex C5 verification lanes (ADR-0092).
//!
//! The enumerated Annex G tables (`annex_g_special_values.rs`), the dispatch
//! totality enumeration (`dispatch_totality.rs`), the algebraic identities
//! (`identities.rs`), and the independent acb differential
//! (`differential_acb.rs`) all draw on the constructors, the IEEE class model,
//! and the deterministic PRNG defined here, so a tightening of a generator
//! improves every lane at once and a divergence localizes to the oracle rather
//! than to the inputs.

// Each lane that `mod common;`s this file uses a different subset; the unused
// subset would otherwise warn under that binary's compilation.
#![allow(dead_code)]

use core::cmp::Ordering;

use pfloat::{BigFloat, RoundingMode, Sign};
use pfloat_complex::Complex;

pub const NE: RoundingMode = RoundingMode::NearestEven;

/// The five IEEE 754 rounding modes, in a fixed order for mode sweeps.
pub const ALL_MODES: [RoundingMode; 5] = [
    RoundingMode::NearestEven,
    RoundingMode::NearestAway,
    RoundingMode::TowardZero,
    RoundingMode::TowardPositive,
    RoundingMode::TowardNegative,
];

// ---------- BigFloat constructors for the special-value classes ----------

/// An exact integer at precision `p`.
pub fn bf(n: i64, p: u32) -> BigFloat {
    BigFloat::try_from_i64_exact(n, p).expect("i64 fits, p >= 1")
}

/// `+0` at precision `p`.
pub fn pz(p: u32) -> BigFloat {
    BigFloat::try_new_zero(Sign::Positive, p).expect("p >= 1")
}

/// `-0` at precision `p`.
pub fn nz(p: u32) -> BigFloat {
    BigFloat::try_new_zero(Sign::Negative, p).expect("p >= 1")
}

/// `+inf` at precision `p`.
pub fn pinf(p: u32) -> BigFloat {
    BigFloat::try_new_infinity(Sign::Positive, p).expect("p >= 1")
}

/// `-inf` at precision `p`.
pub fn ninf(p: u32) -> BigFloat {
    BigFloat::try_new_infinity(Sign::Negative, p).expect("p >= 1")
}

/// A quiet NaN at precision `p`.
pub fn qnan(p: u32) -> BigFloat {
    BigFloat::try_new_quiet_nan(Sign::Positive, p, &[]).expect("p >= 1")
}

/// A signaling NaN at precision `p`.
pub fn snan(p: u32) -> BigFloat {
    BigFloat::try_new_signaling_nan(Sign::Positive, p, &[]).expect("p >= 1")
}

/// A complex value from two integer components at precision `p`.
pub fn cbf(re: i64, im: i64, p: u32) -> Complex<BigFloat> {
    Complex::new(bf(re, p), bf(im, p))
}

// ---------- The IEEE class model (the finite special-value grid) ----------

/// A coarse IEEE 754 class of a real component: the eight values the complex
/// special-value dispatch branches on. `PosFin` / `NegFin` stand for any finite
/// nonzero value of that sign; the enumeration uses a representative.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cls {
    NegInf,
    NegFin,
    NegZero,
    PosZero,
    PosFin,
    PosInf,
    QNan,
    SNan,
}

/// Every class, for the exhaustive dispatch grid.
pub const ALL_CLASSES: [Cls; 8] = [
    Cls::NegInf,
    Cls::NegFin,
    Cls::NegZero,
    Cls::PosZero,
    Cls::PosFin,
    Cls::PosInf,
    Cls::QNan,
    Cls::SNan,
];

/// A `BigFloat` representative of a class at precision `p` (the finite reps are
/// `±2`, chosen so squares/exponentials stay well-defined).
pub fn rep(c: Cls, p: u32) -> BigFloat {
    match c {
        Cls::NegInf => ninf(p),
        Cls::NegFin => bf(-2, p),
        Cls::NegZero => nz(p),
        Cls::PosZero => pz(p),
        Cls::PosFin => bf(2, p),
        Cls::PosInf => pinf(p),
        Cls::QNan => qnan(p),
        Cls::SNan => snan(p),
    }
}

/// The class of a finished `BigFloat` result component. Distinguishes the
/// signed zeros (which IEEE comparison treats as equal), so a wrong-signed zero
/// is caught.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ResultCls {
    NegInf,
    NegFin,
    NegZero,
    PosZero,
    PosFin,
    PosInf,
    Nan,
}

/// Classify a result component for the table assertions.
pub fn classify(v: &BigFloat) -> ResultCls {
    if v.is_nan() {
        ResultCls::Nan
    } else if v.is_infinite() {
        if v.is_sign_negative() {
            ResultCls::NegInf
        } else {
            ResultCls::PosInf
        }
    } else if v.is_zero() {
        if v.is_sign_negative() {
            ResultCls::NegZero
        } else {
            ResultCls::PosZero
        }
    } else if v.is_sign_negative() {
        ResultCls::NegFin
    } else {
        ResultCls::PosFin
    }
}

/// Assert a result component has the expected class, with a descriptive
/// message. `who` names the row (e.g. `"csqrt(-4 + 0i).im"`).
#[track_caller]
pub fn expect_cls(v: &BigFloat, want: ResultCls, who: &str) {
    let got = classify(v);
    assert_eq!(
        got, want,
        "{who}: expected {want:?}, got {got:?} (value {v})"
    );
}

/// Assert a finite result component equals the integer `n` exactly (and so is
/// the correct value, not merely the right class).
#[track_caller]
pub fn expect_int(v: &BigFloat, n: i64, who: &str) {
    let p = v.precision().max(64);
    assert_eq!(
        v.partial_cmp(&bf(n, p)).0,
        Some(Ordering::Equal),
        "{who}: expected exactly {n}, got {v}"
    );
}

/// `|v - target| < 2^(-(p - slack))`, the tolerance check the identity lanes
/// use (each op rounds, so the round-trips hold to a few ulps, not bit-exactly).
pub fn close(v: &BigFloat, target: &BigFloat, p: u32, slack: i64) -> bool {
    let d = v.sub(target, NE).0.abs();
    let tol = bf(1, p).scale_by_pow2(-(p as i64) + slack).0;
    matches!(
        d.partial_cmp(&tol).0,
        Some(Ordering::Less | Ordering::Equal)
    )
}

// ---------- Deterministic PRNG (no `rand` dependency) ----------

/// xorshift64; fixed seeds keep the lanes reproducible.
pub struct Rng(pub u64);
impl Rng {
    pub fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    /// A signed integer in `[-range, range]`.
    pub fn int(&mut self, range: i64) -> i64 {
        (self.next() % (2 * range as u64 + 1)) as i64 - range
    }
}

// ---------- Lefevre-Muller hard-to-round seeding ----------
//
// Hard-to-round scalars produce hard-to-enclose complex components: a
// component value sitting within sub-ULP of a rounding boundary stresses the
// directed-pair enclosure of `csqrt`/`cexp`/`clog` far harder than an integer
// does. The corpus is pfloat-core's committed CORE-MATH / Lefevre-Muller table,
// reused verbatim through `include!` so there is one source of truth and the
// MIT attribution header travels with it (the same path the ball lane uses).

mod lm_corpus {
    #![allow(dead_code)]
    // Inputs are CORE-MATH hard-to-round binary64 cases (MIT, attribution in
    // the included file's header); only the input field is used for seeding.
    include!("../../../tests/differential/lefevre_muller_data.rs");
}

/// The hard-to-round INPUT corpus for a scalar sub-kernel that a complex
/// elementary function composes through (`exp`/`ln`/`sin`/`cos`/`sqrt`). Only
/// the input bits are used (the seed); the complex result is checked
/// independently by the acb oracle, not against the binary64 output bits.
pub fn lm_inputs_for(scalar: &str) -> Option<&'static [(u64, u64)]> {
    use lm_corpus::*;
    Some(match scalar {
        "exp" => EXP_CASES,
        "ln" => LN_CASES,
        "sin" => SIN_CASES,
        "cos" => COS_CASES,
        _ => return None,
    })
}

/// Whether a binary64 bit pattern is finite and non-zero (a usable seed).
pub fn is_finite_nonzero_f64(bits: u64) -> bool {
    let exp_field = (bits >> 52) & 0x7ff;
    let mant = bits & 0x000f_ffff_ffff_ffff;
    exp_field != 0x7ff && !(exp_field == 0 && mant == 0)
}

/// A `BigFloat` at precision `p` (`>= 53`) equal to the binary64 value with bit
/// pattern `bits` EXACTLY: the integer significand times an exact power-of-two
/// scale (the bit-exact route; a decimal round-trip would round at the finer
/// `BigFloat` precision and lose the hard-to-round property). Caller filters
/// non-finite / zero via [`is_finite_nonzero_f64`].
pub fn bf_of_f64_bits(bits: u64, p: u32) -> BigFloat {
    debug_assert!(p >= 53, "binary64 significand needs p >= 53 to stay exact");
    let sign = bits >> 63;
    let exp_field = (bits >> 52) & 0x7ff;
    let mant = bits & 0x000f_ffff_ffff_ffff;
    let (int_mant, scale) = if exp_field == 0 {
        (mant, -1074i64) // subnormal
    } else {
        (mant | 0x0010_0000_0000_0000, exp_field as i64 - 1023 - 52)
    };
    let mut v = BigFloat::try_from_i64_exact(int_mant as i64, p)
        .expect("binary64 significand fits i64 and p >= 53");
    v = v.scale_by_pow2(scale).0;
    if sign == 1 {
        v = v.negated();
    }
    v
}
