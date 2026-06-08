//! Shared generators for the pfloat-ball verification lanes.
//!
//! Lifted verbatim from `property_ftia.rs` (pf-fe5f.1) so the
//! self-consistency lane and the independent Arb containment lane
//! (`differential_arb.rs`) draw on one ball generator and one exact
//! witness reconstruction. A tightening of the generator then improves
//! every lane at once, and a divergence between two lanes localizes to
//! the oracle rather than to the inputs.

// Each lane that `mod common;`s this file uses a different subset; the
// unused subset would otherwise warn under that binary's compilation.
#![allow(dead_code)]

use core::cmp::Ordering;
use pfloat::{BigFloat, Parts, RoundingMode};
use pfloat_ball::{Ball, Mag};

/// The Arb `BRACKET` subprocess driver + exact dyadic codec, gated to the
/// independent containment lane (it drives a python subprocess and needs
/// `std`). property_ftia does not pull `differential-arb`, so it never
/// compiles this submodule.
#[cfg(feature = "differential-arb")]
pub mod arb_bracket;

pub const NE: RoundingMode = RoundingMode::NearestEven;

/// Deterministic xorshift64 PRNG (no `rand` dependency; fixed seeds keep
/// the lanes reproducible).
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

pub fn bf(n: i64, p: u32) -> BigFloat {
    BigFloat::try_from_i64_exact(n, p).unwrap()
}

/// `lower <= x <= upper`.
pub fn contains(b: &Ball<BigFloat>, x: &BigFloat) -> bool {
    b.lower().partial_cmp(x).0 != Some(Ordering::Greater)
        && b.upper().partial_cmp(x).0 != Some(Ordering::Less)
}

/// A random ball at precision `p`: integer midpoint in `[-range, range]`
/// scaled by `2^scale`, radius `2^radexp` (or exact / entire).
pub fn random_ball(rng: &mut Rng, p: u32) -> Ball<BigFloat> {
    let m = rng.int(1 << 20);
    let scale = rng.int(40);
    let (mid, _) = bf(m, p).scale_by_pow2(scale);
    let rad = match rng.next() % 8 {
        0 => Mag::ZERO,
        _ => Mag::from_pow2(scale + rng.int(8) - 30),
    };
    Ball::new(mid, rad).unwrap()
}

/// Witnesses inside `[mid - rad, mid + rad]`, reconstructed EXACTLY (not
/// the outward-rounded `lower()`/`upper()`): `mid` and `mid ± rad·t` for
/// dyadic `t ∈ {0, ±1/2, ±1}` at high precision.
pub fn witnesses(b: &Ball<BigFloat>, work: u32) -> Vec<BigFloat> {
    let mid = b.midpoint().round_to_precision(work, NE).unwrap().0;
    let mut out = vec![mid.clone()];
    if let Mag::Finite { .. } = b.radius() {
        let rad = b
            .radius()
            .to_bigfloat()
            .round_to_precision(work, NE)
            .unwrap()
            .0;
        for &(num, den_pow) in &[(1i64, 0u32), (1, 1)] {
            let (scaled, _) = rad.scale_by_pow2(-(den_pow as i64));
            let (scaled, _) = scaled.mul(&bf(num, work), NE);
            out.push(mid.add(&scaled, NE).0);
            out.push(mid.sub(&scaled, NE).0);
        }
    }
    out
}

/// Binary exponent `floor(log2|v|)` of a finite non-zero `BigFloat`;
/// `i64::MIN` for zero / non-finite.
pub fn bin_exponent(v: &BigFloat) -> i64 {
    match v.parts() {
        Parts::Normal { exponent, .. } => exponent,
        _ => i64::MIN,
    }
}

// ---------- Lefevre-Muller hardest-to-round seeding (pf-vcqh) ----------
//
// Hard-to-round scalars produce hard-to-enclose balls: a midpoint where
// `f(mid)` sits within sub-ULP of a rounding boundary stresses the directed
// endpoints that define the ball radius far harder than a random integer
// midpoint does. The corpus is pfloat-core's committed CORE-MATH / Lefevre-
// Muller table (`tests/differential/lefevre_muller_data.rs`), reused verbatim
// through `include!` so there is one source of truth and the MIT attribution
// header travels with it. ADR-0078's deferred seeding item.

mod lm_corpus {
    #![allow(dead_code)]
    // Inputs are CORE-MATH hard-to-round binary64 cases (MIT, attribution in
    // the included file's header); only the input field is used for seeding.
    include!("../../../tests/differential/lefevre_muller_data.rs");
}

/// The hard-to-round corpus `(input_bits, expected_bits)` for a ball function
/// id, or `None` when the function has no corpus (`cbrt` / `sqrt`, and the
/// inverse-trig edge functions, keep random inputs only). The corpus is
/// per-function: an input hard to round for `f` makes `f(mid)` boundary-close,
/// which is what stresses `f`'s ball endpoints.
pub fn lm_cases_for(fn_id: &str) -> Option<&'static [(u64, u64)]> {
    use lm_corpus::*;
    Some(match fn_id {
        "exp" => EXP_CASES,
        "expm1" => EXPM1_CASES,
        "exp2" => EXP2_CASES,
        "exp10" => EXP10_CASES,
        "sinh" => SINH_CASES,
        "cosh" => COSH_CASES,
        "tanh" => TANH_CASES,
        "atan" => ATAN_CASES,
        "asinh" => ASINH_CASES,
        "sin" => SIN_CASES,
        "cos" => COS_CASES,
        "tan" => TAN_CASES,
        "ln" => LN_CASES,
        "log2" => LOG2_CASES,
        "log10" => LOG10_CASES,
        "log1p" => LOG1P_CASES,
        _ => return None,
    })
}

/// Whether a binary64 bit pattern is finite and non-zero (a usable seed
/// midpoint). The corpus is already finite by construction, but a ±0 or
/// special slipping through is filtered rather than panicking the builder.
pub fn is_finite_nonzero_f64(bits: u64) -> bool {
    let exp_field = (bits >> 52) & 0x7ff;
    let mant = bits & 0x000f_ffff_ffff_ffff;
    exp_field != 0x7ff && !(exp_field == 0 && mant == 0)
}

/// Build a `BigFloat` at precision `p` (which must be `>= 53`) that equals the
/// binary64 value with bit pattern `bits` EXACTLY: the integer significand
/// times an exact power-of-two scale (the bit-exact route; a shortest-decimal
/// round-trip would round at the finer `BigFloat` precision and lose the
/// hard-to-round property). Caller filters non-finite / zero via
/// [`is_finite_nonzero_f64`].
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

/// A ball whose midpoint is the hard-to-round binary64 value `bits`, at
/// precision 53 or 113 (where the 53-bit value is exact), with a small radius
/// bound to the midpoint magnitude (`2^(e-38 ..= e-22)`, the same relative
/// regime as [`random_ball`]) or exact 1/8 of the time. `bits` must be a
/// finite non-zero pattern ([`is_finite_nonzero_f64`]).
pub fn seeded_ball(rng: &mut Rng, bits: u64) -> Ball<BigFloat> {
    let p = if rng.next().is_multiple_of(2) {
        53
    } else {
        113
    };
    let mid = bf_of_f64_bits(bits, p);
    let e = bin_exponent(&mid);
    let rad = match rng.next() % 8 {
        0 => Mag::ZERO,
        _ => Mag::from_pow2(e + rng.int(8) - 30),
    };
    Ball::new(mid, rad).unwrap()
}
