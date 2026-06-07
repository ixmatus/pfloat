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
use pfloat::{BigFloat, RoundingMode};
use pfloat_ball::{Ball, Mag};

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
