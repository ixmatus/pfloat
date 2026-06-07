//! Cross-check `BigFloat::to_f32_round` / `to_f64_round` (the pure-Rust
//! `no_std` conversion in `src/convert.rs`) against the rug/MPFR
//! reference. This is the lane that puts pfloat's OWN `BigFloat`-to-grid
//! rounding under test (the conversion the `pfloat-libm` shell and
//! `pfloat-ball` depend on), as opposed to the kernel-value oracle sweep.
//!
//! Two shapes:
//!
//! - A broad random sweep over off-grid values and all five modes (f32
//!   via `round_f32`, f64 via `round_f64`, both synthesizing
//!   `NearestAway`). A real rounding decision is made on every input.
//! - A boundary-complete corpus (pf-3rtr.8): exact grid points, exact
//!   midpoints (the `NearestEven` / `NearestAway` tie inputs), values a
//!   hair either side of the midpoint, the subnormal floor, and the
//!   overflow threshold. These are the inputs where a directed rounding
//!   decision is most delicate, and where the directed-mode bug-hunt's
//!   named failure mode (a boundary landing a hair the wrong way) would
//!   surface in the conversion itself.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod oracle;

use oracle::convert::{bigfloat_to_rug, round_f32, round_f64};
use pfloat::{BigFloat, RoundingMode};

const MODES: [RoundingMode; 5] = [
    RoundingMode::NearestEven,
    RoundingMode::NearestAway,
    RoundingMode::TowardZero,
    RoundingMode::TowardPositive,
    RoundingMode::TowardNegative,
];

const ITERS: usize = 200_000;

/// Deterministic 64-bit LCG (Knuth MMIX constants). Pure-function
/// generator so the sweep is reproducible without a `rand` dependency.
fn next(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state
}

/// A high-precision off-grid value: the 128-bit quotient of two random
/// finite f64s, which lands off both the f32 and the f64 grid.
fn ratio(a_bits: u64, b_bits: u64) -> Option<BigFloat> {
    let a = f64::from_bits(a_bits);
    let b = f64::from_bits(b_bits);
    if !a.is_finite() || !b.is_finite() || b == 0.0 {
        return None;
    }
    let a128 = BigFloat::from_f64(a)
        .round_to_precision(128, RoundingMode::NearestEven)
        .unwrap()
        .0;
    let b128 = BigFloat::from_f64(b)
        .round_to_precision(128, RoundingMode::NearestEven)
        .unwrap()
        .0;
    Some(a128.div(&b128, RoundingMode::NearestEven).0)
}

#[test]
fn to_f32_round_matches_mpfr() {
    let mut state = 0x1234_5678_9abc_def0u64;
    for _ in 0..ITERS {
        // Two value families: p=53 from a random f64 (off the f32 grid),
        // and a p=128 ratio (off every grid).
        let direct = {
            let x = f64::from_bits(next(&mut state));
            if x.is_finite() {
                Some(BigFloat::from_f64(x))
            } else {
                None
            }
        };
        let ratio = ratio(next(&mut state), next(&mut state));

        for bf in [direct, ratio].into_iter().flatten() {
            for mode in MODES {
                let mine = bf.to_f32_round(mode).0.to_bits();
                let reference = round_f32(&bigfloat_to_rug(&bf), mode)
                    .expect("finite value rounds to Some")
                    .to_bits();
                assert_eq!(
                    mine, reference,
                    "to_f32_round {mode:?}: pfloat 0x{mine:08x} != mpfr 0x{reference:08x}"
                );
            }
        }
    }
}

#[test]
fn to_f64_round_matches_mpfr() {
    let mut state = 0x0fed_cba9_8765_4321u64;
    for _ in 0..ITERS {
        let Some(bf) = ratio(next(&mut state), next(&mut state)) else {
            continue;
        };
        let rug_val = bigfloat_to_rug(&bf);
        // All five modes directly, including a synthesized NearestAway
        // (pf-3rtr.8): the f64 lane no longer leans on the f32 lane for NA.
        for mode in MODES {
            let mine = bf.to_f64_round(mode).0.to_bits();
            let reference = round_f64(&rug_val, mode)
                .expect("finite value rounds to Some")
                .to_bits();
            assert_eq!(
                mine, reference,
                "to_f64_round {mode:?}: pfloat 0x{mine:016x} != mpfr 0x{reference:016x}"
            );
        }
    }
}

// --- Boundary-complete corpus -------------------------------------------

/// `x / 2^k`, exact (a chain of power-of-two divisions only shifts the
/// binary exponent).
fn scale_down(x: &BigFloat, k: u32) -> BigFloat {
    let two = BigFloat::try_from_i64_exact(2, x.precision()).unwrap();
    let mut v = x.clone();
    for _ in 0..k {
        v = v.div(&two, RoundingMode::NearestEven).0;
    }
    v
}

/// Assert pfloat's `to_f32_round` matches the oracle on `input` under
/// every mode.
fn check_f32(input: &BigFloat, label: &str) {
    let rug = bigfloat_to_rug(input);
    for mode in MODES {
        let mine = input.to_f32_round(mode).0.to_bits();
        let want = round_f32(&rug, mode).expect("finite -> Some").to_bits();
        assert_eq!(
            mine, want,
            "{label} to_f32_round {mode:?}: pfloat 0x{mine:08x} != mpfr 0x{want:08x}"
        );
    }
}

/// Assert pfloat's `to_f64_round` matches the oracle on `input` under
/// every mode.
fn check_f64(input: &BigFloat, label: &str) {
    let rug = bigfloat_to_rug(input);
    for mode in MODES {
        let mine = input.to_f64_round(mode).0.to_bits();
        let want = round_f64(&rug, mode).expect("finite -> Some").to_bits();
        assert_eq!(
            mine, want,
            "{label} to_f64_round {mode:?}: pfloat 0x{mine:016x} != mpfr 0x{want:016x}"
        );
    }
}

#[test]
fn convert_boundary_corpus_f32() {
    // Representable f32 bit patterns spanning subnormal, smallest normal,
    // midrange, and the largest finite values (each below the next so the
    // b+1 neighbour exists).
    const REPS: &[u32] = &[
        0x0000_0001, // smallest +subnormal
        0x0000_0002,
        0x0040_0000, // mid subnormal
        0x007F_FFFE,
        0x007F_FFFF, // largest +subnormal
        0x0080_0000, // smallest +normal
        0x0080_0001,
        0x3F80_0000, // 1.0
        0x4049_0FDA, // ~pi neighbour
        0x7F00_0000, // large
        0x7F7F_FFFE, // max - ulp
    ];
    for &b in REPS {
        let lo = BigFloat::from_f32(f32::from_bits(b))
            .round_to_precision(80, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let hi = BigFloat::from_f32(f32::from_bits(b + 1))
            .round_to_precision(80, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let two = BigFloat::try_from_i64_exact(2, 80).unwrap();
        let mid = lo
            .add(&hi, RoundingMode::NearestEven)
            .0
            .div(&two, RoundingMode::NearestEven)
            .0;
        let ulp = hi.sub(&lo, RoundingMode::NearestEven).0;
        let eps = scale_down(&ulp, 30);
        // Exact grid point, exact midpoint (the tie), and a hair either
        // side of the midpoint, for +x and the mirrored −x.
        for base in [
            lo.clone(),
            mid.clone(),
            mid.add(&eps, RoundingMode::NearestEven).0,
            mid.sub(&eps, RoundingMode::NearestEven).0,
        ] {
            check_f32(&base, &format!("f32 b=0x{b:08x}"));
            check_f32(&base.negated(), &format!("f32 -b=0x{b:08x}"));
        }
    }

    // Overflow boundary: max_finite and values straddling max + ulp/2,
    // where nearest rounds to ±∞ and the directed modes split.
    let max = BigFloat::from_f32(f32::MAX)
        .round_to_precision(80, RoundingMode::NearestEven)
        .unwrap()
        .0;
    let ulp_max = BigFloat::from_f32(2f32.powi(104)) // ulp at f32::MAX
        .round_to_precision(80, RoundingMode::NearestEven)
        .unwrap()
        .0;
    for frac_k in [2u32, 1, 0] {
        // max + ulp/4, max + ulp/2 (the tie), max + ulp (= overflow).
        let add = scale_down(&ulp_max, frac_k);
        let v = max.add(&add, RoundingMode::NearestEven).0;
        check_f32(&v, &format!("f32 overflow +ulp/2^{frac_k}"));
        check_f32(&v.negated(), &format!("f32 overflow -ulp/2^{frac_k}"));
    }
}

#[test]
fn convert_boundary_corpus_f64() {
    const REPS: &[u64] = &[
        0x0000_0000_0000_0001, // smallest +subnormal
        0x0000_0000_0000_0002,
        0x0008_0000_0000_0000, // mid subnormal
        0x000F_FFFF_FFFF_FFFF, // largest +subnormal
        0x0010_0000_0000_0000, // smallest +normal
        0x0010_0000_0000_0001,
        0x3FF0_0000_0000_0000, // 1.0
        0x4009_21FB_5444_2D18, // ~pi neighbour
        0x7FE0_0000_0000_0000, // large
        0x7FEF_FFFF_FFFF_FFFE, // max - ulp
    ];
    for &b in REPS {
        let lo = BigFloat::from_f64(f64::from_bits(b))
            .round_to_precision(120, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let hi = BigFloat::from_f64(f64::from_bits(b + 1))
            .round_to_precision(120, RoundingMode::NearestEven)
            .unwrap()
            .0;
        let two = BigFloat::try_from_i64_exact(2, 120).unwrap();
        let mid = lo
            .add(&hi, RoundingMode::NearestEven)
            .0
            .div(&two, RoundingMode::NearestEven)
            .0;
        let ulp = hi.sub(&lo, RoundingMode::NearestEven).0;
        let eps = scale_down(&ulp, 40);
        for base in [
            lo.clone(),
            mid.clone(),
            mid.add(&eps, RoundingMode::NearestEven).0,
            mid.sub(&eps, RoundingMode::NearestEven).0,
        ] {
            check_f64(&base, &format!("f64 b=0x{b:016x}"));
            check_f64(&base.negated(), &format!("f64 -b=0x{b:016x}"));
        }
    }

    // Overflow boundary at f64::MAX, ulp there is 2^(1023-52) = 2^971.
    let max = BigFloat::from_f64(f64::MAX)
        .round_to_precision(120, RoundingMode::NearestEven)
        .unwrap()
        .0;
    let ulp_max = BigFloat::from_f64(2f64.powi(971))
        .round_to_precision(120, RoundingMode::NearestEven)
        .unwrap()
        .0;
    for frac_k in [2u32, 1, 0] {
        let add = scale_down(&ulp_max, frac_k);
        let v = max.add(&add, RoundingMode::NearestEven).0;
        check_f64(&v, &format!("f64 overflow +ulp/2^{frac_k}"));
        check_f64(&v.negated(), &format!("f64 overflow -ulp/2^{frac_k}"));
    }
}
