//! Slice 7d: Karatsuba/Toom-Cook threshold calibration bench.
//!
//! Measures `BigFloat::mul` across operand limb counts spanning the
//! schoolbook -> Karatsuba crossover region so `KARATSUBA_THRESHOLD`
//! (`src/ops/limbs.rs`) can be set from a measurement rather than the
//! MPFR-ballpark guess of 30. The asm read (ADR-0027) established that
//! the schoolbook inner loop is an irreducibly serial scalar
//! multiply-accumulate (not vectorizable) and that the in-tree
//! Karatsuba pays several heap allocations per recursion node, so the
//! real crossover for *this* implementation is an empirical question,
//! host- and arch-dependent.
//!
//! Run: `cargo bench --bench mul_thresholds`. `harness = false`; this
//! is not part of `cargo test`.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use pfloat::{BigFloat, Parts, RoundingMode};

/// Build a finite `BigFloat` whose normalized mantissa occupies
/// exactly `limbs` 64-bit limbs and is densely non-zero (the
/// schoolbook `if ai == 0 { continue }` fast path is never taken, so
/// the bench measures the genuine O(n*m) cost). Built once per size
/// outside the timing loop.
///
/// A decimal string of `20 * limbs + 16` non-zero digits carries well
/// over `64 * limbs` bits (decimal digit ~= 3.32 bits), so rounding to
/// precision `64 * limbs` fills every limb with the top bit set.
fn dense_bigfloat(limbs: usize, seed: u64) -> BigFloat {
    let precision = u32::try_from(64 * limbs).expect("limb count fits u32 precision");
    let n_digits = 20 * limbs + 16;
    let mut s = String::with_capacity(n_digits);
    for i in 0..n_digits {
        // A non-trivial, non-repeating-per-limb digit pattern; first
        // digit forced non-zero so the value has the full width.
        let d = ((i as u64).wrapping_mul(7).wrapping_add(seed * 3 + 1)) % 9 + 1;
        s.push(char::from(b'0' + d as u8));
    }
    let (v, _status) = BigFloat::parse_str(&s, precision, RoundingMode::NearestEven)
        .expect("dense decimal literal parses");
    // Self-validation: a size-N sweep that silently used 1-limb
    // operands would be worthless.
    match v.parts() {
        Parts::Normal { mantissa, .. } => assert_eq!(
            mantissa.len(),
            limbs,
            "operand must occupy exactly {limbs} limbs"
        ),
        _ => panic!("dense_bigfloat produced a non-normal value"),
    }
    v
}

/// Limb counts swept across the crossover region. Dense around the
/// current threshold (30), then geometric out to where Karatsuba is
/// unambiguously winning.
const LIMB_SIZES: &[usize] = &[
    8, 16, 20, 24, 26, 28, 30, 32, 34, 36, 40, 48, 56, 64, 80, 96, 128, 192, 256, 384, 512,
];

/// Equal-size operands: the dispatcher's
/// `a.len().min(b.len()) <= KARATSUBA_THRESHOLD` reduces to the size
/// itself, so this isolates the algorithmic crossover.
fn bench_equal(c: &mut Criterion) {
    let mut group = c.benchmark_group("mul_equal");
    for &n in LIMB_SIZES {
        let a = dense_bigfloat(n, 1);
        let b = dense_bigfloat(n, 2);
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |bench, _| {
            bench.iter(|| black_box(&a).mul(black_box(&b), RoundingMode::NearestEven));
        });
    }
    group.finish();
}

/// Skewed operands (n by n/4): exercises the `min()` in the dispatch
/// and the unbalanced Karatsuba split, which a real workload (a wide
/// accumulator times a narrower factor) hits.
fn bench_skewed(c: &mut Criterion) {
    let mut group = c.benchmark_group("mul_skewed");
    for &n in LIMB_SIZES {
        let small = (n / 4).max(1);
        let a = dense_bigfloat(n, 3);
        let b = dense_bigfloat(small, 4);
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |bench, _| {
            bench.iter(|| black_box(&a).mul(black_box(&b), RoundingMode::NearestEven));
        });
    }
    group.finish();
}

criterion_group!(benches, bench_equal, bench_skewed);
criterion_main!(benches);
