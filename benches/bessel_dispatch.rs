//! Sub-slice 2b.2.a: Bessel asymptotic-vs-Miller dispatch baseline bench.
//!
//! Measures the Bessel kernel quartet (`J₀`, `Y₀`, `I₀`, `K₀`) at
//! input magnitudes straddling the current `bessel_j_threshold`
//! (`src/math/bessel_j.rs:455-464`). All four kernels share the
//! threshold: `Y` / `I` / `K` import it via
//! `use super::bessel_j::bessel_j_threshold` at `bessel_y.rs:255`,
//! `bessel_i.rs:233`, `bessel_k.rs:250`. Tightening the formula
//! changes the dispatch boundary for all four simultaneously.
//!
//! The current `bessel_j_threshold(target_precision)` returns the
//! smallest binary exponent `e` such that `2^e ≥ target_precision +
//! 64`. The strict accuracy law `|x| ≳ (target+64)·ln2/2 ≈
//! 0.347·(target+64)` (`bessel_j.rs:446-448`) is ~2.88× smaller in
//! `|x|`; the conservative cut is the gap sub-slice 2b.2.a measures.
//! Per `bessel_j.rs:450-452`: "the crossover is not perf-tuned
//! without a bench, CLAUDE.md".
//!
//! Cells: `precision ∈ {256, 1024} × |x| ∈ {2^(T−1), 2^T, 2^(T+1)} ×
//! kernel ∈ {J₀, Y₀, I₀, K₀}` = 24 cells. `T(256) = 9` so the
//! straddle is `|x| ∈ {256, 512, 1024}`; `T(1024) = 11` so the
//! straddle is `|x| ∈ {1024, 2048, 4096}`. The `|x|` values are
//! integers parsed exactly at the target precision, so the timing
//! loop sees no upstream conversion work.
//!
//! `p = 4096` is dropped: `T(4096) ≈ 13` puts the straddle at `|x| ≈
//! 8192` and the asymptotic optimal-truncation index `N ≈ √(2|x|) ≈
//! 128` terms × 4 kernels makes the cell prohibitive (>1 second per
//! call at `working_prec` ~4160).
//!
//! Run: `cargo bench --bench bessel_dispatch --features bessel`.
//! `harness = false`; not part of `cargo test`. Save baseline with
//! `-- --save-baseline phase2b-bessel-baseline`; re-run with
//! `--baseline phase2b-bessel-baseline` after tightening to diff.

use core::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use pfloat::{BigFloat, RoundingMode};

/// `(label, target_precision, |x|)` cells. The `|x|` integers
/// straddle `bessel_j_threshold(target_precision)`: one binary
/// exponent below the current cut, one at it, one above.
const CELLS: &[(&str, u32, u32)] = &[
    // p=256: T(256) = 9, |x| ≥ 512 by the current cut.
    ("p256_x256", 256, 256),
    ("p256_x512", 256, 512),
    ("p256_x1024", 256, 1024),
    // p=1024: T(1024) = 11, |x| ≥ 2048 by the current cut.
    ("p1024_x1024", 1024, 1024),
    ("p1024_x2048", 1024, 2048),
    ("p1024_x4096", 1024, 4096),
];

/// Per-call cost at the worst cell (`p=1024, |x|=4096`, asymptotic
/// with `N ≈ 90` terms at `working_prec` ~1088) is order ~100 ms.
/// Below-threshold Miller-regime cells (`p=1024, |x|=1024`) carry a
/// seed index `M` that scales with `target_precision`; per-call cost
/// is similar order. Default `measurement_time = 5 s` would yield
/// very few samples at those cells; bump to 20 s and drop
/// `sample_size` to 20 to keep the run finite. 24 cells × 20 s ≈ 8
/// minutes plus per-cell warmup.
const MEASUREMENT_TIME: Duration = Duration::from_secs(20);
const WARMUP_TIME: Duration = Duration::from_secs(2);
const SAMPLE_SIZE: usize = 20;

/// `Jₙ` / `Yₙ` / `Iₙ` / `Kₙ` at order `n = 0` for each cell. Order
/// 0 picks the simplest kernel of each family; the threshold change
/// affects every order identically (the threshold depends only on
/// `target_precision`, not on `n`).
fn bench_bessel(c: &mut Criterion) {
    let mut group = c.benchmark_group("bessel_dispatch");
    group.measurement_time(MEASUREMENT_TIME);
    group.warm_up_time(WARMUP_TIME);
    group.sample_size(SAMPLE_SIZE);

    for &(label, target, x_int) in CELLS {
        // Integer `|x|` is exact at any precision ≥ ⌈log₂|x|⌉ + 1,
        // so the parse below is bit-exact at the target precision.
        let x_str = x_int.to_string();
        let (x, _) =
            BigFloat::parse_str(&x_str, target, RoundingMode::NearestEven).expect("integer parses");

        // J₀
        let id_j = BenchmarkId::from_parameter(format!("J0_{label}"));
        group.bench_with_input(id_j, &(target, x.clone()), |bench, (t, x)| {
            bench.iter(|| {
                black_box(x)
                    .j0_round(*t, RoundingMode::NearestEven)
                    .expect("target_precision >= 1")
            });
        });

        // Y₀
        let id_y = BenchmarkId::from_parameter(format!("Y0_{label}"));
        group.bench_with_input(id_y, &(target, x.clone()), |bench, (t, x)| {
            bench.iter(|| {
                black_box(x)
                    .yn_round(0, *t, RoundingMode::NearestEven)
                    .expect("target_precision >= 1")
            });
        });

        // I₀
        let id_i = BenchmarkId::from_parameter(format!("I0_{label}"));
        group.bench_with_input(id_i, &(target, x.clone()), |bench, (t, x)| {
            bench.iter(|| {
                black_box(x)
                    .in_round(0, *t, RoundingMode::NearestEven)
                    .expect("target_precision >= 1")
            });
        });

        // K₀
        let id_k = BenchmarkId::from_parameter(format!("K0_{label}"));
        group.bench_with_input(id_k, &(target, x), |bench, (t, x)| {
            bench.iter(|| {
                black_box(x)
                    .kn_round(0, *t, RoundingMode::NearestEven)
                    .expect("target_precision >= 1")
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_bessel);
criterion_main!(benches);
