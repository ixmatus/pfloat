//! Sub-slice 2b.1: Spouge precision-pegging baseline bench.
//!
//! Measures `BigFloat::lgamma_round` at working precisions where the
//! Spouge dispatch fires (`working_prec > STIRLING_REACH_THRESHOLD =
//! 600`, see `src/math/lgamma.rs:260`). The target precisions cover
//! the realistic Spouge-region: `1024` is just past the dispatch
//! threshold; `2048` is mid-range; `4096` is the asymptotic tail
//! where the current `spouge_a_for` formula's margin is largest in
//! absolute coefficient count.
//!
//! The current `spouge_a_for(working_prec) = working_prec/5 + 20`
//! formula (`src/math/gamma_stirling.rs:422`) is asymptotically
//! 1.5-1.6x over the strict Spouge truncation minimum. Sub-slice
//! 2b.1 tightens the margin against `LGAMMA_ERROR_GUARD = 24` plus a
//! small cancellation budget; the bench captures the speedup.
//!
//! Inputs vary `z` to capture different regimes:
//! - `z = 2.5` — small, near the pole at `z=0`, Spouge's primary use
//!   case (the partial-sum convergence is slowest here).
//! - `z = 10` — moderate, well past the pole, typical caller value.
//! - `z = 100` — large, where Stirling would normally take over but
//!   the dispatch still routes to Spouge above the threshold.
//!
//! Run: `cargo bench --bench spouge_lgamma`. `harness = false`; not
//! part of `cargo test`. Per-call cost at `working_prec = 1024` is
//! dominated by the `a ≈ 224` partial-sum-loop divisions; the bench
//! uses a longer `measurement_time` than the criterion default to
//! collect enough samples at the top precisions.

use core::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use pfloat::{BigFloat, RoundingMode};

/// Target precisions to bench. All push `working_prec` past the
/// `STIRLING_REACH_THRESHOLD = 600`, forcing the Spouge dispatch.
const TARGET_PRECISIONS: &[u32] = &[1024, 2048, 4096];

/// Input values of `z` for `ln Γ(z)`. The decimal literal at the
/// matching target precision parses to an exact representation at
/// the listed precision when `z` is rational; non-exact inputs (like
/// `2.5`) exercise the partial-sum convergence with cancellation.
const INPUTS: &[(&str, &str)] = &[
    ("z_2.5", "2.5"),
    ("z_10", "10"),
    ("z_100", "100"),
];

/// Criterion settings. Per-call cost ranges from ~50 µs at p=1024 to
/// ~50 ms at p=4096 (the latter requires ~840 div+add ops at 4096
/// bits each). Default `measurement_time = 5 s` would yield very
/// few samples at p=4096; bump to 20 s and drop `sample_size` to 20
/// so the run finishes in ~20 minutes total.
const MEASUREMENT_TIME: Duration = Duration::from_secs(20);
const WARMUP_TIME: Duration = Duration::from_secs(2);
const SAMPLE_SIZE: usize = 20;

/// `lgamma_round` at each `(target_precision, z)` cell.
fn bench_lgamma(c: &mut Criterion) {
    let mut group = c.benchmark_group("spouge_lgamma");
    group.measurement_time(MEASUREMENT_TIME);
    group.warm_up_time(WARMUP_TIME);
    group.sample_size(SAMPLE_SIZE);

    for &target in TARGET_PRECISIONS {
        for &(label, decimal) in INPUTS {
            // Parse the input at the target precision so the lgamma
            // call does no upstream conversion work inside the
            // timing loop.
            let (z, _) = BigFloat::parse_str(decimal, target, RoundingMode::NearestEven)
                .expect("decimal literal parses");
            let id = BenchmarkId::from_parameter(format!("p{}_{}", target, label));
            group.bench_with_input(id, &(target, z), |bench, (t, z)| {
                bench.iter(|| {
                    black_box(z)
                        .lgamma_round(*t, RoundingMode::NearestEven)
                        .expect("target_precision >= 1")
                });
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_lgamma);
criterion_main!(benches);
