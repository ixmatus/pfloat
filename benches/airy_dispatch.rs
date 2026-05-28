//! Sub-slice 2b.2.b: Airy asymptotic-vs-Maclaurin dispatch baseline bench.
//!
//! Measures `BigFloat::ai_round` at input magnitudes straddling
//! `airy_threshold_exponent` (`src/math/airy.rs:392-406`). The
//! threshold derives from the optimal-truncation accuracy law
//! `e^{−2√ζ}` with `ζ = (2/3)|x|^{3/2}` and grows as the cube of
//! `|x|³` relative to precision (`|x|³ ≥ (9/4)·((p+G)/(2·log₂e))⁴`).
//!
//! Unlike the Bessel `bessel_j_threshold` (which over-cut the strict
//! accuracy bound by ~2.88× in `|x|`, leaving big slack for the
//! sub-slice 2b.2.a per-kernel split), the Airy formula already
//! sits at ~1.5% over strict via the integer approximation
//! `23/8 ≈ 2·log₂e`. No threshold-tightening is mathematically
//! reachable (sub-slice 2b.2.b closes documentation-tier per
//! ADR-0040 GATE A precedent; see ADR-0044). This bench is the
//! durable measurement infrastructure for any future Airy work
//! (`airy_asymptotic_pos`/`neg` boost reduction, `airy_series`
//! guard tightening, boundary-constant memoisation).
//!
//! Cells: `precision ∈ {53, 256} × |x| ∈ {2^T, 2^(T+1)} × Ai` = 4
//! cells, **all on the asymptotic side of the dispatch**. `T(53) =
//! 7` gives the asymptotic-side `|x| ∈ {128, 256}`; `T(256) = 10`
//! gives `|x| ∈ {1024, 2048}`. The `|x|` values are integers parsed
//! exactly at the target precision.
//!
//! Below-threshold Maclaurin cells are **explicitly excluded**:
//! - `p=53, |x|=64` (Maclaurin, working ≈ 549 bits): doesn't
//!   exercise the asymptotic path that's the target of any future
//!   `airy_asymptotic_pos`/`neg` tightening.
//! - `p=256, |x|=512` (Maclaurin, working ≈ 11 000 bits):
//!   prohibitively slow per iteration (~2-3 minutes per call;
//!   the sub-slice 2b.2.b first baseline attempt was killed
//!   mid-cell under contention, and the post-ferrodec re-attempt
//!   for the ADR-0048 baseline took 51 minutes without finishing
//!   this single cell). With `sample_size = 10` criterion needs
//!   ~half-hour to collect samples; not worth the wait for any
//!   change targeting the asymptotic path.
//! - `p=1024, |x|=2048` Maclaurin would have working ≈ 90 000 bits
//!   — fundamentally infeasible.
//!
//! Future Airy work targeting the **Maclaurin path** (e.g.,
//! tightening the `(2/3)|x|^{3/2}·log₂e` cancellation guard or
//! adding the deferred boundary-constant memoisation in
//! `airy_zero_value`) should add a dedicated bench at smaller
//! `|x|` (say `|x| ≤ 32`, where the boost is bounded to a few
//! hundred bits), with `sample_size = 3` or similar to fit the
//! cell into a finite budget.
//!
//! Single kernel (`Ai`): the threshold change affects all four Airy
//! functions identically (one dispatch, four output formulas);
//! benching Ai alone captures the cost profile.
//!
//! Run: `cargo bench --bench airy_dispatch --features airy`.
//! `harness = false`; not part of `cargo test`. Save baseline with
//! `-- --save-baseline phase2b-airy-baseline`; re-run with
//! `--baseline phase2b-airy-baseline` after any future change.

use core::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use pfloat::{BigFloat, RoundingMode};

/// `(label, target_precision, |x|)` cells. All four cells dispatch
/// to `airy_asymptotic_pos` (or `_neg` for negative arguments).
/// One cell at the threshold (`|x| = 2^T`); one above (`|x| =
/// 2^(T+1)`). Per-cell cost is bounded: at `p = 256, |x| = 2048`
/// the asymptotic uses `N ≈ √ζ ≈ √(2/3·2048^{3/2}) ≈ 248` terms at
/// working ≈ 320 bits, ~100 ms per call.
const CELLS: &[(&str, u32, u32)] = &[
    // p=53: T(53) = 7, |x| ≥ 128 enters asymptotic.
    ("p53_x128", 53, 128),
    ("p53_x256", 53, 256),
    // p=256: T(256) = 10, |x| ≥ 1024 enters asymptotic.
    ("p256_x1024", 256, 1024),
    ("p256_x2048", 256, 2048),
];

/// All four cells are asymptotic; per-call cost ranges from ~ms at
/// `p=53` to ~hundreds of ms at `p=256, |x|=2048`. Total bench
/// runs in ~2-3 minutes on a quiet machine.
const MEASUREMENT_TIME: Duration = Duration::from_secs(20);
const WARMUP_TIME: Duration = Duration::from_secs(2);
const SAMPLE_SIZE: usize = 10;

/// `Ai(x)` at each cell. The threshold change would affect all four
/// Airy functions identically (one dispatch site, four output
/// formulas all routed through the same regime selector); a single
/// kernel suffices to characterise the dispatch boundary cost.
fn bench_airy(c: &mut Criterion) {
    let mut group = c.benchmark_group("airy_dispatch");
    group.measurement_time(MEASUREMENT_TIME);
    group.warm_up_time(WARMUP_TIME);
    group.sample_size(SAMPLE_SIZE);

    for &(label, target, x_int) in CELLS {
        // Integer `|x|` is exact at any precision ≥ ⌈log₂|x|⌉ + 1.
        let x_str = x_int.to_string();
        let (x, _) =
            BigFloat::parse_str(&x_str, target, RoundingMode::NearestEven).expect("integer parses");

        let id = BenchmarkId::from_parameter(format!("Ai_{label}"));
        group.bench_with_input(id, &(target, x), |bench, (t, x)| {
            bench.iter(|| {
                black_box(x)
                    .ai_round(*t, RoundingMode::NearestEven)
                    .expect("target_precision >= 1")
            });
        });
    }
    group.finish();
}

criterion_group!(benches, bench_airy);
criterion_main!(benches);
