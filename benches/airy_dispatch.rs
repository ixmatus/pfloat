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
//! Cells: `precision ∈ {53, 256} × |x| ∈ {2^(T−1), 2^T, 2^(T+1)} ×
//! Ai` = 6 cells. `T(53) = 7` so the straddle is `|x| ∈ {64, 128,
//! 256}`; `T(256) = 10` so the straddle is `|x| ∈ {512, 1024,
//! 2048}`. The `|x|` values are integers parsed exactly at the
//! target precision.
//!
//! `p = 1024` is dropped: `T(1024) = 12` puts the straddle at
//! `|x| ≈ 4096` and the below-threshold Maclaurin cell would have
//! working precision `1024 + 64 + (2/3)·2048^{3/2}·log₂e` ≈ 90 000
//! bits (prohibitive per call). Above-threshold p=1024 cells (e.g.
//! `|x| = 8192`) could bench fast but only exercise the asymptotic
//! at extremes; the p=256 cells already cover the asymptotic across
//! representative `|x|` and precision.
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

/// `(label, target_precision, |x|)` cells. The `|x|` integers
/// straddle `airy_threshold_exponent(target_precision)`: one binary
/// exponent below the current cut, one at it, one above.
const CELLS: &[(&str, u32, u32)] = &[
    // p=53: T(53) = 7, |x| ≥ 128 by the current cut.
    ("p53_x64", 53, 64),
    ("p53_x128", 53, 128),
    ("p53_x256", 53, 256),
    // p=256: T(256) = 10, |x| ≥ 1024 by the current cut.
    ("p256_x512", 256, 512),
    ("p256_x1024", 256, 1024),
    ("p256_x2048", 256, 2048),
];

/// The below-threshold Maclaurin path at `p = 256, |x| = 512` has
/// working precision around 11 000 bits with order-`|x|^{3/2}`
/// cancellation; per-call cost is on the order of seconds. Bump
/// `measurement_time` and drop `sample_size` accordingly to keep
/// the total run finite (~3-6 minutes across the 6 cells).
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
