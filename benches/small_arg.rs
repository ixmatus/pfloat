//! ADR-0059: small-argument fast-path calibration bench.
//!
//! The pf-lm3 exhaustive f32 sweep found `atanh` on tiny/subnormal x
//! (|x| < 2^-95) costs ~0.6-1.5 ms/input, 20-40x the normal-magnitude
//! cost, because `atanh_kernel` drives the full
//! `(log1p(x) - log1p(-x))/2` identity through the Ziv loop even though
//! for tiny x the value is x plus a cubic tail far below ULP.
//!
//! This bench measures atanh and the six odd siblings (asinh, atan,
//! sin, tan, sinh, tanh) on tiny |x| against a normal-magnitude control
//! (0.5 = 2^-1) per precision. The tiny-vs-control ratio is both the
//! ADR-0059 evidence for atanh and the measure-then-decide input for the
//! siblings: a sibling lands a `round_with_infinitesimal` short-circuit
//! only if its ratio shows a real win (strict revert stop-loss, ADR-0041
//! precedent); a neutral ratio reverts and the candidate is recorded as
//! measured-rejected.
//!
//! The tiny-x fast-path activates at `x.exponent <= -(target + 2)`, so
//! not every (exp, precision) cell is in-band: at p=24 all four tiny
//! exponents activate; at p=53 only -60/-95/-149; at p=113 only -149.
//! The out-of-band cells (e.g. 2^-30 at p=113) are themselves a slow
//! baseline that the fast-path does NOT touch, which is the intended
//! contrast.
//!
//! Run: `cargo bench --bench small_arg --features trig`. `harness =
//! false`; this is not part of `cargo test`.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use pfloat::{BigFloat, RoundingMode};

/// Binary exponents bracketing the activation band: 2^-30 sits just
/// out-of-band for the larger precisions, 2^-95 is the reported atanh
/// hotspot, 2^-149 is the smallest-subnormal f32 exponent.
const TINY_EXPS: &[i64] = &[-30, -60, -95, -149];

/// Target precisions: f32 (24), f64 (53), and binary128 (113) mantissa
/// widths, so the in-band cell count shifts with the `-(target + 2)`
/// threshold across the sweep.
const PRECISIONS: &[u32] = &[24, 53, 113];

/// `2^exp` at precision `p`, built exactly. A power of two is a
/// single set bit, so repeated multiply/divide by two only shifts the
/// stored exponent and never rounds; the value is exact at any
/// precision. Built once per input outside the timing loop.
fn pow2(exp: i64, p: u32) -> BigFloat {
    let two = BigFloat::try_from_i64_exact(2, p).expect("precision >= 1");
    let mut v = BigFloat::try_from_i64_exact(1, p).expect("precision >= 1");
    for _ in 0..exp.unsigned_abs() {
        v = if exp < 0 {
            v.div(&two, RoundingMode::NearestEven).0
        } else {
            v.mul(&two, RoundingMode::NearestEven).0
        };
    }
    v
}

/// One criterion group per function: a `2^-1` control plus the tiny
/// exponents, swept across `PRECISIONS`. The control and tiny rows at
/// the same precision give the slowdown ratio directly. `$method` is
/// the `BigFloat` kernel under test.
macro_rules! small_arg_group {
    ($name:ident, $method:ident, $group:literal) => {
        fn $name(c: &mut Criterion) {
            let mut g = c.benchmark_group($group);
            for &p in PRECISIONS {
                let control = pow2(-1, p);
                g.bench_with_input(BenchmarkId::new("control_2^-1", p), &p, |b, _| {
                    b.iter(|| black_box(&control).$method(RoundingMode::NearestEven));
                });
                for &e in TINY_EXPS {
                    let x = pow2(e, p);
                    g.bench_with_input(BenchmarkId::new(format!("2^{e}"), p), &p, |b, _| {
                        b.iter(|| black_box(&x).$method(RoundingMode::NearestEven));
                    });
                }
            }
            g.finish();
        }
    };
}

small_arg_group!(bench_atanh, atanh, "atanh_small");
small_arg_group!(bench_asinh, asinh, "asinh_small");
small_arg_group!(bench_atan, atan, "atan_small");
small_arg_group!(bench_sin, sin, "sin_small");
small_arg_group!(bench_tan, tan, "tan_small");
small_arg_group!(bench_sinh, sinh, "sinh_small");
small_arg_group!(bench_tanh, tanh, "tanh_small");

criterion_group!(
    benches,
    bench_atanh,
    bench_asinh,
    bench_atan,
    bench_sin,
    bench_tan,
    bench_sinh,
    bench_tanh
);
criterion_main!(benches);
