# ADR-0047: Bessel Miller seed-index tightening — accepted, precision-gated

- **Status**: accepted
- **Date**: 2026-05-27

## Context

Phase 2b sub-slice 2b.3 (`pf-6fvx` work, branch `phase-2b-perf-5`,
the final sub-slice of Phase 2b before pf-6fvx closes) targeted the
Miller backward-recurrence starting index `M` in `bessel_j_miller`
and `bessel_i_miller`. The pre-2b.3 implementation used:

1. **Exponential search** from `start = max(m_floor, 1 << (e_x + 2))`,
   doubling `big_m` until `satisfies(big_m)` (the DLMF 10.19.1
   criterion `M·(1 + ln(x/(2M))) ≤ −P·ln2`, `P = target + 64`).
2. **Plus a fixed `+8` safety guard** for `lp = 64`-bit-precision
   evaluation noise in the `satisfies` test.

The exponential search alone overshoots by up to 2× (a doubled
boundary); the `+8` guard adds more. Combined, `M_actual` was up to
~2.4× the optimum. The pre-2b.3 doc-comment explicitly named this
slack:

> overshoot ≤ 2× is a deliberate robustness/cost trade; retune only
> with a bench, CLAUDE.md

This is the same "perf-tuning deferred" pattern that yielded the
Bessel threshold win in ADR-0043 (and contrasted with the
first-principles-tight Airy, Si/Ci, erf/erfc thresholds in
ADR-0044/0045/0046). It was therefore expected to harvest a real
win.

## What was tried

**Stage 1: Unconditional binary refine.** After the exponential
search brackets `[lo, hi]`, binary-search inside the bracket for the
smallest `M` satisfying. Reduce the safety guard from `+8` to `+4`
(smaller margin sufficient once the search converges tightly).
Apply to both `bessel_j_miller` and `bessel_i_miller` via a new
shared `pub(super) fn miller_seed_m` helper in
`src/math/bessel_j.rs`.

**Stage 1 bench result** vs the pre-2b.3 clean baseline:

| Cell | Baseline | Stage 1 | Δ | Outcome |
|------|---------:|--------:|---|---------|
| J0_p256_x256 (Miller, p<512) | 0.99 ms | 1.71 ms | **+74%** | REGRESSION |
| J0_p1024_x1024 (Miller, p≥512) | 17.74 ms | 10.52 ms | −41% | win |

The regression on the small-precision Miller cell exposed the
structural cost of the binary refine. Diagnosis:

- Binary search does ~`log₂((hi − lo))` `satisfies` evaluations,
  each at `lp = 64`-bit precision (an `ln + add + mul + div + …`
  chain). On `aarch64-apple-darwin` the fixed overhead is ~0.4 ms
  for typical brackets (~10 binary-search iterations × ~40 µs each).
- Miller iteration cost scales with the recurrence's working
  precision (`target` plus the `|x|·log₂e` cancellation boost) —
  roughly `O(target²)`. At `target = 1024`, ~15 µs per Miller
  iteration; saving 500 iterations recovers 7.5 ms, dwarfing the
  0.4 ms search overhead. At `target = 256`, ~1.5 µs per Miller
  iteration; saving 500 iterations recovers 0.75 ms, less than
  the search overhead.

**Break-even is around `target ≈ 512`.** Below that the binary
search costs more than it saves; above it the savings dominate.

**Stage 2: Precision-gate the binary refine.** Use the binary refine
only when `target_precision >= 512`; below that, retain the
pre-2b.3 `+8` guard with no binary refine.

```rust
let big_m = if target_precision >= 512 {
    // Binary refine inside [lo, hi], +4 safety margin.
    while lo + 1 < hi {
        let mid = lo + (hi - lo) / 2;
        if satisfies(mid) { hi = mid; } else { lo = mid; }
    }
    hi + 4
} else {
    // Pre-2b.3 behaviour: exp search + 8 guard, no binary refine.
    hi + 8
};
big_m.min(cap).max(m_floor)
```

## Decision

**Land the precision-gated Miller seed refine.** The helper lives
in `src/math/bessel_j.rs` as `pub(super) fn miller_seed_m(ax,
target_precision, e_x, m_floor) -> i64`; both `bessel_j_miller` and
`bessel_i_miller` call it.

## Measurement (clean machine)

Methodological note: the initial bench runs (`bti2op662`,
`brk2n5ygn`) were captured under heavy CPU contention (concurrent
`ferrodec` exhaustive verification at load average 14-17). The
contended numbers showed uniform 15-22% drift on asymptotic-only
cells which I attributed to contention noise in the draft.
**After the ferrodec sweep ended, I re-baselined and re-benched on
the quiet machine** (load average ~1.8). The clean numbers
empirically confirm the noise hypothesis: the asymptotic-only drift
disappears, the true silicon noise floor is ±5%, and exactly one
cell shows a real Miller-refine effect.

**Clean comparison** (`buva2ql3n` vs `phase2b-bessel-clean`):

| Cell | Clean baseline | Stage 2 | Δ | Notes |
|------|---------------:|--------:|---|-------|
| **J0_p1024_x1024** (Miller, p≥512) | 17.74 ms | **10.04 ms** | **−43.4%** | **the slice deliverable** |
| J0_p256_x256 (Miller, p<512, no binary refine) | 0.99 ms | 1.04 ms | +5.9% | within noise; gate skips refine |
| J0_p256_x512 | 9.35 ms | 8.73 ms | −6.6% | noise |
| J0_p256_x1024 | 16.99 ms | 16.51 ms | −2.8% | noise |
| J0_p1024_x2048 | 127.49 ms | 125.44 ms | −1.6% | noise |
| J0_p1024_x4096 | 247.33 ms | 250.64 ms | +1.3% | noise |
| Y0/I0/K0 cells (16 total) | various | various | −3.8% to +2.3% | all within noise |

**The headline:** at the canonical high-precision Miller cell
(`J0_p1024_x1024`, where `J` dispatches to Miller via the
conservative `bessel_j_threshold` and the working precision is
large enough that the binary refine pays off), the Miller seed
tightening saves **~43% per call** (17.74 → 10.04 ms). All other
cells stay within the ±5% silicon noise floor of a quiet
`aarch64-apple-darwin` M-series box.

**Verification:** library unit tests 687/687 pass; `differential_jn`
7/7, `differential_ik` 7/7; the in-module Wronskian and recurrence
tests confirm no behavioural change.

## Why I-family cells don't show wins in this bench

`bessel_yik_threshold` (ADR-0043) moves I to asymptotic at the
dispatch boundary used in the bench grid (`p ∈ {256, 1024}`,
`|x| ∈ {T-1, T, T+1}`). The bench cells for I therefore don't
exercise the Miller path. The Miller-seed refine still applies to
I-side calls with `|x|` below the yik threshold (smaller `|x|`
than this bench exercises), but no bench cell currently catches it.
Future I-side perf work could probe `p ≥ 512, |x| ∈ [tiny, T-2]`
cells where `bessel_i_miller` is dispatched and the binary refine
activates.

## Consequences

**Pf-6fvx closes ACCEPTED.** All five Phase 2b sub-areas done:

| sub-slice | outcome | ADR |
|-----------|---------|-----|
| 2b.1 (Spouge precision-pegging) | REJECTED | 0041 |
| 2b.2.a (Bessel threshold per-kernel split) | ACCEPTED | 0043 |
| 2b.2.b (Airy threshold) | doc-tier | 0044 |
| 2b.2.c (Si/Ci threshold) | doc-tier | 0045 |
| 2b.2.d (erf/erfc threshold) | doc-tier | 0046 |
| **2b.3 (Bessel Miller seed)** | **ACCEPTED** | **0047 (this)** |

Plus the pf-1axr root fix (ADR-0042) surfaced mid-stream by the
2b.2.a boundary-input additions. **Phase 2b is done.**

**General lesson (combined with ADR-0044/0045/0046):** Two patterns
recur in numerical kernel perf tuning:

1. **Threshold formulas derived from first-principles accuracy laws
   land at the mathematical optimum.** Airy, Si/Ci, erf/erfc all
   confirmed (ADR-0044/0045/0046). Only Bessel's
   `bessel_j_threshold` (slice 6o "not perf-tuned without a bench")
   had room — captured by ADR-0043's per-kernel split.

2. **Search-based parameter selection often carries `2×` or larger
   overshoot from the choice of step size, with safety guards
   compounding the overshoot.** Bessel's Miller seed `M` exposed
   this (the doc-comment explicitly named the slack). The fix:
   binary refine inside the exponential bracket, gated on the
   per-iteration cost actually exceeding the search overhead (here
   `target_precision >= 512`). Transfers to any similar search-based
   parameter selection where per-call iteration cost dominates the
   search cost.

Both patterns generalise beyond Bessel: any future investigator
tuning a kernel asymptotic threshold or a recurrence depth should
first determine whether the existing formula derives from a
first-principles bound (likely already-optimal, doc-tier closure)
or carries an explicit conservative guard (likely tightenable with
a bench, real-win candidate).

**Methodological note on bench noise floor.** The clean re-bench
captured here on `aarch64-apple-darwin` M-series shows the true
silicon noise floor is ±5% with `criterion`'s `measurement_time=20s
sample_size=20` settings. ADR-0043's earlier characterisation of a
"4-15% drift on bit-identical code paths attributed to natural
baseline-to-baseline variance" was actually CPU contention from a
concurrent ferrodec exhaustive run (load average 14-17); the clean
re-bench evaporates that drift. **For future Phase 2 / Phase 3 perf
work, capture baselines on a quiet machine** to avoid the
contention-noise inflation that ADR-0043 propagated. ADR-0043's
ACCEPT decision is robust (the 89-97% wins on flip cells dwarf any
contention noise by 5-100×); only the noise-floor language needs
the correction this ADR provides.

## Related

- ADR-0043 (Bessel per-kernel threshold split) — moved I-side
  Miller cells in this bench's grid to asymptotic, so the Miller
  refine's effect on I doesn't show up in the current bench;
  future investigation in smaller `|x|` regimes could expose it.
- ADR-0044 / 0045 / 0046 — the three doc-tier "threshold already at
  the optimum" precedents establishing the pattern the general
  lesson section contrasts with.
- `src/math/bessel_j.rs:317-407` — the new `miller_seed_m` helper
  with the precision-gated binary refine.
- `src/math/bessel_j.rs:418-419` — `bessel_j_miller` call site
  (reduced from ~40 lines to 2).
- `src/math/bessel_i.rs:477-479` — `bessel_i_miller` call site
  (reduced from ~30 lines to 2).
- `benches/bessel_dispatch.rs` — the bench used for the
  measurement; cell grid covers Miller and asymptotic regimes at
  `p ∈ {256, 1024}`.
- `tests/differential_jn.rs:36-49`, `tests/differential_ik.rs:78-90`
  — boundary inputs landed by ADR-0042 / ADR-0043; exercise the
  recurrence at the dispatch boundary; pass under this change.
- `pf-6fvx` (Phase 2b kernel-specific perf) — closes with this ADR
  landed.
