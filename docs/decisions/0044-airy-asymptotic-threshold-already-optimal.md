# ADR-0044: Airy asymptotic threshold — already at the mathematical optimum

- **Status**: accepted (documentation-tier)
- **Date**: 2026-05-27

## Context

Phase 2b sub-slice 2b.2.b (`pf-6fvx` work, branch `phase-2b-perf-3`)
targeted the Airy asymptotic dispatch threshold
`airy_threshold_exponent` (`src/math/airy.rs:392-406`) for tightening,
following the per-kernel threshold-split methodology that landed ADR-0043
as the first Phase 2b ACCEPTED ADR. The hypothesis: an analogous over-cut
exists in Airy's formula and can be shrunk for similar wins.

The hypothesis was wrong.

## Analysis

Airy's dispatch boundary derives from the optimal-truncation accuracy law
`e^{−2√ζ}` (DLMF §9.7) where `ζ = (2/3)|x|^{3/2}`. Reaching `target + G` bits
of accuracy requires:

```text
2√ζ · log₂e ≥ target + G
⇒ ζ ≥ ((target + G) / (2·log₂e))²
⇒ (2/3)|x|^{3/2} ≥ ((target + G) / (2·log₂e))²
⇒ |x|³ ≥ (9/4) · ((target + G) / (2·log₂e))⁴
```

The code at `airy_threshold_exponent` requires
`8^{e_x} ≥ ⌈9·(p+32)⁴·8⁴ / (4·23⁴)⌉` using the integer rational
`23/8 = 2.875` for `2·log₂e ≈ 2.885`. Computing the constant:

| Bound | Coefficient on `(p+32)⁴` for `|x|³ ≥ K·(p+32)⁴` |
|-------|-----|
| Strict (true `log₂e ≈ 1.4427`) | `K_strict = 9/(4·(2·log₂e)⁴) ≈ 1/30.83` |
| Code (integer `23/8`) | `K_code = 9·4096/(4·279841) = 9216/1119364 ≈ 1/121.5` |

Wait — that ratio doesn't read right. Re-examining: the `|x|³` lower bound
the code requires is `K_code·(p+32)⁴`, the strict bound is
`K_strict·(p+32)⁴`. If `K_code > K_strict`, the code requires a larger `|x|`
than strict (more conservative). Numerically:

- `K_code = 9216/1119364 ≈ 8.234 × 10⁻³`
- `K_strict = 9/(4·(2.885)⁴) ≈ 9/277.6 ≈ 3.243 × 10⁻²`

Wait, that gives `K_code < K_strict`, which would mean the code requires
a **smaller** `|x|` than strict — non-conservative. That can't be right.
Let me re-derive.

Re-deriving once more without rounding: the strict requirement is
`|x|³ ≥ (9/4)·((p+G)/(2·log₂e))⁴`. Expanding the RHS:

```text
(9/4)·(p+G)⁴ / (2·log₂e)⁴
```

For `log₂e = 1.4427`, `2·log₂e = 2.8854`, `(2.8854)⁴ ≈ 69.30`. So
`K_strict = (9/4) / 69.30 = 2.25 / 69.30 ≈ 0.03247`.

For the integer `23/8 = 2.875`, `(2.875)⁴ ≈ 68.34`. So
`K_code = (9/4) / 68.34 ≈ 0.03293`.

So `K_code / K_strict ≈ 0.03293 / 0.03247 ≈ 1.014`. The code requires
`|x|³` to be ~1.4% larger than strict, i.e., `|x|` ~0.46% larger
(cube root). **The formula is at ~1.5% over the strict accuracy bound** —
essentially at the mathematical optimum.

Compare to Bessel's `bessel_j_threshold` (`2^{e_x} ≥ target+64`), which
over-cut the strict accuracy bound `|x| ≳ 0.347·(target+64)` by **2.88×**
in `|x|`. That gap was the room for the 2b.2.a per-kernel split
(ADR-0043). Airy has no analogous gap.

## What about the `+32` guard?

The `+32` in `(p+32)⁴` budgets above the target accuracy:
- `AIRY_ERROR_GUARD = 24` (Ziv-side calibrated slack per
  `src/math/ziv_calibration.rs:191`)
- Plus 8 bits margin for `log₂(N)` round-off accumulation in the
  asymptotic sum

Shrinking `G` from 32 to 24 (just the Ziv guard) is the analogous
tightening to ADR-0043's `+64 → ⌈(.../2)⌉` for Bessel. Computing the
threshold exponent at the four standard precisions:

| precision | T(G=32) | T(G=24) | Change |
|----------:|--------:|--------:|--------|
| 53        | 7       | 7       | none   |
| 113       | 8       | 8       | none   |
| 256       | 10      | 10      | none   |
| 1024      | 12      | 12      | none   |

The cube-of-`|x|³` growth makes the threshold sticky — even a meaningful
reduction in `G` doesn't move `T(p)` at any of the precisions on the
verification grid. **No threshold-exponent change is reachable from `G`
tightening alone.**

## Decision

**Close sub-slice 2b.2.b documentation-tier per the ADR-0040 GATE A
precedent** (the FFT measurement, which closed because Karatsuba covers
all consumer reach with decimal-orders headroom). No code change to the
threshold function ships.

**The bench infrastructure (`benches/airy_dispatch.rs` + the
`[[bench]]` entry in `Cargo.toml`) lands as durable measurement
infrastructure** for any future Airy perf work — for example:

- `airy_asymptotic_pos` / `airy_asymptotic_neg` use `working = target +
  64`; the asymptotic's round-off is `≈ log₂(N) ≤ 8` bits for any
  practical `|x|`, so a reduction to `working = target + 32` (just
  `AIRY_ERROR_GUARD + 8`) might save ~30% on every asymptotic call.
  Not pursued here; would be its own sub-slice.
- The `airy_series` boost `(2/3)|x|^{3/2}·log₂e` (`src/math/airy.rs:741-752`)
  could be re-examined for a tighter bound, though the Maclaurin path
  is only ever entered below the threshold where the boost is bounded.
- The `airy_zero_value` boundary constants (`src/math/airy.rs:326-370`)
  are recomputed on every `x=0` call; memoisation is explicitly deferred
  pending a bench. Could be the most impactful single change but is
  orthogonal to the threshold question.

**No baseline timings captured.** The bench was attempted but the
session's hardware was contended by a concurrent `ferrodec` exhaustive
verification run (load average 14-17, ferrodec at ~96% on most cores).
Under contention the worst Maclaurin cell at `p=256, |x|=512` (working
precision ~11 000 bits) stretched into tens of minutes per iteration
with noisy timings; the run was killed rather than yielding contention-
corrupted numbers. A clean baseline can be captured later when the
machine is quiet — the bench file is in tree and ready.

## Consequences

**No behavioural change to `Ai`/`Bi`/`Ai'`/`Bi'`.** The threshold and
all four kernels run unchanged. Verification (the Airy in-module tests,
`differential_ai`, `differential_bi`, `property_ai`) is unaffected.

**Phase 2b's ratio of ACCEPTED-to-doc-tier-to-rejected ADRs continues
to refine the measurement-as-deliverable discipline:**

| ADR | Outcome | Slice |
|----:|---------|-------|
| 0037 | REJECTED | SmallVec mantissa swap |
| 0040 | doc-tier (GATE A) | FFT measurement |
| 0041 | REJECTED | Spouge precision-pegging |
| 0042 | ACCEPTED | pf-1axr trig + bessel_y recurrence root fix |
| 0043 | ACCEPTED | Bessel per-kernel threshold split |
| **0044** | **doc-tier** | **Airy threshold — already optimal** |

The pattern: bench-first as truth (sometimes the bench confirms there's
nothing to do, sometimes it confirms a clean win, sometimes it surfaces
an asymmetric outcome requiring re-scoping). Each ADR records a
non-trivial finding the next investigator inherits.

**A note on the noise-floor characterisation in ADR-0043:** the 4-15%
drift on cells with bit-identical code paths was attributed to "natural
baseline-to-baseline variance on this hardware (aarch64-apple-darwin)".
More accurately: those measurements were taken under concurrent
`ferrodec` load (started ~3:14 PM the same day, ran for the duration
of the post-pf-1axr baseline and both comparison runs). The drift was
CPU-contention noise, not silicon variance; the actual quiet-machine
noise floor is likely smaller. ADR-0043's accept decision is robust —
the 89-97% wins on flip cells dwarf any contention noise by 5-100× —
but the noise-floor language could be tighter. Not amended; flagging
here for the next investigator reading the two ADRs together.

**Sub-slices remaining on `pf-6fvx`:**

- 2b.2.c Si/Ci asymptotic-series cutoffs — their existing tables already
  cover the asymptotic regime at TRANSCENDENTAL_PRECISIONS.
- 2b.2.d erf/erfc asymptotic-series cutoffs — needs added boundary
  inputs analogous to Bessel's `(1025, 1)` / `(2049, 1)` cases.
- 2b.3 Bessel Miller-recurrence depth — smallest scope; closes pf-6fvx.

The Si/Ci formula `2^{e_x} ≥ ⌈(p+32)·6932/10000⌉` (`src/math/si.rs:164-174`)
and erf/erfc formula `4^{e_x} ≥ ⌈(p+32)·16/23⌉`
(`src/math/erf.rs:181-192`) should be analytically vetted for the
same "already-optimal" structure before any bench is built. The Airy
finding suggests that asymptotic-threshold formulas derived from
first-principles accuracy laws (rather than the chunky conservative cut
that Bessel inherited from slice 6o's "not perf-tuned without a bench"
deferral) tend to be at or near the mathematical optimum. If 2b.2.c and
2b.2.d are similarly already-optimal, the remaining Phase 2b sub-slices
that ship code changes will be limited to 2b.3 (Miller depth) and any
follow-up Airy asymptotic-boost tightening per the §"Bench infrastructure"
list above.

## Related

- ADR-0023 / 0024 / 0025 (Bessel J / Y / I+K design) — the slices that
  landed the conservative `bessel_j_threshold` formula. The Airy
  threshold landed via ADR-0021 (the slice 6n Airy kernel) and was
  derived from first principles directly, hence its tighter starting
  point.
- ADR-0040 (FFT GATE A) — the doc-tier closure precedent for "bench
  measurement confirms no perf win available".
- ADR-0042 (pf-1axr) — the latent kernel defect surfaced by 2b.2.a's
  boundary-input additions; mentioned here because the Airy investigation
  does NOT add new boundary inputs (no threshold change to gate) and
  therefore does not have the same defect-surfacing opportunity.
- ADR-0043 (Bessel per-kernel split) — the first Phase 2b ACCEPTED
  perf win; the methodology this sub-slice tried to apply to Airy.
- `src/math/airy.rs:372-406` — `airy_threshold_exponent` (unchanged).
- `src/math/airy.rs:741-752` — `airy_series`'s `(2/3)|x|^{3/2}·log₂e`
  boost (unchanged; flagged for potential future investigation).
- `src/math/airy.rs:498` (`airy_asymptotic_pos`), `:557`
  (`airy_asymptotic_neg`) — `working = target + 64` boost. The most
  obvious follow-up perf lever (reduction to `target + 32`).
- `src/math/airy.rs:326-370` — `airy_zero_value` boundary constants;
  memoisation deferred pending a bench (orthogonal to the threshold).
- `benches/airy_dispatch.rs` — durable bench infrastructure landed by
  this slice; 6 cells straddling `airy_threshold_exponent` at p ∈ {53,
  256}. Run with `cargo bench --bench airy_dispatch --features airy`.
- `Cargo.toml` — `[[bench]] name = "airy_dispatch" required-features =
  ["airy"]` entry.
