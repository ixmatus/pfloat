# ADR-0046: erf/erfc asymptotic threshold — already at the mathematical optimum

- **Status**: accepted (documentation-tier)
- **Date**: 2026-05-27

## Context

Phase 2b sub-slice 2b.2.d (`pf-6fvx` work, batched onto branch
`phase-2b-perf-4` alongside 2b.2.c) targeted the shared erf/erfc
asymptotic dispatch threshold `asymptotic_threshold_exponent`
(`src/math/erf.rs:169-192`, imported by `erfc.rs:36`) for tightening
following the per-kernel-split methodology of ADR-0043 (Bessel). ADR-0044
(Airy) and ADR-0045 (Si/Ci) established the pattern that
first-principles asymptotic threshold formulas land at or near the
mathematical optimum; 2b.2.d closes the third such finding on the
prior-existing erf/erfc threshold.

## Analysis

erf and erfc share the same threshold (`erfc.rs:36` imports via
`use super::erf::asymptotic_threshold_exponent`). The DLMF §7.12
asymptotic expansion of `erfc(x)` has smallest-term truncation error
`≈ e^{−x²}` (the Gaussian decay of the prefactor times the divergent
series sums to a value proportional to `e^{−x²}`). Reaching `target + G`
bits of accuracy needs `x²·log₂ e ≥ target + G`, i.e.:

```text
|x|² ≥ (p + G) / log₂ e = (p + G) · ln 2 ≈ (p + G) · 0.69315
```

The code at `src/math/erf.rs:181-192` uses `|x|² ≥ 4^{e_x}` and the
integer rational `16/23 ≈ 0.69565` for `ln 2 ≈ 0.69315`. The
over-approximation factor is `0.69565 / 0.69315 ≈ 1.0036` — about
**0.36% over the strict accuracy bound**. The doc-comment names this
explicitly: *"use `log₂ e ≈ 23/16`"*. The formula is essentially at the
mathematical optimum modulo rational-approximation precision.

## The `+32` guard cannot move the threshold at any standard precision

The `+32` budget absorbs `ERF_ERROR_GUARD = ERFC_ERROR_GUARD = 24`
(Ziv-side calibrated per Phase 1g, `src/math/ziv_calibration.rs:168,171`)
plus 8 bits of margin. Shrinking `G` even to zero does not move the
threshold at any TRANSCENDENTAL_PRECISIONS value:

| precision | `T(G=32)` | `T(G=0)` | `|x|` cut | moves? |
|----------:|----------:|---------:|----------:|--------|
| 53        | 3         | 3        | ≥ 8       | no |
| 113       | 4         | 4        | ≥ 16      | no |
| 256       | 4         | 4        | ≥ 16      | no |
| 1024      | 5         | 5        | ≥ 32      | no |

The binary-exponent quantization is even more stringent than Si/Ci's:
each exponent step changes `4^{e_x}` by 4× (vs Si/Ci's `2^{e_x}` 2×
steps). At `p=1024`, dropping `T` from 5 to 4 would require
`need ≤ 256`, i.e. `(1024+G)·16/23 ≤ 256` → `1024+G ≤ 368` →
`G ≤ -656` — nonsensical. At `p=256`, similarly impossible.

## Decision

**Close sub-slice 2b.2.d documentation-tier per the ADR-0040 GATE A,
ADR-0044, and ADR-0045 precedent.** No code change to the threshold
function ships.

**No new bench infrastructure lands.** The same reasoning as ADR-0045:
the threshold is structurally locked, and the existing
`tests/differential_erf.rs` / `tests/differential_erfc.rs` already
exercise both regimes (small-`|x|` Maclaurin via `next_i64_in(state,
-10, 10)`, and `|x| > 30` asymptotic via the larger fixed test points
when present). If future erf/erfc perf work targets something other
than the threshold (the `erf_maclaurin` cancellation guard or the
`erfc_asymptotic` working-precision boost), a fresh bench can be built
then.

**Honest caveat: the existing differential test grid does not exercise
the asymptotic dispatch at high precision.** As noted in
project-state-pointer memory `feedback_precision_gated_verification_surface`
and the ADR-0044 §"What lands" section, `differential_erf` and
`differential_erfc` use a random sweep capped at `|x|=10`, which sits
below the asymptotic threshold at every precision (T(53)=3 ⇒ |x|≥8 is
the lowest cut). The test grid is adequate to verify the *current*
behaviour (which uses the asymptotic at |x|≥8, well inside the
random sweep at p=53), but if a future investigator changes the
threshold formula they should add explicit boundary inputs analogous to
the Bessel `(1025, 1)` / `(2049, 1)` cases per the sub-slice 2b.2.a
methodology. This slice ships no threshold change, so no boundary
inputs are added; the gap is flagged here for the next investigator.

## Consequences

**No behavioural change to `erf(x)` or `erfc(x)`.** The threshold and
both kernels run unchanged. Verification (in-module tests,
`differential_erf`, `differential_erfc`) is unaffected.

**Pattern fully confirmed for first-principles asymptotic thresholds.**
All three remaining 2b.2 sub-slices found the existing formula already
at the mathematical optimum:

| sub-slice | function family | finding |
|-----------|-----------------|---------|
| 2b.2.b (ADR-0044) | Airy | ~1.5% over strict; cube-of-`|x|³` growth locks T |
| 2b.2.c (ADR-0045) | Si, Ci | ~0.007% over strict; `2^e` quantization locks T |
| **2b.2.d (this ADR)** | **erf, erfc** | **~0.36% over strict; `4^e` quantization locks T** |

Only Bessel's `bessel_j_threshold` had room (the "deliberately
conservative" cut inherited from slice 6o's "not perf-tuned without a
bench" deferral, which ADR-0043 split per-kernel for the 89-97% wins
on Y/I/K log-series flip cells). The general lesson, durable for any
future asymptotic-threshold investigation: a formula derived directly
from the kernel's accuracy law tends to land at the optimum; only
formulas explicitly carrying conservative slack (the "perf-tuning
deferred" pattern) have room for a sub-slice to harvest.

**ADR ledger after 2b.2.d:**

| ADR | Outcome | Slice |
|----:|---------|-------|
| 0037 | REJECTED | SmallVec mantissa swap |
| 0040 | doc-tier (GATE A) | FFT measurement |
| 0041 | REJECTED | Spouge precision-pegging |
| 0042 | ACCEPTED | pf-1axr trig + bessel_y recurrence root fix |
| 0043 | ACCEPTED | Bessel per-kernel threshold split |
| 0044 | doc-tier | Airy threshold — already optimal |
| 0045 | doc-tier | Si/Ci threshold — already optimal |
| **0046** | **doc-tier** | **erf/erfc threshold — already optimal** |

The only Phase 2b sub-slice remaining on `pf-6fvx` is **2b.3 Bessel
Miller-recurrence depth** — the last sub-area expected to be a real
code-change perf sub-slice. After 2b.3 closes (whatever outcome),
`pf-6fvx` closes, and Phase 2b is done.

## Related

- ADR-0040 (FFT GATE A) — first doc-tier closure precedent.
- ADR-0043 (Bessel per-kernel split) — the methodology this slice
  attempted to apply; found inapplicable because the formula is
  already at the optimum.
- ADR-0044 (Airy threshold already optimal), ADR-0045 (Si/Ci threshold
  already optimal) — the two prior doc-tier ADRs establishing the
  pattern.
- `src/math/erf.rs:169-192` — `asymptotic_threshold_exponent` (the
  threshold function shared by erf and erfc; unchanged). The
  doc-comment names `log₂ e ≈ 23/16` explicitly as the integer
  rational, which is the strongest evidence the formula was derived
  with mathematical-optimum intent.
- `src/math/erf.rs:132`, `src/math/erfc.rs:135` — dispatch sites.
- `src/math/erfc.rs:36` — `use super::erf::{asymptotic_threshold_exponent, …}`
  showing erfc's reuse.
- `src/math/ziv_calibration.rs:168,171` — `ERF_ERROR_GUARD` and
  `ERFC_ERROR_GUARD` both at `DEFAULT_ERROR_GUARD = 24` per Phase 1g
  (ADR-0039); the `+32` in the threshold formula matches with 8 bits
  margin.
- `tests/differential_erf.rs`, `tests/differential_erfc.rs` — existing
  differential tests using `next_i64_in(state, -10, 10)`. Flagged as
  inadequate for verifying a *changed* threshold formula in the §"Honest
  caveat" section; sufficient for verifying the *current* unchanged
  formula.
