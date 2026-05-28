# ADR-0045: Si/Ci asymptotic threshold — already at the mathematical optimum

- **Status**: accepted (documentation-tier)
- **Date**: 2026-05-27

## Context

Phase 2b sub-slice 2b.2.c (`pf-6fvx` work, branch `phase-2b-perf-4`)
targeted the shared Si/Ci asymptotic dispatch threshold
`asymptotic_threshold_exponent` (`src/math/si.rs:158-174`) for
tightening, following the per-kernel-split methodology that landed
ADR-0043 (Bessel) as the first Phase 2b ACCEPTED ADR. The same
methodology immediately found in ADR-0044 (Airy) that the Airy
threshold was already at the mathematical optimum because of the
cube-of-`|x|³` growth. ADR-0044 flagged this as a likely pattern:

> The Airy finding suggests that asymptotic-threshold formulas derived
> from first-principles accuracy laws (rather than the chunky
> conservative cut that Bessel inherited from slice 6o's "not
> perf-tuned without a bench" deferral) tend to be at or near the
> mathematical optimum.

Sub-slice 2b.2.c confirms the pattern for Si/Ci.

## Analysis

Si and Ci share the same threshold function (`ci.rs:43` imports
`super::si::asymptotic_threshold_exponent`). The DLMF 6.12.3 / 6.12.4
asymptotic auxiliaries `f(x)` and `g(x)` have optimal-truncation
error `≈ e^{−|x|}` (the smallest term of the divergent series sits
near `k ≈ |x|`). Reaching `target + G` bits of accuracy needs
`|x|·log₂ e ≥ target + G`, i.e.:

```text
|x| ≥ (p + G) · ln 2 ≈ (p + G) · 0.69315
```

The code at `src/math/si.rs:164-174`:

```rust
pub(super) fn asymptotic_threshold_exponent(target_precision: u32) -> i64 {
    let bits_needed: u64 = u64::from(target_precision) + 32;
    let need: u64 = (bits_needed * 6932).div_ceil(10000);
    // ... smallest e with 2^e ≥ need
}
```

uses the integer rational `6932/10000 = 0.6932` for `ln 2`. The
over-approximation factor is `0.6932 / 0.69315 ≈ 1.00007` — about
**0.007% over the strict accuracy bound**. The formula is at the
mathematical optimum essentially to within rational-approximation
precision.

## The `+32` guard cannot move the threshold at any standard precision

Beyond the already-tight constant, the `+32` guard adds slack
absorbing `SI_ERROR_GUARD = CI_ERROR_GUARD = 24` (Ziv-side
calibrated, `src/math/ziv_calibration.rs:177,180`) plus 8 bits
margin. Shrinking `G` even to zero does not move the threshold at any
of the standard verification precisions:

| precision | `T(G=32)` | `T(G=24)` | `T(G=0)` | moves? |
|----------:|----------:|----------:|---------:|--------|
| 53        | 6         | 6         | 6        | no |
| 113       | 7         | 7         | 7        | no |
| 256       | 8         | 8         | 8        | no |
| 1024      | 10        | 10        | 10       | no |

The binary-exponent quantization is the structural lock: at `p=256`,
`(256+G)·0.6932` must stay below 128 for `T(256)` to drop from 8 to 7,
requiring `G ≤ -71.3` — impossible. At `p=1024`, even `G=0` gives
`1024·0.6932 ≈ 710` and the smallest `e` with `2^e ≥ 710` is still 10.

In words: the formula was tight in the constant (rational
approximation to `ln 2`), and the precision range we care about lands
each `(p+G)·ln 2` value comfortably inside its binary octave. There
is no slack to harvest by tightening `G`.

## Decision

**Close sub-slice 2b.2.c documentation-tier per ADR-0040 GATE A
and ADR-0044 precedent.** No code change to the threshold function
ships.

**No new bench infrastructure lands.** Unlike the Bessel and Airy
cases where a `benches/*_dispatch.rs` file was useful as durable
infrastructure for potential future re-measurement, Si/Ci's threshold
is so structurally locked that there is no realistic future
re-measurement scenario the bench would help. The existing
`tests/differential_si.rs::SI_TABLE` and `tests/differential_ci.rs::CI_TABLE`
already span both regimes (small-`|x|` log-series and large-`|x|`
asymptotic, with table entries at `x ∈ {½, 1, 2, 3, 5, 8, 13, 50, 100,
500, 2000, 5000}` for Si and the analogous Ci coverage). If future Si/Ci
perf work targets something other than the threshold (the
`si_series`/`ci_series` cancellation guard `≈ |x|·log₂ e` at
`src/math/si.rs:184-194` and `src/math/ci.rs:170-180`, or the shared
`si_ci_f`/`si_ci_g` auxiliaries themselves), a fresh bench can be
built then.

## Consequences

**No behavioural change to `Si(x)` or `Ci(x)`.** The threshold and
both kernels run unchanged. Verification (the in-module tests,
`differential_si`, `differential_ci`) is unaffected.

**Pattern confirmed for first-principles asymptotic thresholds.** Two
of the three asymptotic-threshold sub-slices examined since ADR-0043
(Airy, Si/Ci) have found the existing formula already at the
mathematical optimum. The remaining sub-slice 2b.2.d (erf/erfc) uses
`4^{e_x} ≥ ⌈(p+32)·16/23⌉` at `src/math/erf.rs:181-192` — the same
shape (integer rational approximation to a tight accuracy law). It is
very likely also already-optimal; the analytical vet for 2b.2.d
should examine whether `(p+G)·16/23` (a tight approximation to
`(p+G)/log₂ e`) similarly locks the threshold across the
TRANSCENDENTAL_PRECISIONS grid.

**ADR ledger update:**

| ADR | Outcome | Slice |
|----:|---------|-------|
| 0037 | REJECTED | SmallVec mantissa swap |
| 0040 | doc-tier (GATE A) | FFT measurement |
| 0041 | REJECTED | Spouge precision-pegging |
| 0042 | ACCEPTED | pf-1axr trig + bessel_y recurrence root fix |
| 0043 | ACCEPTED | Bessel per-kernel threshold split |
| 0044 | doc-tier | Airy threshold — already optimal |
| **0045** | **doc-tier** | **Si/Ci threshold — already optimal** |

Sub-slices remaining on `pf-6fvx`: 2b.2.d (erf/erfc — likely doc-tier
per the pattern above) and 2b.3 (Bessel Miller-recurrence depth — the
remaining sub-area expected to be a real code-change candidate).

## Related

- ADR-0040 (FFT GATE A) — doc-tier closure precedent for
  "measurement confirms no perf win available".
- ADR-0043 (Bessel per-kernel split) — the methodology this slice
  attempted to apply.
- ADR-0044 (Airy threshold already optimal) — the precedent for the
  analytical-vet-only doc-tier closure when a threshold is at the
  mathematical optimum. ADR-0044 explicitly anticipated this finding
  for Si/Ci and erf/erfc.
- `src/math/si.rs:158-174` — `asymptotic_threshold_exponent` (the
  threshold function shared by Si and Ci; unchanged).
- `src/math/si.rs:135`, `src/math/ci.rs:144` — Si/Ci dispatch sites
  using the shared threshold.
- `src/math/ci.rs:43` — `use super::si::{asymptotic_threshold_exponent, …}`
  showing Ci's reuse.
- `tests/differential_si.rs::SI_TABLE`, `tests/differential_ci.rs::CI_TABLE`
  — existing rich boundary coverage spanning both regimes; sufficient
  for any future Si/Ci perf change verification without new bench
  infrastructure.
- `src/math/ziv_calibration.rs:177,180` — `SI_ERROR_GUARD` and
  `CI_ERROR_GUARD` both at `DEFAULT_ERROR_GUARD = 24` per
  Phase 1g (ADR-0039); the `+32` in the threshold formula matches
  with 8 bits margin.
