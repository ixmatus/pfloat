# ADR-0017: Transcendental constants computed on the fly via AGM

Status: accepted (slice 7b shipped)

## Context

ADR-0014's slice 6h status update flagged "transcendental
precisions capped at 256 bits" as the second of three structural
limitations the MPFR differential lane exposed. The diagnosis
attributed the cap to pfloat's hardcoded 1024-bit reduction
constants (`LN2_LIMBS_1024`, `PI_LIMBS_1024`,
`TWO_OVER_SQRT_PI_LIMBS_1024`, `LN_2PI_LIMBS_1024`) running out of
faithful bits above the 960-bit working-precision threshold that
the 64-bit guard imposed.

Slice 7b investigated the failure mode in preparation for lifting
the cap and surfaced a sharper finding: `LN2_LIMBS_1024` is
encoded incorrectly past bit ~450. The hardcoded mantissa agrees
with the mathematical value of `ln(2)` to about the first 450
bits, then diverges. Authoritative reference values (Brent's
1976 paper, MPFR's reference output, and an independently derived
1100-decimal-digit constant) all match each other; the
pfloat-embedded constant matches only the top of those values.
The slice-6h diagnosis "constants run out of bits above 960" was
misattributing the symptom: the constant runs into a real
encoding defect at ~450 bits, well below the diagnosed boundary.
`PI_LIMBS_1024` was cross-checked and is correct to >1000 bits.

The slice's primary task — adding on-the-fly AGM-based
computation so the transcendental kernels can produce results at
any working precision — is unchanged; the secondary task of
fixing the broken constant became the immediate-priority half of
the slice.

## Decision

1. **Add `src/math/agm_constants.rs`.** The module exports a
   minimal high-precision constant kit:

   - `pi_via_agm(prec)` uses Brent–Salamin's iteration. `a_0 = 1`,
     `b_0 = 1/√2`, `t_0 = 1/4`, `p_0 = 1`; iterate
     `a_{n+1} = (a_n + b_n)/2`, `b_{n+1} = √(a_n · b_n)`,
     `t_{n+1} = t_n − p_n · (a_n − a_{n+1})²`,
     `p_{n+1} = 2 · p_n`; the limit is
     `π = (a_∞ + b_∞)² / (4 t_∞)`. Quadratic convergence; `O(log
     p)` iterations.
   - `ln_2_via_atanh(prec)` uses the identity
     `ln(2) = 2 · atanh(1/3)`. The atanh series converges
     `log₂(9) ≈ 3.17` bits per term; `O(p)` terms.
   - `ln_10_via_atanh(prec)` uses `ln(10) = 3·ln(2) + 2·atanh(1/9)`
     which converges `log₂(81) ≈ 6.34` bits per term.
   - `two_over_pi_via_agm`, `two_over_sqrt_pi_via_agm`, and
     `ln_2pi_via_agm` compose `π` and `ln(2)` through the existing
     arithmetic kernels.

   The module is gated `cfg(feature = "exp-log")` and the
   `exp-log = ["big", "agm"]` widening (ADR-0015 anticipated this)
   makes the AGM kernel available to any code that compiles in
   transcendentals.

2. **Per-constant dispatch in `src/math/mod.rs`.** Each
   `*_at(prec)` accessor picks the table fast-path or the AGM
   slow-path based on the table's verified-correct precision:

   | Accessor | Table cap | Above the cap |
   | --- | --- | --- |
   | `ln_2_at` | 256 (per `LN2_TABLE_PRECISION_CAP`) | `ln_2_via_atanh` |
   | `pi_at` | 1024 | `pi_via_agm` |
   | `pi_over_2_at` | 1024 (exponent shift on the table) | `pi_via_agm`, then exponent decrement |
   | `two_over_pi_at` | 4096 | `two_over_pi_via_agm` |
   | `two_over_sqrt_pi_at` | 1024 | `two_over_sqrt_pi_via_agm` |
   | `ln_2pi_at` | 1024 | `ln_2pi_via_agm` |
   | `ln_10_at` | always AGM | `ln_10_via_atanh` |

   `LN2_TABLE_PRECISION_CAP` is conservatively set to 256 (well
   below the ~450-bit boundary where the hardcoded mantissa starts
   diverging). Regenerating the hardcoded table to faithfully
   represent ln(2) to 1024 bits is a Phase 7 polish item; the
   current AGM path is correct at every precision.

3. **Lift the `.min(1024)` working-precision cap in all 22
   transcendental kernels.** The cap was a workaround for the
   constants' precision ceiling. With the constants now correct at
   any precision, the cap can come off cleanly. Each kernel's
   `working_prec = target_precision.saturating_add(N)` now scales
   with the caller's target instead of saturating at 1024.

4. **Hold the differential lane at `TRANSCENDENTAL_PRECISIONS =
   [53, 113, 256]`.** The constants are now correct at any
   precision, but each call recomputes them; at `p ≥ 1024` an exp
   call costs on the order of seconds because `ln(2)` runs ~360
   atanh terms at working precision 1088. Lifting the lane to
   include `p = 1024` makes the per-op sweep prohibitive
   (≥30 minutes per op). The lift waits on memoization (see
   Consequences).

## Consequences

- The pre-slice-7b `BigFloat::ln`, `exp`, `pow`, and `trig`
  transcendentals at any target precision above ~450 bits were
  silently using a faulty `ln(2)` for argument reduction. Slice
  7b's fix retroactively closes that correctness gap. No public
  API change; the affected output values now match the
  mathematical truth.
- The slice-6h status update in ADR-0014 attributed
  `differential_exp` at `p > 256` divergence to the constants'
  bit budget. The actual root cause is the broken `LN2_LIMBS_1024`
  encoding. ADR-0014's status update should be revisited once
  this slice merges.
- The `LN2_TABLE_PRECISION_CAP = 256` is a conservative
  threshold. A future polish slice can regenerate the table
  faithfully to 1024 bits and lift the threshold, which gains
  back the cheap fast-path for the common `p ≤ 1024` range. The
  conservative-and-correct path ships first.
- The transcendental kernels now correctly produce results at
  any target precision, but the per-call cost at high precision
  is dominated by AGM constant recomputation. A `thread_local!`
  memoization slice (Phase 7 follow-up) lifts this restriction:
  the most-recent `(prec, BigFloat)` for each constant is
  cached, and consecutive calls at the same precision (the
  common case in differential testing and hot loops) hit O(1).
  After memoization, `TRANSCENDENTAL_PRECISIONS` can include
  `1024` without the time penalty.
- The `agm` feature is now transitively present whenever
  `exp-log` is enabled. Embedded users targeting `--features=big,
  exp-log` get the AGM kernel even if they would not have opted
  in directly; the AGM module is small (under 300 lines of
  source) and adds negligible binary footprint.

## References

- Plan: `let-s-review-the-backlog-vast-harbor.md` — slice 7b.
- ADR-0014 — MPFR differential gating. The slice-6h status
  update flagged the cap that this slice partially closes; the
  remaining piece (caching, differential lift) is the open work.
- ADR-0015 — AGM kernel formulation. Anticipated this slice's
  widening of `exp-log` to depend on `agm`.
- Brent, R. P. *Multiple-precision zero-finding methods and the
  complexity of elementary function evaluation.* In *Analytic
  Computational Complexity* (J. F. Traub, ed.), pp. 151–176,
  Academic Press, 1976. — The Brent–Salamin iteration's
  derivation and convergence analysis.
- Borwein, J. M. and Borwein, P. B. *Pi and the AGM.* Wiley,
  1987. — The reference text on AGM-based constant computation;
  the `ln(x) ≈ π / (2 · AGM(1, 4/x))` family and its variants
  are catalogued in Chapter 11.
