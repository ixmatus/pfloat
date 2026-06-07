# ADR-0080: directed-mode rounding for the asymptote-saturating kernels

- **Status**: accepted
- **Date**: 2026-06-07

## Context

The Phase 4 directed-mode bug-hunt (ADR-0079, the lane in
`tests/directed_mode_bug_hunt.rs`) found a real, shipped, range-affecting
directed-rounding defect shared by five kernels that approach a nonzero,
on-grid asymptotic limit they never reach:

- `tanh(x) → ±1`,
- `erf(x) → ±1`,
- `erfc(x) → 2` as `x → −∞`,
- `expm1(x) → −1` as `x → −∞`,
- `zeta(s) → 1` from above as `s → +∞`.

Each computes through the Ziv interval-test driver (ADR-0022). For `|x|`
large enough that the residual to the limit underflows every working
precision the Ziv loop reaches (it doubles its guard up to
`target + ZIV_GUARD_CAP`), the kernel's evaluator returns the exact limit
at every iteration. The interval test can then never separate the two
sides of the limit, so the loop exhausts `ZIV_MAX_ITERS` and the
cap-fallback returns the on-grid limit itself. That value is correctly
rounded under `NearestEven`, `NearestAway`, and the directed mode that
rounds away from the interior, but it is one ULP wrong under `TowardZero`
and the inward directed mode, where the correctly-rounded result is the
interior neighbour (`±(1 − ulp)`, `2 − ulp`, `1 + ulp`). The defect spans
an unbounded input tail, not a measure-zero set: `tanh(8066)`,
`erf(−4.5e8)`, `erfc(−30)`, `expm1(−1500)`, and `zeta(5000)` all reproduce
it. It violates the surface's correct-rounding-under-every-mode claim, and
it is the kind of one-sided directed error `pfloat-ball`'s directed-pair
enclosure is most exposed to.

This is the saturation analogue of the documented cancellation-to-zero
hole (the tanh tiny-x case, ADR-0050): in both, the evaluator collapses to
an exact representable value while the true value is strictly to one side,
and the Ziv half-width, being relative to that value, cannot see the
residual. Kernels whose limit is off-grid (`atan → π/2`,
`log2`/`log10 → ±∞`) are unaffected, because the final directed round of
an off-grid value moves off it on the correct side.

## Decision

Each affected kernel gains a large-`|x|` short-circuit, mirroring the
tiny-x short-circuit already present in `tanh` and `expm1`. Once the
argument's binary exponent crosses a precision-dependent threshold, the
kernel returns `crate::rounding::round_with_infinitesimal(limit, sign,
subtracts_magnitude, target, mode)` rather than entering the Ziv loop. The
infinitesimal carries the residual's direction, so the result is the limit
under the nearest and outward modes and the interior neighbour under the
inward directed mode and `TowardZero`, exactly as correct rounding
requires.

The threshold is the smallest argument exponent for which the residual to
the limit is below half a `target`-ULP, so the infinitesimal model is
valid. Two helpers in `src/math/mod.rs` cover the two decay rates on the
surface:

- `saturation_threshold_exponent` for residuals decaying like `2^-|x|`
  (`zeta`), `e^-|x|` (`expm1`), or `e^-2|x|` (`tanh`). For these the span
  between where the residual first drops below half a ULP and where it
  underflows the guard cap is wide, so a single linear threshold lands
  safely inside it.
- `saturation_threshold_exponent_gaussian` for residuals decaying like
  `e^-x^2` (`erf`, `erfc`), which saturate at a much smaller `|x|`. The
  threshold exponent is halved accordingly.

Per kernel the limit and residual direction are: `tanh`/`erf` round `±1`
with a magnitude-shrinking infinitesimal; `erfc` rounds `2` from below;
`expm1` rounds `−1` from above (magnitude shrinking); `zeta` rounds `1`
from above (magnitude growing). The forced `INEXACT` each kernel already
applies to a finite-normal transcendental result (ADR-0063) is preserved.

## Consequences

- The five kernels are correctly rounded under all five modes across the
  saturation tail, verified by `saturation_limit_reproducers` (the five
  reproducers) and `saturation_directed_sweep` (the boundary and tail for
  the rug-orable kernels, all modes), with no regression in the 826 library
  unit tests, the `erf`/`erfc`/`zeta` differential lanes, or the small
  -argument lane.
- The short-circuit also makes the common large-`|x|` path cheaper: it
  skips the Ziv loop's five capped iterations for a single infinitesimal
  round.
- **Honest range caveat.** The threshold is keyed to the argument's binary
  exponent, a power-of-two-coarse quantity. For the linear-decay kernels
  the valid window is wide and the threshold is correct across the whole
  supported precision range. For the Gaussian-decay kernels (`erf`,
  `erfc`) the window between validity and the guard cap narrows as the
  target precision grows; the threshold is verified correct through the
  precisions the directed-mode lanes exercise and well beyond, but at
  extreme target precision (past roughly two thousand bits) a narrow band
  of `|x|` just beyond the cap could still reach the Ziv fallback. A
  residual-probe formulation would remove the caveat at the cost of an
  extra evaluation on every call; the threshold form is preferred for the
  precision range pfloat's directed-mode claim covers, and the caveat is
  recorded here rather than hidden.

## Related

- Plan: `~/.claude/plans/crystalline-seeking-pancake.md`.
- Beads: `pf-3rtr.11` (this fix), discovered from `pf-3rtr.2` (the
  bug-hunt); epic `pf-3rtr`.
- Commits: the bug-hunt lane and reproducer landed first; this commit
  carries the kernel fixes and activates the reproducer.
- Other ADRs: ADR-0079 (the directed-mode verification posture this fix
  serves), ADR-0050 (the tiny-x cancellation-to-zero analogue), ADR-0022
  (the Ziv interval test whose cap the saturation tail exhausts), ADR-0063
  (the forced `INEXACT` for transcendental results preserved here).
