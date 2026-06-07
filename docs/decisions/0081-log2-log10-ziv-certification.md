# ADR-0081: Ziv-certify the directed rounding of log2 and log10

- **Status**: accepted
- **Date**: 2026-06-07

## Context

`log2` and `log10` were the only transcendental kernels on the surface
still rounding by a fixed working-precision guard: each composed
`ln(x) / ln(2)` (resp. `ln(10)`) at `target + 64` bits under `NearestEven`
and applied the caller's mode in a single final round, with no Ziv
interval-test certification. The differential harness's own
`NEAREST_EVEN_ROUNDING_MODES` documentation records the failure shape of
that pattern: at an input whose true value lies within the guard band of a
target-grid boundary, the directed final round can land one ULP off the
correctly-rounded result.

The Phase 4 directed-mode bug-hunt (ADR-0079) exercised `log2` and `log10`
under all five modes over a general sweep and an adversarial
near-power-of-two sweep (where the output sits next to an integer grid
point, the directed hard-to-round zone) and found them empirically correct
across the sampled precisions. So this is not a defect fix: no triggering
input was found. It is the certification gap that remains. The cross-check
record corroborates the structural gap: `log2`'s status row reports
`skipped_trace_not_final` on essentially every input, because `log2` left
only its inner `ln`'s Ziv trace and none of its own, the mechanical
signature of a kernel with no outer interval test.

Leaving one kernel family uncertified while the rest of the surface is
per-input certified is exactly the kind of fact that gets rediscovered by
every future reviewer. Certifying it is entropy reduction.

## Decision

Route the positive-finite-normal path of both kernels through the shared
`ziv_round` driver (ADR-0022). The eval closure composes
`ln(x) / ln(const)` at the driver's working precision under `NearestEven`;
the outer interval test then settles the caller's directed mode per input,
doubling its guard until both ends of the evaluation-error interval round
to the same value. Each kernel gets a named calibrated bound in
`ziv_calibration` (`LOG2_ERROR_GUARD`, `LOG10_ERROR_GUARD`), at the
elementary template value, since the composition adds only a constant
materialisation and one divide over `ln`'s own op count.

Two existing behaviours are preserved deliberately:

- **The exact-input dispatch stays ahead of the loop.** A power of two
  (resp. ten) returns its integer logarithm exactly before the Ziv loop
  runs. This is the exact-value-defeats-Ziv guard: the composition of two
  working-precision approximations returns `k + epsilon`, and under a
  directed mode the loop would round that off the exact integer.

- **The non-finite and non-positive inputs keep the composition.** A naive
  Ziv wrap would discard the eval's status, losing `ln`'s `INVALID`
  (from `x < 0`, `NaN`, `-inf`) and `DIV_BY_ZERO` (from `+/-0`) flags.
  Those inputs instead fall to the original composition, which carries the
  correct `+/-inf` or `NaN` value and the flag. Their results are exact,
  so they need no interval-test certification.

The forced `INEXACT` on a finite-normal result (the pf-njs5 over-report:
`log2`/`log10` of a non-exact-set input is irrational) is retained.

## Consequences

- `log2` and `log10` are now per-input correctly-rounded-certified under
  every mode by the interval test, on the same footing as `ln`, `exp`,
  `pow`, and the rest of the transcendental surface, rather than relying
  on a fixed guard that the sampled sweep happened not to break.
- The `pf-tqzz` cross-check will now see `log2`/`log10` leave their own
  outer-Ziv trace, so their status rows should report real passes instead
  of skipping essentially every input.
- The directed-mode bug-hunt lane's `log2`/`log10` checks (general and
  near-power-of-two, all five modes) stay green after the change, and the
  `inexact_fidelity` lane confirms the INEXACT behaviour is unchanged.
- A small cost: the certified path runs the Ziv loop (one or two
  iterations on ordinary inputs) instead of a single fixed-guard round.
  The exact-input and non-finite fast paths are unaffected.

## Related

- Plan: `~/.claude/plans/crystalline-seeking-pancake.md`.
- Beads: `pf-3rtr.3` (this slice); epic `pf-3rtr`.
- Other ADRs: ADR-0022 (the Ziv interval test), ADR-0039 (per-kernel
  `ZIV_ERROR_GUARD` calibration), ADR-0079 (the directed-mode verification
  posture this completes), ADR-0080 (the saturation fix that shipped first
  in this phase).
