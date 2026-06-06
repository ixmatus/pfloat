# ADR-0066: Defensive INEXACT guard on the gamma family and zeta

Status: Accepted (2026-06-05)

## Context

ADR-0064 extended the INEXACT-flag fix to the proven-transcendence
special functions. The gamma family (gamma, lgamma, digamma, beta) and
zeta were deferred because their values at dyadic arguments are believed
irrational but with named open problems: the irrationality of Euler's
constant `γ = −ψ(1)`, of `ζ(5)`, and the transcendence of `Γ` at dyadic
denominators ≥ 8.

Before forcing INEXACT on them, a feasibility sweep (a flag-and-value
comparison against a p=200 oracle, roughly 6200 input × mode checks) was
run. It found:

- gamma, lgamma, digamma, zeta: zero defects (0 under-reports, 0
  over-reports, 0 value errors). They already report INEXACT correctly on
  every finite-normal fall-through, because the under-report phenomenon
  does not arise: their exact dyadic outputs are all dispatched (gamma's
  factorials, zeta's `−1/2` and trivial zeros), digamma has no dyadic
  outputs at all, and zeta's only collapse case (`ζ(s) → 1` for large `s`)
  is a true INEXACT since `ζ(s) > 1` strictly for `s > 1`.
- beta: a real correctness bug, fixed separately (ADR-0065).

## Decision

Add the ADR-0060 force-INEXACT to the finite-normal Ziv fall-through of
gamma, lgamma, digamma, zeta, and beta's non-integer fall-through, guarded
on `Class::Normal` so poles, exact non-finite limits, and dispatched exact
outputs keep their status. For beta the force sits only on the Ziv
fall-through; the positive-integer and case-4 construct-and-check
dispatches (ADR-0065) return earlier with the correct flag and are
untouched, so a dyadic beta output is never forced.

This is a structural guard, not a fix. The sweep showed every one of these
paths already flags INEXACT, so the force is presently a no-op; it hardens
the invariant against a future regression (a new collapse-prone code path
that lands a transcendental result on a grid value with the flag cleared),
the same role the force plays for the proven families in ADR-0064.

## Soundness: a conditional guarantee, stated honestly

For the proven families (ADR-0064) the force is unconditionally sound:
outside the dispatched set the result is provably irrational. For the
gamma family and zeta that proof is not available. The force is sound iff
every dyadic-output input is dispatched first, which reduces to the
fall-through values being non-dyadic. They are believed irrational, but
the strongest statements have open cases:

- digamma: `ψ(1) = −γ`; the irrationality of the Euler-Mascheroni constant
  `γ` is an open problem.
- zeta: `ζ(5)` (and the irrationality of `ζ` at a general non-special
  argument) is open. `ζ(3)` is irrational (Apéry), but `ζ(5)`, `ζ(7)`, …
  are not individually settled.
- gamma, lgamma, beta: the transcendence of `Γ` at dyadic rationals with
  denominator ≥ 8 (for example `Γ(1/8)`) is not established.

So a false positive, forcing INEXACT on a genuinely exact result, would
require one of these believed-irrational values to be exactly dyadic at
the target precision, which would resolve a famous open problem in the
most surprising way. The realized risk is nil: the feasibility sweep found
no such input, and for any concrete precision the computed value's bits
are checkable and do not terminate. A false negative (the guard failing to
add INEXACT) cannot happen here, since the sweep shows the underlying Ziv
path already sets it. We accept the conditional worst-case soundness and
record the open problems explicitly; if any is ever resolved against the
believed answer, the affected force is the one to revisit.

## Consequences

- INEXACT fidelity across the gamma family and zeta is now structurally
  guaranteed, matching the proven families. Value-preserving (metadata
  only): the differential-mpfr lanes for gamma, lgamma, digamma, zeta, and
  beta still match MPFR. Point tests in `tests/inexact_fidelity.rs`
  (`cfg(specials)` / `cfg(zeta)`) assert the dispatched exact outputs clear
  and the transcendental / rational-non-dyadic results set INEXACT.
- This completes pf-umlm and, with ADR-0064 and ADR-0065, the whole
  pf-cjnk special-function INEXACT program.

## Related

- pf-umlm (this work); ADR-0064 (the proven families); ADR-0065 (the beta
  bugfix, which makes this force sound for beta); ADR-0060 and ADR-0063
  (the original pattern). On the open problems: the irrationality of `γ`;
  Apéry's theorem that `ζ(3)` is irrational and the open status of
  `ζ(2k+1)` for `k ≥ 2`; the transcendence of `Γ` at small-denominator
  rationals (Chudnovsky) and the open status at general dyadic
  denominators.
