# ADR-0065: Beta exact-input construct-and-check (correctness bugfix)

Status: Accepted (2026-06-05)

## Context

The pf-cjnk feasibility sweep (a flag-and-value comparison against a
high-precision oracle across the gamma family and zeta, run before
extending the INEXACT-flag fix to those families) found gamma, lgamma,
digamma, and zeta already correct on the flag everywhere it checked, but
surfaced a pre-existing correctness defect in `beta`.

For positive-integer arguments, `β(a,b) = (a−1)!(b−1)!/(a+b−1)!` is
rational, and the kernel computed it through the `exp(lgamma)`
composition under the Ziv driver. At the exactly-dyadic outputs
(`β(1,1) = 1`, `β(1,2ᵏ) = β(2ᵏ,1) = 2⁻ᵏ`) the composition's rounding
noise defeats the Ziv interval test: the same exact-value-defeats-Ziv
failure mode that `try_gamma_pos_integer_exact` was added to cure for
gamma. The sweep counted 78 INEXACT over-reports and, worse, 30
directed-rounding value errors. For example `β(1,2)` under TowardPositive
returned the successor of `0.5`, and `β(1,8)` under TowardZero /
TowardNegative the predecessor of `0.125`. Returning a value one ULP off
an exactly-representable true value violates the correctly-rounded
contract, so this is a correctness bug in the released v1.1.0 beta, not
only a flag-fidelity nicety.

## Decision

Add a pre-Ziv exact dispatch for the rational beta cases, mirroring
gamma's factorial fast path. Build the exact rational and divide once at
the target precision; the division's INEXACT flag is the exact dyadicity
verdict, so a dyadic output returns the exact value with INEXACT clear
and a non-dyadic output the correctly-rounded value with INEXACT set.

- **Positive-integer arguments (cases 1 and 2).** With `(s, l) =
  (min, max)`, `β(a,b) = (s−1)! / ∏_{j=0}^{s−1}(l+j)`. The numerator and
  denominator are built as exact integers; the denominator accumulator
  seeds from the input `l` itself, so the `s = 1` family (`β(1, b) = 1/b`,
  which holds every dyadic case `β(1, 2ᵏ) = 2⁻ᵏ`) needs no large
  factorial and works for an arbitrarily large `b`.
  `num.div_round(&den, target, mode)` returns the result.
- **Case 4 (pole cancellation).** `B(−n, m) = (−1)^m (m−1)!(n−m)!/n!` is
  likewise rational and sometimes dyadic (`B(−1,1) = −1`); it gets the
  same construct-and-check. The sign is applied to the numerator before
  the division so the directed-rounding boundary lands on the signed
  value.

Both dispatches return `None` (falling back to the existing `exp(lgamma)`
Ziv path) when a factorial overflows the exact-integer build precision or
factor cap. An output too large to build that way is far below an ULP and
non-dyadic, so the Ziv rounding is correct there; the case-4 fallback
also preserves the O(1) lgamma-factorial form that survives a
caller-supplied huge order (ADR-0030).

## Consequences

- The 78 over-reports and 30 value errors go to zero. `β` is now
  correctly rounded with a correct INEXACT flag at every positive-integer
  and case-4 input under all five rounding modes. A regression test
  (`tests/beta_exact_fidelity.rs`) asserts the dyadic cases exact, the
  non-dyadic cases INEXACT, and runs an oracle sweep (p=53 against the
  kernel's own p=200 result rounded down) over the integer and case-4
  grids with zero mismatches.
- Value-correcting, not value-preserving: the dispatched values change
  for the previously mis-rounded directed-mode cases, toward the correct
  value. The `differential_beta` MPFR lane still passes, because the new
  values match MPFR where the old directed-mode ones were a ULP off.
- gamma, lgamma, digamma, and zeta were verified defect-free by the same
  sweep and are untouched here. The companion defensive INEXACT guard on
  their finite-normal fall-throughs (currently a no-op, with a worst-case
  soundness that is conditional on open irrationality problems) is
  ADR-0066.

## Related

- pf-umlm (this work); pf-cjnk and ADR-0064 (the feasibility sweep that
  found this defect); ADR-0030 (the beta case classification and the
  case-4 closed form); gamma's `try_gamma_pos_integer_exact` positive
  integer fast path and the exact-value-defeats-Ziv lesson it encodes.
