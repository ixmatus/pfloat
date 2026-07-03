# ADR-0118: gamma/beta exact-dispatch precision — walk factorials at build precision, not at target

- **Status**: accepted
- **Date**: 2026-07-03

## Context

Three exact-dispatch defects in `gamma`/`beta` from the 2026-06 review
(epic pf-8iji). All share a root cause: an *exact* integer/rational
walk performed at the **caller's target (or operand) precision**, where
the intermediate integers round.

1. **pf-8scf — `gamma` positive-integer walk at target precision.**
   `try_gamma_pos_integer_exact` accumulated `1·2·…·(n−1)` **and the loop
   counter** at `target_precision`. At a tiny target the counter itself
   rounds: for `gamma(4)@p1` the increment `2 → 3` rounded to `4`, so the
   loop saw `k ≥ x` one step early and returned `Some(2)` — a wrong value
   with OK, where `Γ(4) = 6` is not representable at p1.

2. **pf-ihsp — `beta` case-4 closed form at operand precision.**
   `try_beta_case4_exact` built `n` and `m` from the operands and computed
   `n − m` (and the factorial arguments) at OPERAND precision: for
   `B(−128@p1, 1@p1)` the subtraction `128 − 1 = 127` rounded to `128` at
   p1, giving `(128)!/128! · 0! = 1` and returning `−1` with OK where
   `−1/128` is due.

3. **pf-k5ll — `beta` unit case missed the non-integer dyadic.**
   `B(1, b) = 1/b` exactly (`Γ(1)=1`, `Γ(1+y)=y·Γ(y)`). The
   positive-integer dispatch covered integer `b` (`B(1, 8) = 1/8`) but
   `B(1, 1/8)` — `b = 2^-3`, not an integer — fell through to the
   exp(lgamma) Ziv path and returned `8` with INEXACT over-reported (a
   directed-mode ulp risk on an exact dyadic output).

## Decision

1. pf-8scf: walk at a **build precision** `target + 64` (which holds every
   factorial that could fit at target — a factorial that fits at target is
   `≤ target` bits — and the counter, bounded by the 512-iteration cap),
   so the walk is exact; a product past the build width means the
   factorial exceeds it (and hence target), so bail. At the end,
   round the exact factorial to `target_precision`; return it only if that
   round is exact, else `None` (the Ziv envelope then correctly rounds the
   non-representable factorial — e.g. `Γ(4)@p1` goes to Ziv and truncates
   6 to 4 under TowardZero).

2. pf-ihsp: lift `n` and `m` to `BETA_EXACT_BUILD_PREC` (4096) before the
   integer arithmetic; a non-representable operand rounds inexactly there
   and bails to the Ziv path (whose O(1) lgamma-factorial form survives a
   huge argument).

3. pf-k5ll: add a construct-and-check reciprocal dispatch —
   `B(1, b) = 1/b` (and `B(a, 1) = 1/a`) built by one division at target
   under the caller's mode. The reflection/pole cases upstream have peeled
   off every argument where `1+y` is a Γ pole, so `1/y` is the finite
   value here, and the division's INEXACT flag is the exact-dyadicity
   verdict: a power-of-two operand returns the exact reciprocal with OK,
   any other operand the correctly-rounded inexact reciprocal — still the
   *exact value of B*, so this supersedes the Ziv path for the unit case
   rather than merely fast-pathing it. Extends the ADR-0065
   construct-and-check from the integer `b` to the dyadic `b = 2^-k`.

## Consequences

- `gamma` and `beta` return correct values (and honest flags) at tiny and
  operand-mismatched precisions, closing the last of the exact-walk
  precision-leak family. The common (moderate-precision) path is
  unchanged: the fast paths already fit at `target + 64` there, and the
  reciprocal dispatch only fires when an operand is exactly 1.
- Verified against MPFR: `differential_gamma`, `differential_beta`,
  `beta_exact_fidelity`, `property_gamma`, `property_digamma_beta`
  unchanged.

### Inversion (failure paragraphs considered)

- *"Fix pf-8scf by only lifting the counter."* Insufficient in general:
  the accumulator also rounds at tiny target, and the clean end-state is
  a single exact walk plus a target-representability check, which also
  correctly *declines* non-representable factorials (the Ziv path then
  owns the hard-to-round tie, e.g. `Γ(4)@p1 = 6`).
- *"pf-k5ll only needs the power-of-two case."* The reciprocal is the
  exact value of `B(1, b)` for EVERY `b`, so the construct-and-check is
  correct (and more direct than exp(lgamma)) for all operands; gating it
  to powers of two would leave the general unit case on the noisier path
  for no benefit.
- *"gamma(4)@p1 should return exactly 6."* 6 is not representable at p1;
  the exact dispatch must decline and let Ziv round (a tie under nearest,
  4 under TowardZero) — returning a representable value, never the exact-
  but-unrepresentable 6.

## References

- pf-8scf, pf-ihsp, pf-k5ll (epic pf-8iji); reproducers H1, H2, H4 in the
  review harness.
- ADR-0065 (beta construct-and-check the reciprocal extends), ADR-0039
  (the exact-value-defeats-Ziv dispatch family).
- `src/math/gamma.rs` (`try_gamma_pos_integer_exact`), `src/math/beta.rs`
  (`try_beta_case4_exact`, `try_beta_one_reciprocal`).
- `tests/regression_review_2026_06_10.rs`:
  `gamma_positive_integer_exact_walk_at_build_precision`,
  `beta_case4_exact_lifts_to_build_precision`,
  `beta_one_reciprocal_is_exact_for_powers_of_two`.
