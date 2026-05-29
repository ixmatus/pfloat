# ADR-0050: tanh stable `expm1` form replaces the cancelling composition and its tiny-x short circuit

- **Status**: accepted
- **Date**: 2026-05-29

## Context

`tanh` evaluated `tanh(|x|) = (1 − e^{−2|x|}) / (1 + e^{−2|x|})` at the
Ziv driver's working precision, plus a tiny-`|x|` short circuit that
returned `|x|` directly (slice p1.4, closing pf-7d7). The short circuit
existed because the bare numerator `1 − e^{−2|x|}` collapses for tiny
`|x|`: `e^{−2|x|}` rounds to exactly `1`, the numerator becomes exactly
`0`, and the Ziv interval test certifies the false `0` (since
`half_width(0) = 0`). Returning `|x|` sidestepped that for round-to-
nearest, where `tanh(|x|) ≈ |x|` to within half a ULP.

The pf-hcz4 cross-check local pre-flight (the directed-mode-unswept
MPFR shards, run before the EC2 campaign) flagged `tanh` with a
violation on essentially every directed-mode subnormal input
(0 violations in NE/NA, 65535 each in TZ/TP/TN over a 65536×5 sweep).
Investigation (recorded in pf-zhcy):

- `tanh` is **correctly rounded** in every mode. A direct probe of the
  smallest subnormals confirmed TZ/TN round down (`0` for `2^−149`,
  `x − 1 ULP` above it) and TP/NE round to `x`, matching the correct
  directed rounding of `tanh(x) < x`.
- The cross-check assertion
  `|eval_w − midpoint| ≤ 2^(error_guard − working)·|midpoint|` is
  exactly the soundness condition for the Ziv interval test's
  half-width. Its failure means `error_guard = 24` did **not establish**
  interval-test soundness for tiny-`x` `tanh`. Cause: the short circuit
  returns the grid point `|x|`, so under directed rounding the symmetric
  Ziv interval around `|x|` always straddles a rounding boundary; the
  driver never converges there and climbs past the short-circuit-disable
  threshold (`working ≈ 320`) into the bare composition, where
  `1 − e^{−2|x|}` cancels and loses ~148 bits. At `working = 536` the
  converged intermediate held only ~388 bits, far short of the ~512 the
  bound assumes. Correctness survived only by the large boundary margin
  (`x³/3`), not by the calibration.

So the kernel was correct but its calibrated error guard was not a true
bound on its internal error for tiny `x`.

## Decision

Evaluate the numerator through `expm1`:

```text
tanh(|x|) = −expm1(−2|x|) / (2 + expm1(−2|x|))
```

This is algebraically identical to the previous identity
(`expm1(−2|x|) = e^{−2|x|} − 1`, so `−expm1(−2|x|) = 1 − e^{−2|x|}` and
`2 + expm1(−2|x|) = 1 + e^{−2|x|}`), but `expm1` computes the small
difference without the catastrophic cancellation. The numerator is
accurate to working precision for all `|x|`; the denominator lies in
`(1, 2]` and never cancels.

The tiny-`|x|` short circuit is **removed**. Its sole purpose was to
avoid the `1 − 1 = 0` collapse, which `expm1` obviates: `expm1(−2|x|)`
preserves the tiny value `≈ −2|x|` rather than rounding it to `0`, so
the numerator is `≈ 2|x| ≠ 0`. With the short circuit gone the Ziv
driver works on the accurate composition at every working precision,
the error guard becomes a true bound, and the interval test is sound.

`expm1` is an existing Ziv-driven kernel (gated behind `exp-log`, the
same feature as `tanh`); the change reuses it rather than adding a
hand-rolled Taylor path.

## Consequences

- **error_guard = 24 now holds for tanh tiny-x.** The cross-check
  assertion passes legitimately rather than by silencing; the interval
  test's half-width is a genuine error bound.
- **pf-7d7 stays closed** by construction (no cancelling numerator that
  can collapse to zero), now without a special case.
- **Simpler kernel.** The short-circuit threshold derivation and the
  `Class::Normal` exponent branch are gone; `tanh_at_w` is the stable
  three-operation composition plus the sign restore.
- **No correctness change**, verified (see below). The set of returned
  f32 values is unchanged; what changes is that the working-precision
  intermediate is now accurate, so directed-mode tiny-`x` no longer
  forces the Ziv driver into the cancelling path.
- **Generalises as a pattern.** Any small-argument composition of the
  form `1 − e^{−t}` (or `1 ± (1 ± ε)`) should route through `expm1` /
  the corresponding `…m1` primitive rather than the bare difference;
  the cross-check will flag the others the same way if they exist.

## Verification

- `cargo test --features exp-log --lib tanh`: 25/25 pass, including
  `tanh_tiny_input_round_to_nearest_returns_input`.
- pf-tqzz cross-check (`pf_tqzz_sweep --fn-id tanh --modes all`):
  **0 violations** (was 12285 on the 4096-subnormal sample), all
  20710 cells pass.
- Correct-rounding vs MPFR (`verify_input`, all 5 modes, 19773 inputs
  spanning dense tiny subnormals, a stride across the full subnormal
  range, normals, and negatives): **0 mismatches** over 96473
  certifiable checks (2392 oracle-inconclusive at MAX_PREC, kernel-
  independent).

## References

- pf-zhcy — the cross-check finding and full diagnosis.
- ADR-0049 — pf-hcz4 cross-check sweep; its no-auto-widen triage
  protocol scoped this as a separate kernel slice.
- ADR-0022 — the Ziv driver and interval test this kernel runs under.
- pf-7d7 — the original tiny-`x` collapse the removed short circuit
  addressed.
- `feedback_cross_check_trace_fragility` (active memory) — the local
  pre-flight that surfaced this.
