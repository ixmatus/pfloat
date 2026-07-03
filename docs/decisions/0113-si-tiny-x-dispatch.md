# ADR-0113: tiny-x dispatch for Si (the pf-31ql directed-mode 1-ulp fix)

- **Status**: accepted
- **Date**: 2026-07-03

## Context

`Si(x) = ∫₀ˣ (sin t)/t dt = x − x³/18 + x⁵/600 − …`, so for `x > 0` the
value is strictly below `x` (the leading correction `−x³/18` opposes
`x`'s sign; Si is odd, so `|Si(x)| < |x|` on both sides). For tiny `x`
the `si_series` evaluation collapses onto `x` once the `−x³/18` term
falls below the working precision, and `Si` relied on the Ziv driver to
recover the directed-mode side. That works while the correction is
within reach of the guard cap, but the correction sits `~2·|e_x|` bits
below `x`, and past the Ziv guard cap (`target + 1024`) the driver can
never resolve it: it collapses to `x` at every rung and rounds `x` as if
exact, so **every mode returned `x`** — 1 ulp wrong under TowardZero and
TowardNegative, where the correct result is `pred(x)` (`Si(x) < x`).
`Si(2^-1000)` at target 53 (correction ~2004 bits down, past the
`53 + 1024` cap) returned `x` under TN/TZ (pf-31ql).

`Si` is the direct analogue of the ADR-0059/0104 tiny-x family
(`asinh(x) = x − x³/6`, same magnitude-shrinking shape); it simply never
received the dispatch.

## Decision

Add the ADR-0059/0104 tiny-x dispatch to `si_kernel`, after the
special-case peel-off and before the regime/Ziv split, identical in form
to `asinh`:

```rust
let e = /* exponent_of(x) */;
if e <= -(i64::from(target_precision) + 2)
    && e.saturating_mul(-2) >= i64::from(x.precision).saturating_add(6)
{
    return round_with_infinitesimal(x, x.sign(), /*subtracts_magnitude=*/true,
                                    target_precision, mode);
}
```

The first condition (`e ≤ -(target+2)`) puts the `−x³/18` correction far
below the target rounding boundary; the second (ADR-0104, `-2e ≥ p+6`)
clears the **input's** grid so a high-precision `x` parked next to a
rounding-change point is routed to the driver's deep rung instead. Both
conditions together give `2·|e_x| > max(p, target) + 3`, so the true
correction lies strictly below the residue position
`round_with_infinitesimal` uses, making its magnitude-shrinking round the
correctly-rounded `Si(x)`: NE/NA/TowardPositive → `x`, TowardZero/
TowardNegative → `pred(x)` (or `succ(x)` on the negative side by odd
symmetry). `round_with_infinitesimal` raises INEXACT (Si at a nonzero
algebraic argument is transcendental).

`1/18 < 1/6`, so `Si`'s correction is even smaller than `asinh`'s at the
same `e`; the shared thresholds are conservative for Si.

## Consequences

- Positive: `Si` is correctly rounded under every mode at every tiny `x`,
  including the deep region past the Ziv cap where the driver could not
  recover; the fast path also avoids the driver's guard-doubling
  escalation for moderately tiny `x` (a latency win).
- The dispatch fires for all `x` with `e ≤ -(target+2)` and a clear input
  grid — wider than the buggy deep region — but its answer is the
  correctly-rounded value there too (verified: `Si(2^-60)@53` was already
  correct via the driver and stays correct via the dispatch), so widening
  is sound.
- No change for `x` above the threshold: the series/asymptotic Ziv path
  is untouched (the `differential_si` lane stays on it).

### Inversion (failure paragraphs considered)

- *"The bug reproduces at x = 2^-60."* Refuted by run: the Ziv driver
  escalates its guard and resolves the correction at ~2^-60 (correction
  only ~184 bits down, inside the `53 + 1024` cap), returning the correct
  directed value pre-fix. The defect needs the correction past the guard
  cap — `x = 2^-1000` at target 53 — which the reproducer uses.
- *"round_with_infinitesimal could mis-round a high-precision tiny x."*
  Guarded by the ADR-0104 second condition, which routes such inputs to
  the driver; without it a high-precision `x` next to a boundary would
  cross a rounding-change point the residue never reaches.

## References

- pf-31ql (this defect), epic pf-8iji; deferred from R3.3 (ADR-0110
  disposition) as the mechanical tail.
- ADR-0059 (the tiny-x round_with_infinitesimal family: atanh/asinh/
  sinh/tanh), ADR-0104 (`-2e ≥ p+6` input-grid clearance), ADR-0107
  (`round_with_infinitesimal` rim hardening).
- `src/math/si.rs` (`si_kernel`),
  `tests/regression_review_2026_06_10.rs`:
  `si_tiny_x_directed_modes_shrink_toward_zero`.
