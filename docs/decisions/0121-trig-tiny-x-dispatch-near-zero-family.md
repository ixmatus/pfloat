# ADR-0121: tiny-x dispatch for the near-zero `x ± c·x³` trig/inverse family (sin, tan, asin)

- **Status**: accepted
- **Date**: 2026-07-03

## Context

pf-7nnw: for deeply tiny `x`, the trig and inverse-trig kernels reduce
to `x ± c·x³ + …` near 0, and past the Ziv guard cap (`target + 1024`)
the reduced Taylor collapses onto `x` — the cubic correction sits
`~3·|e_x|` bits below `x`, unreachable — so the driver rounds `x` as if
exact and directed modes returned `x` itself: `sin(2^-600)` under
TowardZero returned `x` where `pred(x)` is due (and the same through the
`pfloat-libm` shell, `f64::sin_round(2^-540, TowardZero)`). This is the
Si defect (ADR-0113) one kernel family over; `atan` already carried the
ADR-0059 dispatch, the rest did not.

The claimed family in the review was
`sin/cos/tan/atan/asin/exp/cosh/sec/csc/cot`. It splits by Taylor shape:

- **Near-0, result ≈ `x`** (`x ± c·x³`): `sin` (`−x³/6`, shrinks), `tan`
  (`+x³/3`, grows), `asin` (`+x³/6`, grows), `atan` (already fixed). This
  ADR fixes `sin`, `tan`, `asin`.
- **Near-1, result ≈ `1`** (`1 ± c·x²`, or `1 + x` for `exp`): `cos`,
  `exp`, `cosh`, `sec`. A DISTINCT shape — base `1` not `x`, correction
  `x²`/linear not `x³`, so the `e ≤ −(target+2)` threshold does not apply
  and needs a separate derivation. Filed as pf-767j.
- **`≈ 1/x`** (huge, not tiny): `csc`, `cot`. Probably a different or no
  defect (the result is large, not `≈ x`); probe under pf-767j.

## Decision

Add the ADR-0059/0104 tiny-x dispatch to `sin`, `tan`, `asin`, after the
special-case peel-off and before the reduction/Ziv driver, identical in
form to Si (ADR-0113):

```rust
if e <= -(target_precision as i64 + 2)
    && e.saturating_mul(-2) >= x.precision as i64 + 6
{
    return round_with_infinitesimal(x, x.sign(), SUBTRACTS, target, mode);
}
```

with `SUBTRACTS = true` for `sin` (shrinks), `false` for `tan`/`asin`
(grow). The two-part depth clears both the target ulp and the input grid
(ADR-0104), so the round is correctly rounded (shrink: TZ/TN → `pred(x)`,
else `x`; grow: TP → `succ(x)`, else `x`) and INEXACT. `|x| ≤
2^-(target+2)` is far below any `π`/`π/2` multiple, so `sin`/`tan` need
no argument reduction and `asin` is trivially in-domain; arm-failing
inputs fall through to the Ziv driver unchanged.

## Consequences

- `sin`, `tan`, `asin` are correctly rounded under every mode at deeply
  tiny `x`, including past the Ziv cap; the fast path also skips the
  reduction + Ziv escalation for tiny inputs. The reduction/Taylor path
  is untouched above the threshold (the differential lanes stay on it).
- The near-1 family (`cos`, `exp`, `cosh`, `sec`) and the reciprocals
  (`csc`, `cot`) are a **filed follow-up (pf-767j)**, deliberately scoped
  out: their tiny-x behaviour is a different Taylor shape that a wrong
  threshold would MIS-round (worse than the current 1-ulp), so it gets
  its own derivation and reproducers rather than an unverified extension
  here.

### Inversion (failure paragraphs considered)

- *"Apply one dispatch to the whole named family."* Refuted: `cos`/`exp`/
  `cosh` are near-1 with an `x²`/linear correction, so the `x³` threshold
  `e ≤ −(target+2)` is not the right depth and would round the wrong
  neighbour; `csc`/`cot` are `≈ 1/x`, not `≈ x`, and may have no defect of
  this shape. Fixing only the verified same-shape subset is the safe move.
- *"The bug reproduces at moderately tiny x."* Like Si, the defect needs
  the cubic past the Ziv guard cap; the reproducer uses `x = 2^-1000` at
  target 53 (cubic ~3000 bits down, past `53+1024`).

## References

- pf-7nnw (this defect, near-0 family), pf-767j (the near-1/reciprocal
  follow-up filed here), epic pf-8iji.
- ADR-0059 (the tiny-x `round_with_infinitesimal` family), ADR-0104
  (input-grid clearance), ADR-0113 (the Si sibling).
- `src/math/sin.rs`, `src/math/tan.rs`, `src/math/asin.rs`;
  `tests/regression_review_2026_06_10.rs::trig_tiny_x_directed_modes_shrink_and_grow`.
