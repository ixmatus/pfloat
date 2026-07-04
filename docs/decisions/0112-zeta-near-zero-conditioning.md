# ADR-0112: a near-zero dispatch for ζ(s), s → 0⁻ (the pf-qt7v 0·∞ form)

- **Status**: accepted
- **Date**: 2026-07-03

## Context

Found by the R3.1 (pf-hkoj, ADR-0109) adversarial verifier. `ζ(0) = -1/2`
is dispatched exactly for `±0`, but a tiny negative normal `s = -2^-k`
reaches `zeta_fe`, where the functional equation
`ζ(s) = 2·(2π)^{s-1}·sin(πs/2)·Γ(1-s)·ζ(1-s)` is a **0·∞ indeterminate
form**: `sin(πs/2) → 0` while `ζ(1-s)` approaches the pole at 1. The two
factors are computed separately. `zeta_fe` boosts its working precision by
`pole_proximity_depth(s) ≈ k` but caps the boost at
`target + 4096 + s.precision()`. For `k` past that cap, `1-s` rounds onto
exactly 1, `zeta_borwein`'s defensive belt returns an **honest NaN**
(never a certified-wrong finite — the pf-hkoj/ADR-0109 soundness holds),
and for `k` just below the cap the FE path is slow-but-correct.

The certified-wrong class was already closed; this is the honest-NaN
(and slow-path) quality gap. Reproduced: `ζ(-2^-8000)` at target 53
returned NaN under all five modes where a representable correct value is
due.

The result is representable and cheap. Near 0,
`ζ(s) = -1/2 - (1/2)ln(2π)·s + O(s²)` (the Taylor expansion at 0;
`ζ'(0) = -(1/2)ln(2π)`). For `s < 0` the linear term is positive, so
`ζ(s) = -1/2 + δ` with `δ > 0`. Verified against mpmath 1.4.1 across
`s ∈ [-2^-4, 0)`: `δ = ζ(s) + 1/2 ∈ (0.86, 0.92)·|s|`, hence
`0 < δ < |s|` on the whole range.

## Decision

Add a **near-zero dispatch** to `zeta_kernel`, after the ζ(1) pole, the
trivial-zero, and the large-s dispatches and before the Ziv driver:

```rust
if matches!(x.sign(), Sign::Negative)
    && exponent_of(x) <= -(i64::from(target_precision).saturating_add(4))
{
    let two = ci(2, target_precision);
    let (half, _) = ci(-1, target_precision).div(&two, NearestEven);
    return round_with_infinitesimal(&half, Sign::Negative, /*subtracts_magnitude=*/true,
                                    target_precision, mode);
}
```

**Soundness.** The trigger `e_s ≤ -(target+4)` gives
`|s| < 2^(e_s+1) ≤ 2^-(target+3)`, and since `e_s ≤ -5` (any target ≥ 1)
we have `|s| < 2^-4`, inside the verified range, so `0 < δ < |s| <
2^-(target+3)`. That is strictly below the half-ulp `2^-(target+2)` of
`nextUp(-1/2) = -1/2 + 2^-(target+1)`, so `ζ(s) = -1/2 + δ` lies strictly
inside `(-1/2, nextUp(-1/2))` at less than half the gap from `-1/2`. The
correctly-rounded value is therefore a pure function of the mode:

| mode | result |
|---|---|
| NearestEven, NearestAway, TowardNegative | `-1/2` |
| TowardPositive, TowardZero | `nextUp(-1/2)` |

`round_with_infinitesimal(&(-1/2), Negative, subtracts_magnitude=true, …)`
represents `-(1/2 - ε) = -1/2 + ε` with `ε` placed strictly below the
target rounding boundary and produces exactly that table (the rim-hardened
helper, ADR-0107; the same tool the large-s ζ→1 case uses). It raises
INEXACT (ζ at a nonzero dyadic is irrational).

**Completeness.** The FE path NaNs only for `k > target + 4096 +
s.precision()`; the dispatch fires for all `k ≥ target + 4`, which is far
shallower, so every NaN'ing input is covered. For `k` in the overlap
`[target+4, target+4096+prec]` the FE path was slow-but-correct and the
dispatch returns the identical mode-neighbour, now DoS-free (no big
evaluation). For `k < target + 4` (`δ` resolvable, above the ulp) the
dispatch does not fire and the FE path resolves the real value unchanged
(the `ζ(-2^-53)@53` control).

**Disjointness.** `|s| < 2^-4` is disjoint from the trivial-zero
neighbourhoods `s ≈ -2n` (pf-hkoj) and from the `s > 0` Borwein region;
`±0` and the exact trivial zeros are dispatched upstream. The dispatch
touches only tiny-negative normals.

## Consequences

- Positive: the honest-NaN region becomes correct and fast; the
  slow-but-correct overlap region becomes fast. No new kernel, no cost
  (a division and a rim-hardened round), DoS-free by construction.
- The fix is a *value* dispatch justified by a *strict* magnitude bound,
  not a transcendence claim, so it does not re-introduce the ADR-0105
  irrationality overclaim: `δ ∈ (0, 2^-(target+3))` is a rigorous interval
  from the mpmath-checked bound `δ < |s|`, independent of whether `δ` is
  rational.

### Inversion (failure paragraphs considered)

- *"Raise/coordinate the zeta_borwein cap so the FE path resolves the
  pole (the bead's literal suggestion)."* Rejected: to resolve `1-s` from
  the pole needs working ≥ k, unbounded for `s = -2^-(2^60)` — a DoS. Any
  capped variant still needs a fallback for the deep tail, and that
  fallback is exactly this dispatch, which alone is complete. Coordinating
  the cap would add cost and code for no case the dispatch does not
  already handle correctly.
- *"The neighbour table is wrong for TowardZero."* Checked explicitly:
  `ζ(s) = -1/2 + δ` is negative with magnitude `1/2 - δ < 1/2`; TowardZero
  rounds to the smaller magnitude `1/2 - 2^-(target+1)`, i.e.
  `nextUp(-1/2)`. Confirmed by the five-mode reproducer.
- *"δ could be ≥ |s| for some s in range, breaking the half-ulp bound."*
  Refuted by the mpmath sweep over `[-2^-4, 0)`: `δ/|s| ∈ (0.86, 0.92)`,
  with the worst case at the shallow end `s = -2^-4` (ratio 0.86). The
  trigger keeps `|s| ≤ 2^-5`, strictly inside.
- *"The dispatch could fire where the FE path would have given a
  different (resolvable) value."* Refuted: firing requires `δ <
  2^-(target+3) <` half-ulp, so both the dispatch and a correct FE
  evaluation round to the same mode-neighbour; the `ζ(-2^-53)@53` control
  guards the just-below-threshold boundary where the value is still
  resolvable and the dispatch must not fire.

## References

- pf-qt7v (this defect), discovered-from pf-hkoj (ADR-0109).
- ADR-0109 (zeta_fe trivial-zero conditioning), ADR-0098 (pole-proximity
  boost), ADR-0107 (`round_with_infinitesimal` rim hardening), ADR-0105
  (conditional-soundness posture the value bound respects).
- `src/math/zeta.rs` (`zeta_kernel` near-zero dispatch),
  `src/rounding.rs` (`round_with_infinitesimal`).
- `tests/regression_review_2026_06_10.rs`:
  `zeta_near_zero_negative_is_the_correctly_rounded_neighbour`.
- Oracle: mpmath 1.4.1, `ζ(s) + 1/2 ∈ (0.86, 0.92)·|s|` on `[-2^-4, 0)`.
