# ADR-0095: agm convergence floor is relative to the iterate magnitude

- **Status**: accepted
- **Date**: 2026-06-11

## Context

The 2026-06-10 workspace review (epic pf-8iji, finding pf-ddfl) confirmed by
reproducer that `BigFloat::agm(2^-300, 3·2^-302)` at p53 returned exactly the
first arithmetic mean `(a+b)/2` with Status OK: a 0.5% relative error,
certified. mpmath 1.4.1 at 4000 bits gives
`4.273399828000648542805471530695713670719e-91`; the kernel returned
`4.2954567821355107e-91`, bitwise equal to the AM.

Root cause: the Gauss-iteration convergence test in `agm_kernel` compared the
gap exponent of `|a_n − b_n|` against the absolute floor `−(w + 4)` (with `w`
the working precision). The gap of small operands starts below that floor, so
the loop "converged" before its first iteration. The same absolute floor was
unreachable for large operands (a gap of magnitude-`2^300` iterates never
falls below `2^-(w+4)` at working precision), so those inputs silently ran all
64 iterations and exited on the iteration cap; the returned average was
correct, the budget spent was not. Only operands of magnitude near 1, where
absolute and relative coincide, behaved as designed — and that is the regime
the Ziv error guard (`AGM_ERROR_GUARD`) was calibrated in.

The Ziv driver cannot catch this class: its first-iteration error model is
relative half-width `2^-(w-guard)`, and the returned AM varies smoothly with
`w`, so the interval test certifies it (the review's Theme 1 mechanism).

## Decision

Make the convergence criterion relative to the iterate magnitude, and derive
the iteration budget instead of assuming it.

**Criterion.** Converged when the gap is exactly zero, or when the gap falls
below 4 ulps of `a_n`: `exponent(|a_n − b_n|) < exponent(a_n) − w + 3`,
computed with `saturating_sub` (saturation errs toward more iterations, never
premature convergence). Derivation, using pfloat's convention
`exponent(v) = floor(log2 |v|)` so `2^e ≤ |v| < 2^(e+1)`, and
`ulp(a_n) = 2^(ae − w + 1)`: the Gauss iterates satisfy `b_n ≤ AGM ≤ a_n`,
and the kernel returns the midpoint `m = (a_n + b_n)/2`, so
`|m − AGM| ≤ (a_n − b_n)/2 < 2 ulp = 2^(ae − w + 2)`, a relative error under
`2^-(w-3)` — more than `2^20` inside the `AGM_ERROR_GUARD = 24` half-width
`|y|·2^-(w-24)` the Ziv driver charges. The canonicalization swap
(`a_n ≥ b_n`) becomes load-bearing: the test reads `a_n`'s exponent as the
magnitude, and AM ≥ GM preserves the ordering across iterations.

The threshold deliberately sits a few ulps ABOVE the working grid. Two
distinct w-bit values can never differ by less than `2^(ae − w)` relative
(one ulp of the lower binade), so any relative floor at or below one ulp —
including this ADR's first draft at `2^-(w+4)`, and the old absolute floor
read at magnitude ~1 — is unsatisfiable: the loop then terminates only by
bit-equality or the iteration cap, and iterates that settle into a
persistent 1-ulp oscillation (an AM half-ulp tie that rounds away from the
GM, roughly half of inputs) burn the whole cap for a measured ~10× cost.
The independent adversarial verification of this slice caught that dead
floor by replicating the loop and watching which exit fired.

**Iteration budget.** `max_iter` rises from 64 to 104, now a backstop rather
than the expected exit, derived in two regimes. Far regime
(`R = a_n/b_n > 2`): one step maps `R` to at most `sqrt(R)`
(`a'/b' = (a+b)/(2·sqrt(ab)) ≤ sqrt(a/b)`), halving `log2 R`; exponents span
i64, so `log2 R < 2^65` needs at most 65 halvings. Near regime (`R ≤ 2`):
the relative gap squares per step, so reaching the few-ulp floor for any
supported `w < 2^32` needs at most `log2(2^32) + 1 = 33` steps. With 6
slack: 104. A cap exit (now reachable only by exponent-saturation
degeneracies) still averages a pair within a few ulps, inside the guard.

## Consequences

- `agm` returns the correctly rounded AGM at every operand scale whose loop
  arithmetic stays inside the i64 exponent range (see the exclusion below);
  the regression lane `tests/regression_review_2026_06_10.rs` pins the
  defect input against the mpmath reference, a large-operand control, and
  precision-refinement self-consistency. The adversarial verification
  additionally confirmed bit-exact agreement with mpmath on directed modes,
  `agm(2^1000, 2^-1000)`, `agm(2^-300, 2^300)`, `agm(2^(2^62), 2^-(2^62))`,
  and p ∈ {1, 2, 2000}.
- Both former cap-burn classes go away: small operands no longer converge
  prematurely, and inputs whose iterates oscillate within an ulp (formerly
  104 iterations × up to 5 Ziv retries, measured ~10×) now exit on the
  few-ulp floor in `O(log w)`.
- **Known exclusion (pre-existing, out of scope here):** operand exponents
  within ~2^62 of the i64 rim make the loop's `mul`/`sqrt` saturate their
  result exponents; those statuses are discarded inside the eval closure, so
  the corrupted iterate is certified — worst case
  `agm(2^(2^62+1000), 3·2^(2^62+1000))` returns a value ~10^31 wrong with
  Status OK. Verified pre-existing at the pre-fix commit (the fix neither
  caused nor fixes it; it did extend correctness down to |exponent| ≈ 2^62
  on the small side). Filed as a discovered-from bead into the
  pf-a77o/pf-kh3z saturation arc; a candidate fix is exploiting AGM's
  degree-1 homogeneity to normalize operands toward exponent 0 and scale
  back through `scale_by_pow2`'s honest saturation contract.
- Failure modes considered (inverted view): (1) a wrong inequality direction
  or off-by-one in the criterion would either reintroduce premature
  convergence or never converge; the first draft of this very ADR had the
  second failure (a sub-ulp floor), caught by the adversarial loop replica,
  and the few-ulp bound is now derived against the grid-quantization
  argument above. (2) `max_iter` too small for extreme exponent spreads
  would fall back to averaging un-converged iterates; the far-regime halving
  argument bounds the spread by the i64 exponent domain itself.
  (3) Exponent arithmetic near i64::MIN could wrap; `saturating_sub`
  degrades toward extra iterations only.
- The criterion still reads `Class::Normal` exponents directly; if either
  iterate degenerates, the match arm falls through to "not converged" and
  the cap exit applies — sound, not silent.

## Related

- Issues: pf-ddfl (closed by this ADR), epic pf-8iji; pf-a77o/pf-kh3z
  (adjacent, separate arc).
- Review: `~/.claude/plans/pfloat-workspace-review-2026-06-10.md` Theme 1
  item 6; reproducer check I1 in `~/.claude/plans/pfverify-harness/`.
- Other ADRs: ADR-0015 (Gauss iteration choice), ADR-0038 (agm under the
  shared Ziv driver), ADR-0039 (equal-operand exact dispatch).
