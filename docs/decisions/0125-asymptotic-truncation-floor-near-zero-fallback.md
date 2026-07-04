# ADR-0125: near a zero, a divergent asymptotic hands off to a convergent series (the pf-1vzg floor)

- **Status**: accepted
- **Date**: 2026-07-04

## Context

The R4.12 completeness probe (pf-b8w1) filed pf-1vzg: the large-argument
oscillatory paths `ci_asymptotic`, `airy_asymptotic_neg`,
`bessel_j_asymptotic`, and `bessel_y_asymptotic` compute `f·sin ∓ g·cos`
(and the `cos φ·P ± sin φ·Q` / `[±cosω, ±sinω, …]` analogues) directly
under `ziv_round` with no `cancellation_boosted` wrapper, unlike their
series siblings. Near a function zero the two terms cancel
catastrophically, and past the Ziv guard cap the kernel certifies a wrong
value. The bead prescribed the sibling fix: wrap each asymptotic
evaluation in `cancellation_boosted`.

Implementing that fix falsified it for three of the four kernels. A
DIVERGENT asymptotic series summed to optimal truncation carries an
IRREDUCIBLE error floor equal to its smallest retained term, `≈ e^{−c·x}`,
and no working precision lowers it (truncation error is precision
independent). Near a zero the true `|f(input)|` can fall far below that
floor; the asymptotic then returns `≈ floor`, a wrong value.
`cancellation_boosted` only grows working precision, so it converges to
the floor and certifies it, now fast, having removed the Ziv escalation
that at least presented as a hang. The decisive rule:

  **`cancellation_boosted` rescues a near-zero asymptotic only when the
  truncation floor is DEEPER than `ZIV_GUARD_CAP` (1024).**

Floors measured at target 53 with deep reproducers
(`tests/regression_review_2026_07_04_r5.rs`, oracles
`scratchpad/gen_r51_oracles_v2.py`, mpmath 1.4.1 on bit-identical
dyadics):

| kernel | floor at the test `x` | vs cap 1024 |
|---|---|---|
| Airy Ai, `x ≈ −140` | `2^-1593` (`e^{−(2/3)|x|^{3/2}}`) | deeper: boost suffices |
| Bessel J₀, `x ≈ 150` | `2^-547` | shallower: boost is useless |
| Bessel Y₀, `x ≈ 148` | `2^-547` | shallower |
| Ci, `x ≈ 100.5` | `2^-148` (`e^{−x}`) | shallower |

The defect is reachable only through the explicit low-target public API
(`*_round(low_target, deep_input)`) near a zero. A normal `.ci()`/`.j0()`
on a deep input uses target equal to the input precision, whose dispatch
picks the convergent path already; only the low-target form takes the
asymptotic there.

## Decision

Detect the near-zero, below-floor condition and hand off to a convergent
method that has no floor.

1. **Each asymptotic returns `(value, op_scale, floor_exp)`.** `op_scale`
   is the largest cancelling term's exponent (the realised cancellation
   scale, ADR-0097); `floor_exp` is the smallest retained term's exponent,
   both lifted to the result scale. `si_ci_f`/`si_ci_g` expose their
   smallest term for Ci; the Bessel kernels already track `prev_mag`.

2. **A shared driver `ziv::asymptotic_reliable`** grows the working
   precision (like `cancellation_boosted`) until the result RESOLVES, to
   its true value or to the floor it cannot cross, then applies the
   soundness test. A single low-working probe is insufficient: a deep
   near-zero input is truncated there to a shallow, large-magnitude result
   that looks resolvable (this is exactly why the J₀ deep case first
   passed the reliability test and still returned the floor). Sound to
   round via the asymptotic iff the realised cancellation `C = op_scale −
   result_exp` resolves before the working precision reaches the floor:
   the Ziv driver needs `target + C` working bits and stays sound only
   while working `≤ (result_exp − floor_exp) + guard`; combining and
   cancelling `guard` gives `2·C + target ≤ op_scale − floor_exp`. A
   result stuck at the floor has `C = op_scale − floor_exp` and fails.

3. **Reliable ⇒ `cancellation_boosted(asymptotic)`, not plain `ziv`.**
   The asymptotic evaluated at the Ziv working precision loses `C` bits to
   the cancellation against a fixed internal guard, so plain `ziv`
   resolves only `C < ~guard_cap`. At large `x` the floor is deep, the
   reliable-`C` ceiling is high, and a moderate near-zero can exceed the
   cap; `cancellation_boosted` (with the result above its floor, so it
   converges to the true value) resolves it.

4. **Unreliable ⇒ `cancellation_boosted(convergent)`.** The fallback is
   `ci_series` (Ci), `bessel_j_tiny` (the Maclaurin series, no fixed cap),
   and `bessel_y_series` (the log series). Each now returns its true
   `op_scale` equal to its largest partial term. For large `x` the terms
   peak at `≈ 2^{x·log₂ e}`, so the prior hardcoded `op_scale = 4` on
   `ci_series` undercharged the cancellation; correcting it folds the
   R5.2 series-boost concern (pf-6naq) into the same change for these
   three kernels. The convergent series has no truncation floor, so the
   boost drives it to any near-zero depth (cost input proportional, the
   DoS-budget posture; the deep rows are release gated).

5. **Airy keeps the bare `cancellation_boosted` wrapper.** Its floor
   (`≈ 2^-1593` at `|x| ≥ 128`, the asymptotic-regime minimum) is deeper
   than the cap, so the wrapper correctly resolves every reachable
   near-zero in the gap `(1024, floor)`; the R5.1 Airy deep reproducer
   (`D = 1294`) passes.

## Consequences

- The certified-wrong-near-a-zero family is closed for Ci, Bessel J, and
  Bessel Y at ANY depth: their convergent fallback has no floor.
- **Named boundary (competence limit).** Airy retains its asymptotic
  floor; an input DEEPER than `≈ 2^-1593` (a `> 1593`-bit input within
  `2^-1593` of an Airy zero at `|x| ≥ 128`) still certifies the floor
  value. This is the same class at an extreme depth; closing it needs the
  same pattern (`airy_uv_sums_split` floor, `airy_series` op_scale,
  reliability plus fallback) and an expensive `~3200`-bit-working
  reproducer. Tracked as a follow-up (discovered-from pf-1vzg). Accepting
  it here is a named decision, not a gap.
- The first-draft fix (the bare wrapper for all four) was a
  MISDIAGNOSIS. It read as a fix because it removed the hang symptom while
  leaving the wrong value; for Ci and Bessel it made the certified-wrong
  faster and quieter. The "consumed kernel lying at its own cap" shape:
  the asymptotic returns its truncation floor as if it were the answer.
- **Inversion: the reproducer must be DEEP, past the floor, not merely
  past the cap.** A shallow near-zero within the floor is slow-but-correct
  (resolvable), not wrong; a first reproducer picked at `D ≈ 245`
  (below Airy/Bessel floors, above Ci's) mislabels the fix. The deep rows
  parse the zero `> floor` bits deep and certify bit-exact against mpmath.
- **Inversion: reliability at a single working precision is a trap.** The
  probe must resolve the input first (grow working to the true value or
  the floor); a low-working probe of a deep input reports a shallow,
  reliable-looking result.
- Cost: the reliability decision and the convergent fallback each grow
  working toward the floor depth near a zero, seconds at the deep rows.
  Away from zeros the decision resolves at the base working precision in
  one or two evaluations, so the common path is unchanged.

## References

- pf-1vzg (the P1 defect); pf-b8w1 (the R4.12 probe that filed it).
- ADR-0097 (charging the realised operand scale); ADR-0110 (the geometric
  `cancellation_boosted`); the 2026-05-29 review root cause 2 (the
  cancellation-near-a-zero family this closes for the asymptotic paths).
- pf-6naq (the fixed-cap series boost), folded here for Ci/Bessel.
- DLMF 6.12 (Si/Ci auxiliaries), 9.7 (Airy asymptotics), 10.17 (Bessel
  Hankel asymptotics), 6.6 / 10.2 / 10.8 (the convergent series).
