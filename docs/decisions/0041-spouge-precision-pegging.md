# ADR-0041: Spouge precision-pegging — measured, rejected

- **Status**: rejected at Phase 2b sub-slice 2b.1. The bench-targeted
  hypothesis (tighten `spouge_a_for`'s safety margin from ~700 bits
  to ~50 bits, gaining 30-44% wall-clock on direct `lgamma_round`
  calls at p ≥ 1024) wins the targeted bench but breaks downstream
  verification: `differential_zeta` fails at `ζ(-1/2)` at p=1024 in
  NearestEven mode, with the result off by ~170 bits at the mantissa
  level — orders of magnitude beyond any Ziv-driver guard could
  absorb. The slice reverts per the CLAUDE.md strict revert
  stop-loss; the bench harness and this ADR are the deliverable.
- **Date**: 2026-05-27

## Context

`src/math/gamma_stirling.rs` implements Spouge's approximation for
`ln Γ(z)` (Spouge, J.L. "Computation of the Gamma, Digamma, and
Trigamma Functions" SIAM J. Numer. Anal. 31:1, 1994). The lgamma
kernel dispatches to Spouge when `working_prec > 600` per
`STIRLING_REACH_THRESHOLD` at `src/math/lgamma.rs:260`.

The pre-Phase-2b `spouge_a_for(working_prec) = working_prec / 5 + 20`
formula (slice pf-l6s5) was deliberately conservative: at
`working_prec = 1024` it returned `a = 224`, yielding
`(a − 1/2) · log_2(a) ≈ 1745`, a margin of 721 bits over the strict
Spouge truncation requirement `(a − 1/2) · log_2(a) ≥ working_prec`.
The doc on the old formula recorded the trade-off as known: "the
margin is asymptotically wasteful but cost is linear in `a`; this
trades CPU for confidence in the bit-exactness gate."

Phase 1g (ADR-0039) closed the verification-architecture gap that
motivated the conservatism. Every kernel now carries a calibrated
`error_guard` constant (`LGAMMA_ERROR_GUARD = 24`); per-release
pf-tqzz cross-check empirically guards the bound across the full
f32 grid. **The hypothesis under test:** the pre-Phase-1g
700-bit Spouge margin is no longer load-bearing in the verification
chain; it can tighten without weakening correctness.

The Phase 2 sequencing decision in
`project_perf_before_full_sweep` named "Spouge precision-pegging for
lgamma" as the first sub-slice of Phase 2b. ADR-0027 (Karatsuba
calibration), ADR-0040 (FFT decision-gate measurement), and ADR-0037
(SmallVec rejection) jointly set the discipline: measure-before-
shipping, strict revert stop-loss on a regression, ADR is the
deliverable regardless of outcome.

### What was tried

Replace `spouge_a_for` with the smallest `a ≥ 20` satisfying
`(2a − 1) · floor_log_2(a) ≥ 2 · (working_prec + 50)`, where the
50-bit guard was budgeted as:
- `LGAMMA_ERROR_GUARD = 24` (the Ziv-side calibrated slack)
- `O(log a) ≈ 16` bits of cancellation in the partial-sum loop and
  the leading-factor logarithms
- 10 bits safety margin

`floor_log_2` (computed as `31 - a.leading_zeros()`) is a
conservative under-approximation of `log_2`, keeping the formula
pure-integer per the existing `z_min_for_target` no_std pattern.

The new formula reduced the coefficient count substantially:

| working_prec | pre-2b.1 `a` | proposed `a` | reduction |
|-------------:|-------------:|-------------:|----------:|
|          600 |          140 |          109 |       22% |
|         1024 |          224 |          154 |       31% |
|         2048 |          429 |          263 |       39% |
|         4096 |          839 |          512 |       39% |

### Bench measurement

`benches/spouge_lgamma.rs` measured `BigFloat::lgamma_round` at
`target ∈ {1024, 2048, 4096}` with `z ∈ {2.5, 10, 100}` on
`aarch64-apple-darwin`. Baseline saved as
`phase2b-spouge-baseline`; the proposed-formula run diffed against
it via criterion's `--baseline` flag.

| target | z=2.5 | z=10 | z=100 |
|-------:|---------------:|--------------:|---------------:|
|   1024 | 29.4 ms (-14.5%) | 30.8 ms (-21.2%) | 24.6 ms (-16.8%) |
|   2048 |  101 ms (-44.0%) |  141 ms (-14.5%) |  149 ms (-16.4%) |
|   4096 |  877 ms (-15.5%) |  618 ms (-32.6%) |  886 ms (-54.3%) |

Every cell improved with criterion-reported `p < 0.05`. The targeted
bench would have justified the change. The decision-gate, however,
is downstream verification, not the targeted bench.

### Downstream verification failure

`differential_lgamma`, `differential_gamma`, and `differential_beta`
all passed under the proposed formula. This was a **false
positive**: those tests sweep at `precisions ∈ [53, 113, 256]`, and
the Spouge dispatch fires only when `working_prec > 600` (i.e. for
`target_precision ≥ 537`). At the test grid those tests use, the
modified code path is never executed. The verification cycle that
the slice plan called for explicitly did exercise this gap, but the
sub-slice missed naming `differential_zeta` (which sweeps
`TRANSCENDENTAL_PRECISIONS = [53, 113, 256, 1024]`) as a required
check before declaring the verification clean.

`differential_zeta` did fail under the proposed formula:

```
assertion `left == right` failed: ζ(-1/2) at p=1024, mode=NearestEven
  left:  -2.078862249773545...887888090998015812616591311245238302601219171687564907285478829e-1
  right: -2.078862249773545...887888090998009499795715903243876671640429050998872118365581306e-1
```

The two values agree to ~257 decimal digits (~853 bits) then diverge.
At `p_target = 1024` that puts the disagreement ~170 bits above the
mantissa LSB, an error of order `2^170` ULPs — far beyond any
Ziv-driver guard could absorb. The result is not "1 ULP off near a
boundary"; it is a fundamental accuracy regression in the
composition chain (zeta-FE composes lgamma via
`Γ(1-s) = exp(lgamma(1-s))`, then multiplies by `sin(πs/2)`, `2^s`,
`π^(s-1)`, `ζ(1-s)`).

After reverting `spouge_a_for` to the pre-Phase-2b formula on the
same branch tip, `differential_zeta` passes cleanly (43 minutes,
all 300 dyadic-input cells across 5 modes × 4 precisions).
**Diagnosed: the proposed-formula change is the proximate cause.**

### Why the analysis was wrong

The error is ~170 bits, far larger than any of the cancellation
budgets the 50-bit guard absorbed. Several hypotheses about why the
Spouge truncation bound `|ε| ≤ a^(1/2-a)` does not translate cleanly
to the predicted accuracy at `working_prec`:

1. **The published bound may be a relative error on `Γ(z+1)`, not
   `ln Γ(z+1)`.** Even taking the log to convert relative-to-absolute
   does not explain a 170-bit gap; some other factor is at play.
2. **NE-rounding accumulation in the partial-sum loop interacts
   with the loss-of-precision in the final `ln(S(z, a))` step.**
   The partial-sum is dominated by `√(2π)` (~2.5) with small
   alternating-sign corrections; `ln` of that has its own
   round-off behaviour.
3. **The composition chain (lgamma → exp → multiply by sin / 2^s /
   π^(s-1) / ζ(1-s)) amplifies the Spouge-side error by some factor
   not captured by the per-op `½ ULP` count.** Catastrophic
   cancellation at specific input values (like `s = -1/2`) could
   make worst-case error orders of magnitude larger than typical.
4. **The pre-Phase-2b 700-bit margin was empirically chosen to
   absorb downstream composition effects, not documented as such.**
   The "asymptotically wasteful" framing in the old doc was
   misleading: the margin was functional, not slack.

A correct precision-pegging derivation requires either (a) a
provable bound on the full lgamma → composition chain for each
composing kernel (gamma, beta, zeta-FE), or (b) an empirical
bisection: at each working_prec, find the smallest `a` such that
every composing kernel still passes its differential test.

Sub-slice 2b.1 does not have budget to do either. The new bead
`<post-v1.0 Spouge investigation, ID assigned at slice landing>`
captures the unfinished work.

## Decision

**Revert the `spouge_a_for` change.** The pre-Phase-2b formula
`(working_prec / 5).saturating_add(20).max(20)` stays in tree. The
strict revert stop-loss (CLAUDE.md, ADR-0027, ADR-0037) fires
because:

1. The change does not pass downstream verification under any
   formulation of the per-release oracle gate. `differential_zeta`
   is not a flaky test; the failure reproduces cleanly.
2. A revised formula that absorbs the missing slack would either
   need a much larger guard (effectively reverting the win) or a
   per-composing-kernel calibration (significantly larger scope
   than this sub-slice). Both deserve separate ADRs after empirical
   investigation.

The bench harness `benches/spouge_lgamma.rs` lands in tree
unchanged as durable infrastructure. Future investigation diffs
against the saved `phase2b-spouge-baseline` measurement; the bench
is reusable without modification.

## Consequences

**No `spouge_a_for` change ships.** The `lgamma` / `gamma` / `beta`
/ `zeta` perf profile is bit-identical to the pre-2b.1 state on
main (`b4cf6a3`). The pre-Phase-2b 700-bit Spouge margin stays;
v1.0 ships with the same lgamma performance Phase 1g closed at.

**The bench is the durable artifact.** The criterion baseline
`phase2b-spouge-baseline` captures the pre-change cost profile at
`target ∈ {1024, 2048, 4096}` × `z ∈ {2.5, 10, 100}` on
`aarch64-apple-darwin`. Future investigators diff against it
without re-running the baseline.

**The verification gap exposed by this slice is a separate finding
worth recording.** `differential_lgamma`, `differential_gamma`, and
`differential_beta` do not exercise the Spouge dispatch (they sweep
precisions strictly below the Stirling-vs-Spouge threshold). For
Phase 2b sub-slices touching Spouge or any high-precision path, the
verification surface must include the kernels that compose through
the affected code at `p ≥ 1024` (`differential_zeta` is the
existing high-precision composition test). A `feedback_` MEMORY
entry captures the lesson; see the post-slice memory commit for
the cross-reference.

**Post-v1.0 follow-up bead opened.** The unfinished work — an
empirical bisection to derive a correct precision-pegging formula
that accounts for downstream composition cancellation — is filed
as a P3 bead, discovered-from this slice. Not a v1.0 blocker.

**ADR-0037 / pf-cvs / ADR-0040 precedent maintained.** Three
slices in pfloat's history have now landed as measurement-deliverable
rejection ADRs: ADR-0037 (SmallVec field-swap), ADR-0040 (FFT
decision-gate), and ADR-0041 (Spouge precision-pegging). The
project's measure-before-shipping discipline is load-bearing in
the v1.0 verification posture; each rejection ADR records a
non-trivial finding the next investigator inherits.

## Related

- ADR-0027 — Karatsuba threshold calibration. Methodology template
  (asm read, strict revert stop-loss on neutral measurements,
  host-dependent comment on chosen constants).
- ADR-0037 — `SmallVec` for mantissa, rejected. The
  measure-before-shipping; ADR-as-rejection precedent.
- ADR-0039 — Phase 1g verification architecture closure. The
  per-kernel `error_guard` calibration this slice attempted to
  benefit from. The actual calibrated bounds remain at
  `DEFAULT_ERROR_GUARD = 24`; the Spouge-margin compression that
  would have made tighter per-kernel calibration matter is
  rejected here.
- ADR-0040 — Phase 2a FFT measurement. The bench-first discipline
  this slice continues. Two consecutive slices now produce
  rejection ADRs as the deliverable; the pattern is durable.
- `~/.claude/projects/-Users-parnell-Development-pfloat/memory/project_perf_before_full_sweep.md`
  — the Phase 2 sequencing decision that named Spouge precision-
  pegging as a sub-slice. The naming stands; the implementation
  hypothesis is rejected for empirical reasons.
- `src/math/gamma_stirling.rs:422` — the `spouge_a_for` function
  this slice tried (and failed) to replace; reverted to the
  pre-Phase-2b formula.
- `src/math/lgamma.rs:260` — `STIRLING_REACH_THRESHOLD = 600`,
  the dispatch boundary above which Spouge fires.
- `tests/differential_zeta.rs` — the test that caught the
  regression; sweeps `TRANSCENDENTAL_PRECISIONS = [53, 113, 256,
  1024]` and exercises Spouge at p=1024 via the FE-branch
  composition `Γ(1-s) = exp(lgamma(1-s))`.
- `tests/differential_lgamma.rs` / `differential_gamma.rs` /
  `differential_beta.rs` — the tests that passed under the
  proposed formula but did NOT exercise the modified code path
  because their precision grids stay below the Stirling-vs-Spouge
  dispatch threshold.
- Spouge, J.L. (1994), *op. cit.* — primary source for the
  approximation and its truncation bound.
- Pugh, G.R. (2004), *An Analysis of the Lanczos Gamma
  Approximation*, PhD thesis, UBC, §3 — error analysis for Spouge.
