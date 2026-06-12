# ADR-0104: the tiny-x family gains input-precision arms and depth hints

- **Status**: accepted
- **Date**: 2026-06-12

## Context

The R2.1 grid derivation (ADR-0102) found, and run-verified, a defect
in every shipped ADR-0059 tiny-x fast path (pf-fbjn): the triggers
tested only the target-side depth (`e ≤ −(target + 2)`), but the
`round_with_infinitesimal` boundary-gap argument also needs the series
correction to clear the *input's* grid. A high-precision input parked
next to a rounding-change point — `x = 2^-600·(1 + 2^-53 − 2^-2000)`
at p2001, sitting 2^-2600 below the p53 midpoint while atanh's +x³/3
correction (≈ 2^-1801.6) carries the truth above it — makes the fast
path round the wrong side: the residue (position `e − p − 2`) stays on
the input's side of the boundary the correction crosses. One wrong
ulp, INEXACT, in atanh/asinh/sinh/tanh (cubic corrections, position
`3e`) and expm1/log1p (quadratic, position `2e`).

R2.1 also established that an arm alone is no fix: the Ziv
fall-through rounds the input collapsed onto the same midpoint
(working capped below the input precision) and returns the identical
wrong side. The repair had to wait for ADR-0103's depth-hint driver —
a deep rung that takes the input at full precision and certifies the
true boundary side. The same machinery closes the band/hole residues
ADR-0102 recorded for its own three kernels (atan's
precision-past-the-arm band, hypot's big-operand hole, atan2's
inexact-ratio deep zone).

## Decision

**Arms.** Each family trigger gains the input-precision conjunct
derived from its correction's position: the cubic sites require
`2|e| ≥ p + 6` (correction `< 2^(3e+2)` strictly inside the
boundary-free zone `≥ 2^(e − p − 1)` needs `2|e| > p + 1`; two-plus
bits of slack), the quadratic sites `|e| ≥ p + 3` (correction
`< 2^(2e+1)` against the same zone needs `|e| > p`). Arm-failing
inputs fall through to the driver. The existing target-side arm is
unchanged, so the fast paths never fire on anything they previously
refused.

**Hints.** Every family driver call passes the lazy ADR-0103 hint
`max(correction-depth, p) + 64`: the deep rung's working precision
then covers the input exactly and resolves the tail-versus-correction
comparison that decides the boundary side. The correction-depth is
`5|e|` for the cubic sites and atan, and `3|e|` for the quadratic
sites — one series term DEEPER than the leading correction. The first
draft used the leading-term depths (`2|e|`/`|e|`) and the slice's
adversarial verification refuted it with a cheaply constructible
family: when the leading correction is exactly dyadic (atanh's `M³/3`
terminates whenever `3` divides the right mantissa structure — a
third of p53 midpoints; expm1's `x²/2` always terminates), an
adversary cancels it exactly inside a modest input precision and the
boundary side falls to the NEXT term, one rung past the hint; the
fallback returned NE/NA one ulp wrong (honest INEXACT) on inputs the
old unsound fast path happened to get right. The deeper floors cover
that whole first rung (run-verified on the verifier's three
reproducers, now matching mpmath in every mode); rung `k ≥ 2` needs
compounding divisibility coincidences (`15 | …`, `105 | …`) plus
growing precision, and stays on the documented honest-INEXACT floor —
the ladder is unbounded in the limit, so no finite hint closes it.
atan2's hint comes from the ratio exponents and the operands'
combined precision (right half-plane only); hypot's is
`max(2·gap + 2·p_small, 2·p_big) + target + 64`, justified by the
exact dyadic `truth² − B²` whose bottom bit those quantities bound.

**The expm1/log1p internal boost caps are removed.** Deriving the
arm's failure modes found that the arm alone would have *resurrected
the certified-zero class*: an arm-failing very deep input
(`p ≥ |e| − 2` with `|e|` past the closures' internal `+1024`
cancellation-boost cap) reaches the driver, the capped composition
collapses to exactly 0, and `half_width(0) = 0` certifies it with
Status OK at the first rung — the deep rung never runs, because
certification (of the collapsed zero) succeeds. The caps predate the
ADR-0059 fast paths; with the arms in place every driver-reached tiny
input has `|e| ≤ max(p + 2, target + 2)` (deeper inputs take the fast
path), so the uncapped boost is input-proportional by construction
and the caps were pure hazard. Removed, with the lane row pinning the
zone (`expm1`/`log1p` at `(1 + 2^-2997)·2^-3000`, p2998: nonzero,
correct, INEXACT, every mode).

**hypot's exactness audit.** Wiring the hole row exposed a second lie
the hint cannot reach: under the nearest modes a collapsed eval
*legitimately certifies* — the interval sits strictly inside the
half-ulp, the value is even correct — but the status claims OK while
the truth is inexact (the interval test certifies values, not
statuses). hypot is algebraic, so its exact set is decidable by
construct-and-check: on a claimed-OK result, verify `v² = x² + y²`
with exact dyadic arithmetic at an audit width covering every
genuinely-exact case's span; any nonzero residue (or any audit step
rounding — unreachable for true exacts) forces the honest INEXACT.
Runs only on the rare claimed-exact path. This is the β(1, 2^k)
construct-and-check pattern (ADR-0065) on hypot's exact set.

## Consequences

- The lane pins six rows, every one red at the pre-slice commit by
  run and mpmath-pinned (8000–12000 bits) before code: atanh and
  asinh parked at the midpoint from both sides (the bead's verified
  constructions), expm1 and log1p (the quadratic pair), atan's band
  (`2^-600 + 2^-2199` at p1600, TowardZero now pred), and hypot's
  hole (`1 + 2^-4990` at p5000 against `2^-700`: TowardPositive now
  `nextUp(1)` and — the audit's row — NE keeps the value 1 but
  finally flags INEXACT where `(1, OK)` was the falsely-exact
  defect).
- sinh and tanh carry the identical arm and hint with no dedicated
  lane row (the mechanism is shared line-for-line with atanh/asinh;
  the slice verifier probes them independently).
- Costs: the fast paths are unchanged on every input they previously
  served correctly; arm-rejected inputs pay the legacy burn plus one
  deep evaluation proportional to their own precision and depth —
  only inputs that actually park structure that deep pay it.
- Failure modes considered (inverted):
  1. **The arm could reject inputs the fast path handled correctly**
     — it does, by design (correct-but-unprovable cases inside the
     gray zone now go to the driver); the deep rung returns the same
     correct answer with certification instead of luck.
  2. **The hint could under-reach the parked depth** — it DID, in
     the first draft: the terminating-correction ladder above, found
     by the slice's verifier with run reproducers (`p ≈ 1360` at
     `|e| = 600` parking the decision at `5|e| ≈ 3003`). The deeper
     floors close rung 1; the residual ladder and the genuinely
     Diophantine parks (`δ` matching the transcendental correction
     beyond any input-proportional budget) exhaust the deep rung and
     fall back 1-ulp-honest INEXACT — never flipped-with-OK, which
     the verifier confirmed across the family. The floor is
     documented in ADR-0103; no effective bound exists.
  3. **The audit could flag a genuinely exact hypot as INEXACT** —
     impossible by construction: a true exact has `x² − v² = −y²`
     with span inside the audit width, every audit operation exact,
     residue exactly zero. The Pythagorean lane controls (3-4-5 and
     the gap-29 triple) pin this.
  4. **The rim composition** (`round_with_infinitesimal`'s residue
     saturation near `i64::MIN`, pf-a77o) is deliberately NOT
     guarded at the six family sites: the defect there is
     pre-existing and recorded on pf-a77o with reproducers, and
     R2.5 fixes it at the root for every caller — adding six interim
     guards that R2.5 immediately deletes would be churn. ADR-0102's
     three guards stay until then.

## Related

- Issues: pf-fbjn (closed by this ADR), epic pf-8iji; pf-a77o (the
  rim root, R2.5), pf-7nnw (Opus arc).
- Other ADRs: ADR-0059 (the fast paths gaining the arms), ADR-0102
  (the derivation that found the family hole and the band/hole
  residues this closes), ADR-0103 (the depth-hint machinery),
  ADR-0065 (the construct-and-check precedent).
