# ADR-0098: input-structure-aware resolution for zeta, trig reduction, and beta

- **Status**: accepted
- **Date**: 2026-06-11

## Context

The 2026-06-10 workspace review (epic pf-8iji; pf-gg96, pf-k68i, pf-pdda)
confirmed three certified-wrong-answer defects with one shape: the input
encodes proximity — to the ζ pole at 1, to a multiple of π/2, to a Γ pole of
the beta denominator — finer than the resolution the kernel works at, the
proximity collapses to exact coincidence, and the resulting special value is
certified (`half_width` of a zero or infinity is 0, so the Ziv interval test
accepts immediately).

1. **zeta** (pf-gg96): `ζ(1 + 2^-5000)` at p5001 → 53 returned **+∞ with
   Status OK**. The conditioning probe rounded `s` to target+8 bits — the
   probe itself collapsed — so the working precision got no boost; the
   working round of `s` then made `1 − 2^{1−s}` exactly 0; and the
   `DIV_BY_ZERO` from `η/0` was discarded on the `let (zeta, _) = …` line.
   The truth `≈ 2^5000 + γ` is representable (it rounds to exactly `2^5000`
   at 53 bits).
2. **trig reduction** (pf-k68i): `sin(RN(π, 2048))` → 53 returned **−0 with
   Status OK**. `reduce()`'s product width clamped to [2048, 4096];
   `y = x·(2/π)` rounded to exactly 2.0, the residual subtracted to exact
   zero. Affects sin/cos/tan/cot/sec/csc through the shared reduction.
3. **beta** (pf-pdda): `B(0.5 + 2^-60 @p61, −3.5 @p3)` returned **+0 with
   Status OK** where the truth is **negative** (−2.4913e-18; the review's
   claimed truth 0.089 was another reviewer-oracle error — the corrected
   value makes the defect worse). The ADR-0030 case-5 dispatch classified
   the pole on `a+b` rounded to max(operand precisions), where
   `−3 + 2^-60` ties-evens to exactly −3.

## Decision

Resolve input-encoded proximity exactly, or grow resolution to the input's
precision; never let a collapsed special reach the interval test.

**zeta: exact conditioning probe, input-scaled cap, loud belt.** The probe
computes `s − 1` at `max(precision(s), target + 8)` bits — Sterbenz-exact
near 1, and elsewhere only the exponent matters. The working-precision cap
gains a `+ precision(s)` term: the conditioning boost is bounded by the
input-encoded proximity, itself bounded by the input precision, so cost
stays proportional to what the caller supplied (the old fixed cap would
re-collapse deep inputs). A defensive belt refuses a zero `1 − 2^{1−s}`
factor with a NaN the driver cannot certify instead of dividing
(unreachable with the fix; the discarded-`DIV_BY_ZERO` certified-inf path
is the named trap).

**trig reduction: grow the product width on realized collapse.** The
residual `y − q` must clear the product's rounding-noise floor
`2^(e_y + 1 − mul_prec)` by `working + 8` bits; an exactly-zero residual is
always unresolved (`x·(2/π)` is irrational for dyadic `x ≠ 0`). On failure
the width doubles, up to a cap derived from the irrationality measure of π
(`μ(π) ≤ 7.103205334137`, Zeilberger–Zudilin 2020,
`docs/references/zeilberger-zudilin-pi-2020.md`): a px-bit dyadic with
exponent `e_x` cannot sit closer to `q·π/2` than `~2^-(μ·(px + e_x) + c)`,
so `8·(px + e_x) + working + 256` always resolves (`8 > μ` with slack).
Generic inputs resolve on the first pass at the old cost; `two_over_pi_at`
already computes live past its 4096-bit table. The cap fall-through —
unreachable by the bound — returns `None` (the kernels' existing
NaN + INVALID range path): a loud refusal instead of a certified −0.

**lgamma/digamma reflections: pole-proximity pre-boost.** The adversarial
verification of this slice found the family's missing fourth member,
pre-existing and newly load-bearing through beta's exact-sum handoff: the
negative-axis reflections compute `π·x` at the working width, so
input-encoded proximity to the *poles* (the integers) collapses before
sin/cos see it, and the realised-cancellation probe never fires because
the result is O(depth), not near zero (`lgamma(−3 + 2^-80 @p84) → 53` was
wrong from bit ~41; `B(0.5 + 2^-120, −3.5)` came back with the sign of Γ
flipped). The depth to the nearest integer is exactly computable on `x`'s
own grid (`pole_proximity_depth`), and both reflection closures now
pre-boost their working precision by it — the ADR-0097 asin/ln pattern.

**beta: classify the pole on the exact sum.** `a + b` is computed exactly
at its bit-span precision (`top exponent − bottom ulp position + 2` bits;
the add at that width is exact by construction, debug-asserted), within a
budget `2(pa + pb) + 2·target + 4096` that keeps the allocation
proportional to caller-supplied sizes (the bignum DoS-budget posture). Past
the budget the rounded sum is sound: a dominant negative-integer operand
was already dispatched as a Γ-operand pole earlier in the decision table,
and a dominant non-integer stays non-integer under a below-half-ulp
perturbation. The exact sum also feeds the signed-lgamma evaluation, so a
near-pole `a + b` reaches lgamma's reflection — and its ADR-0097
cancellation boost — unsnapped.

## Consequences

- The regression lane pins all three: ζ(1+2^-5000) = 2^5000 at 53 bits
  INEXACT (34s — the 5000-bit conditioning boost is real work, paid only on
  deep inputs); sin(RN2048(π)) = −9.1948…e-618 bit-exact; the beta
  reproducer negative and bit-exact at p61.
- cos/tan/cot/sec/csc inherit the reduction fix; gamma/zeta-FE
  compositions reach beta/lgamma with exact near-pole sums.
- Failure modes considered (inverted): (1) the trig resolution test could
  accept a residual that only looks resolved — it demands clearance above
  the noise floor of the *product*, the only rounding in the chain, and
  the subtraction below it is Sterbenz-exact; (2) the beta span budget
  could mis-classify past its edge — the two sound-fallback arguments
  above are each guarded by the earlier pole dispatches, and the budget
  scales with the target so the perturbation sits below anything the
  caller can observe; (3) the zeta belt could mask a real pole — `s = 1`
  is dispatched exactly upstream, so a zero factor is only reachable
  through a resolution failure, which the NaN surfaces as uncertifiable
  rather than silently wrong.
- The μ(π) and μ(ln 2) caps (this ADR and ADR-0096) share the same
  derivation shape; the references registry carries both bounds with
  factor-of-safety notes so neither decimal is load-bearing.

## Related

- Issues: pf-gg96, pf-k68i, pf-pdda (closed by this ADR), epic pf-8iji.
- Review: `~/.claude/plans/pfloat-workspace-review-2026-06-10.md` Theme 1
  items 4, 5, 7; reproducer checks F4/J1 and the corrected beta truth in
  `~/.claude/plans/pfverify-harness/`.
- References: `docs/references/zeilberger-zudilin-pi-2020.md`,
  `docs/references/marcovecchio-log2-2009.md`.
- Other ADRs: ADR-0026 (zeta algorithm), ADR-0030 (beta domain decision
  table), ADR-0050 (cancellation-to-zero class), ADR-0095/0096/0097 (this
  arc).
