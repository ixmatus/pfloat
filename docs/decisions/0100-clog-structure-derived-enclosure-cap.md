# ADR-0100: clog's enclosure guard derives its cap from the input structure

- **Status**: accepted
- **Date**: 2026-06-11

## Context

The 2026-06-10 workspace review (epic pf-8iji, pf-qm8a) confirmed that
`clog(1 + 2^-545 i)` at p64 returns `re = +0` with INEXACT where the truth
`≈ 2^-1091·…` (precisely `3.7694…e-329`) is representable. The real part
`ln(hypot(re, im))` is enclosed by a directed pair whose guard schedule
(`GUARDS`, topping out at `p + 1024` working bits) cannot resolve a
`|z| = 1` straddle deeper than ~512 input bits — for `z = 1 + εi` the real
part is `~ε²/2` — and the exhausted loop silently returned the unconverged
bracket end. ADR-0091 recorded a measure-zero caveat for the enclosure
fall-through; this band is constructible and positive-measure, exceeding
that caveat.

## Decision

`ln_hypot_enclosure` keeps the GUARDS schedule for everything it already
handled, then continues doubling the guard up to a cap derived from the
input structure: `a² + b² − 1`, if nonzero, is a dyadic value no smaller
than `2^bot` with `bot = min(2(e_a − p_a + 1), 2(e_b − p_b + 1))` (the
squared components' lowest possible set bit), so a guard of `−bot + 64`
always either resolves the bracket or reaches the exactly-Pythagorean case
`a² + b² = 1`, where the directed hypot pair collapses onto exactly 1 and
`ln(1) = 0` converges exactly. The residual fall-through past the derived
cap — unreachable for finite-precision inputs by the bound — reports
INEXACT through `resolve_bracket`'s non-converged arm (never OK), which is
the honest status available inside the frozen 1.0 Status surface
(ADR-0093); the ADR-0091 caveat text now describes only that genuinely
unreachable residue.

This mirrors the arc's scalar-side pattern (ADR-0097/0098): the collapse
depth is input-encoded and exactly computable, so resolution grows to meet
it instead of stopping at an output-sized schedule.

**The false-convergence trigger.** The slice's adversarial verification
refuted the first draft at depth ≥ 576: the inner scalar `hypot_round`
exhausts ITS OWN Ziv cap (the eval's NE rounding absorbs `1 + 2^-2d` to
exactly 1 whenever `2d ≥ working`) and returns a falsely-exact `(1, OK)`
(filed as pf-71u2), so the outer bracket "converged exactly" on `[0, 0]`
at the FIRST guard and the depth-scaled growth never ran — with the re
component's isolated status degraded to OK, worse than the original
defect. The repair rests on a number-theoretic fact the verifier
supplied: no nontrivial dyadic point lies on the unit circle
(`x² + y² = 2^2k` has only trivial Gaussian-integer solutions), so a
hypot bracket whose BOTH ends are exactly 1 with both components nonzero
is ALWAYS a lie. The loop treats that shape as unresolved and keeps
growing; once the widened target lets the scalar hypot resolve genuinely,
the bracket becomes honest. The first draft also fed a negative value
through `u32::try_from(...).unwrap_or(u32::MAX)` when `bot > 64` (a
latent unbounded-cost arm), now clamped. A cleaner future shape the
verifier sketched: form `δ = a² + b² − 1` exactly (it is a dyadic) and
drive the real part through `½·log1p(δ)` — no straddle at all — but the
extreme-exponent-gap case needs a sparse-sum representation pfloat does
not have, so the trigger-plus-growth stays.

## Consequences

- `pfloat-complex/tests/regression_review_2026_06_10.rs` pins the
  reproducer bit-exactly against mpmath (both components) plus a shallow
  control. The deep case costs what the depth costs (~4s at depth 545 in
  the debug lane) and only when the input actually straddles `|z| = 1`
  that deeply; the schedule is unchanged otherwise.
- Other `GUARDS` consumers (the divide and csqrt enclosures) keep the
  static schedule: their brackets converge at output-scaled precision
  (no input-encoded straddle of a log singularity); if a sibling band is
  found it should get the same derivation, not a copied constant.
- Failure modes considered (inverted): (1) a wrong `bot` bound (forgetting
  the squares double the ulp depth) would under-cap and reintroduce the
  silent fall-through — the lane's depth-545 reproducer needs the doubled
  term; (2) the false-convergence trigger could mask a genuine
  exactly-on-the-circle input looping to the cap — no such input exists:
  nontrivial dyadic points on the unit circle are impossible (the
  Gaussian-integer argument above), and the trivial ones (±1, 0i),
  (0, ±1i) carry a zero component, which the trigger's nonzero guard
  exempts (and the kernel's special-value dispatch owns them anyway);
  the cap fall-through still returns INEXACT, never OK, as the total
  backstop; (3) cost on
  adversarial inputs is capped by the input's encoded depth, which the
  EXPONENT carries (an i64, nearly free to write): a depth-2d straddle
  costs O(d)-bit arithmetic regardless of the input's precision. That is
  the honest price of resolving the band (the lane's depth-2000 case runs
  ~37s in the debug lane); a cheaper escape needs the pf-71u2 scalar fix
  or the sparse log1p design above.

## Related

- Issues: pf-qm8a (closed by this ADR), epic pf-8iji; pf-71u2 (the
  falsely-exact hypot, filed from this slice), pf-e2ow (atan2 directed
  tiny-ratio, filed from this slice).
- Review: `~/.claude/plans/pfloat-workspace-review-2026-06-10.md` Theme 1
  item 8; harness check X1.
- Other ADRs: ADR-0091 (the enclosure design and its caveat, now
  narrowed), ADR-0092/0093 (verification posture and the frozen surface),
  ADR-0097/0098 (the scalar-side siblings of this pattern).
