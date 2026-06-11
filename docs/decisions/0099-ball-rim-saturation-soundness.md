# ADR-0099: Ball widening at the exponent rim and the parse-budget collapse

- **Status**: accepted
- **Date**: 2026-06-11

## Context

The 2026-06-10 workspace review (epic pf-8iji; pf-m37w, pf-1bqy) confirmed
two Law-1 (containment) unsoundness classes in pfloat-ball:

1. **Bracket saturation** (pf-m37w): pfloat's scalar mul/div saturate the
   result exponent at the i64 rim (the no-emax contract) and flag
   OVERFLOW/UNDERFLOW — but `Ball::mul`/`div` discarded the directed
   brackets' statuses. All three directed products clamp to the same
   finite value, the bracket spread vanishes, the propagated radius of a
   point ball is zero, and the resulting "exact" ball excludes the truth.
   Correcting the record on the reproducer: the review's BL3 input
   (`k = 4611686018427387900`) never crosses the rim — its square has
   exponent `i64::MAX − 7`, is exactly representable, and the exact ball
   it produces is CORRECT (the harness comment claiming `2k > i64::MAX`
   was wrong by 7; the same misadjudication family as the review's asin
   and beta oracles). The class is real from `k = 2^62` upward, where the
   scalar ops genuinely saturate and flag.
2. **Parse-budget collapse** (pf-1bqy): `Ball::parse_decimal` already
   routes through a directed-pair interval, but pfloat's parser saturates
   past-budget tiny magnitudes to ±0 in every rounding mode (the pf-mw6u
   defect, separate arc) while flagging UNDERFLOW; both interval ends came
   back +0 and the ball collapsed to an exact [0, 0] excluding the
   positive truth (`0.1e-1000000`).

## Decision

**Saturation-driven widening in mul/div.** The directed-bracket statuses
are now captured and merged. On OVERFLOW saturation the radius goes
`Mag::INFINITY` (the truth may exceed every representable; an unbounded
end is the only sound enclosure), and OVERFLOW is surfaced in the returned
Status. On UNDERFLOW saturation the exponent clamps UP toward the floor,
overstating the true magnitude, so the radius is widened by the clamped
midpoint's own magnitude — the ball then reaches zero on the far side of
the truth — and UNDERFLOW is surfaced. Ordinary products and quotients are
untouched (the widening keys off flags the scalar ops only raise at the
rim).

**Containing interval for the flagged parse collapse.** When either
directed parse flags UNDERFLOW and the corresponding interval end
collapsed to zero, that end is replaced by a cheap power-of-two bound on
the literal's magnitude: `|value| < 10^(E + D + 1) ≤ 2^bound`, with the
`log2(10)` rational chosen per sign so integer arithmetic never rounds the
bound down. Loose but sound; tightness returns when pf-mw6u lands a
mode-aware parse in pfloat itself (this fix is deliberately ball-local so
the frozen pfloat 1.x parse surface is untouched in this arc).

**The radius pipeline gets the same treatment.** The slice's adversarial
verification refuted the first draft one level down: the propagated-radius
terms (`up()`'s upper bounds in add/sub/mul/div, and div's
denominator/quotient pair) still discarded scalar statuses, and two
run-verified Law-1 violations followed with Status OK — a mul whose
radius term `|a|·rb` overflow-clamped DOWN (under-sizing the radius while
the midpoint brackets were exact), and a div whose denominator LOWER
bound underflow-clamped UP (over-stating the denominator, under-sizing
the propagated quotient, excluding a representable member of the quotient
set). The direction analysis: every radius term is an upper bound, so
OVERFLOW (clamp-down) is the unsound direction and triggers
`Mag::INFINITY`, while UNDERFLOW (clamp-up) over-estimates and is sound;
the denominator is the one lower bound, so its directions invert. The
returned Status is NOT augmented for radius-only blowups (the midpoint
computation was clean; the infinite radius — `is_entire()` — is the
signal). Both reproducers are pinned in the lane.

**Known residual, blocked on pf-kh3z:** `Ball::add`/`sub` cannot get the
same treatment because the scalar add/sub *silently* saturate without
flagging (pf-kh3z, the flag-fidelity arc). Once add/sub flag their
saturation, the identical widening pattern applies; recorded here so the
dependency is explicit rather than discovered again.

## Consequences

- `pfloat-ball/tests/regression_review_2026_06_10.rs` pins both classes:
  genuine rim crossings widen with surfaced flags (mul overflow/underflow,
  div overflow), the review's non-crossing input stays exactly what it
  was (an exact, correct ball), and the past-budget parse contains its
  truth on both signs with an in-budget control.
- Enclosure tightness at the rim is sacrificed for soundness (ADR-0087's
  posture: soundness first, tightness measured): an OVERFLOW-saturated
  ball is effectively half-unbounded, and the parse bound is orders of
  magnitude loose. Both are honest about it in their flags.
- Failure modes considered (inverted): (1) widening on flags the scalar
  ops raise spuriously would destroy tightness everywhere — the flags are
  raised only by the documented rim-saturation paths (probed directly for
  the midpoint ops; the radius-pipeline flags were NOT probed in the
  first draft, which is exactly where the verifier found the two holes —
  the lane now pins both, and the verifier confirmed nine ordinary
  mul/div cases bit-for-bit unchanged including Status);
  (2) the parse bound could round DOWN through integer division — the
  per-sign rational choice is the guard, and the lane's witness assert
  would catch an excluding bound; (3) the underflow widening could
  under-reach if the truth's sign differed from the midpoint's — mul/div
  sign is exact, so the truth lies between zero and the clamped value on
  the midpoint's side.

## Related

- Issues: pf-m37w, pf-1bqy (closed by this ADR), epic pf-8iji; pf-mw6u
  (the parse root cause, separate arc), pf-kh3z (the add/sub dependency),
  pf-06lk (the agm consumer of the same scalar rim class).
- Review: `~/.claude/plans/pfloat-workspace-review-2026-06-10.md` Theme 4;
  harness checks BL3/BL4 (BL3's input corrected here).
- Other ADRs: ADR-0077 (ball arithmetic radius soundness), ADR-0086/0087
  (the 1.0 freeze and accuracy posture — this is an internal-behavior
  fix), ADR-0095..0098 (this arc).
