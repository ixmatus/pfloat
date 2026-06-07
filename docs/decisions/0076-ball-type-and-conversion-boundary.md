# ADR-0076: `Ball<T>`, the in-tree spec, and the conversion boundary

- **Status**: accepted
- **Date**: 2026-06-07

## Context

With `Mag` (ADR-0074) and `RealScalar` (ADR-0075) in place, the ball type
itself can land: `Ball<T> { mid: T, rad: Mag }`. Two things need pinning
down beyond the struct. First, the enclosure contract — what a ball
denotes and which laws every operation must uphold — must live *in tree*
as a pfloat document, not in a reviewer's memory or implicitly in Arb's
runtime behaviour (the no-human-review posture makes the absent reviewer
and the future maintainer the same party). Second, the
ball-to-endpoints and endpoints-to-ball conversions are not symmetric,
and the asymmetry is a soundness trap: a centred-midpoint
`from_interval` can produce a ball that *excludes* the interval it was
built from.

## Decision

Ship `Ball<T>` with a finite-midpoint invariant: the fallible
constructors (`new`, `point`, `from_interval`) reject a NaN or `±∞`
midpoint (`BallError::NonFiniteMidpoint`). Radius non-negativity needs no
check — it is a `Mag` type fact. `rad = 0` is an exact ball (`{mid}`);
`rad = +∞` is the entire real line. Equality is structural, matching
pfloat's scalar equality.

Write the enclosure contract down in `spec.rs` as five named laws (FTIA
enclosure soundness; the directed-pair radius-soundness law with *no*
Ziv-residual term; exact-in-exact-out; the asymmetric conversion
boundary; status as the secondary channel). The module is prose, the
durable in-tree artifact the verification lanes discharge against.

Implement the conversion boundary per Law 4:

- **`lower`/`upper` are exact.** `lower = mid ⊖ rad` toward `−∞`,
  `upper = mid ⊕ rad` toward `+∞`, via the directed scalar kernels
  (`RealScalar::sub`/`add` with `radius_to_scalar(rad)`). These are the
  tightest representable endpoints; an exact ball collapses both to
  `mid`, an entire ball gives `∓∞`.
- **`from_interval` is sound but inflating.** The midpoint `(lo + hi)/2`
  is computed at working precision (an exact halving of a possibly-rounded
  sum), so it is *not assumed centred*. The radius is
  `rad ≥ round_up(max(mid − lo, hi − mid))`, with each difference bounded
  above toward `+∞` and narrowed up to `Mag` (so it stays an upper bound
  even when the midpoint rounding makes a difference momentarily negative).
  This contains both endpoints unconditionally:
  `lower ≤ mid − (mid − lo) = lo` and `upper ≥ mid + (hi − mid) = hi`.

Unbounded and half-bounded input intervals (an `±∞` endpoint) are
rejected (`BallError::NonFiniteEndpoint`); they are the IEEE 1788
interval face, a separate later crate. A reversed interval (`lo > hi`) is
rejected (`BallError::ReversedInterval`).

## Consequences

- The enclosure law every later slice must satisfy has a single written
  home, so the radius-soundness obligation is checkable rather than
  folkloric. The directed-pair law's "no Ziv-residual term" caveat is
  recorded next to the surfaced-half-width route it must not be conflated
  with.
- `from_interval` cannot silently exclude its input, the classic
  centred-midpoint bug, because the radius is built from the actual
  rounded midpoint, not an assumed-centred one. This is property-tested
  for containment and gets adversarial verification in slice 10.
- The finite-midpoint invariant keeps `lower`/`upper` total and the
  type's denotation always a real interval. Operations that leave the
  reals (division by a zero-containing ball, overflow) are handled in the
  arithmetic slice by widening to an entire-style result plus the IEEE
  status flag, not by storing a non-finite midpoint.
- Reserving "lossless" for the ball-to-endpoints direction (and "sound"
  for the reverse) is now enforced by the code shape, not just the prose.

## Related

- Plan: `~/.claude/plans/melodic-knitting-ripple.md` (slice 5); design
  reference `~/.claude/plans/pfloat-phase4-scoping.md` (the asymmetric
  conversion boundary and the slice-5 `from_interval` ADR obligation).
- Beads: `pf-icgj.5` (under epic `pf-icgj`).
- Other ADRs: consumes `Mag` (ADR-0074) and `RealScalar` (ADR-0075); the
  directed-pair radius-soundness law is implemented by the arithmetic
  kernels (ADR-0077). The directed-kernel correctness it rests on is the
  pfloat `feedback_directed_pair_for_libm` lesson.
