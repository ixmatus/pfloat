# ADR-0087: pfloat-ball enclosure-accuracy posture

- **Status**: accepted
- **Date**: 2026-06-08

## Context

Slice B0.3 of the ball 1.0 ceremony publishes the ball's accuracy
declaration, the artifact ADR-0086 points at for the tightness dimension it
deliberately leaves out of the freeze. The scalar crate publishes a
per-function correct-rounding table (`docs/rounding-status.md`, generated
from the oracle status records and CI-checked by
`scripts/rounding-status-table.sh`). The ball needs the analog, but a ball
makes a different kind of claim and the analog cannot be a copy.

A scalar function returns one value, so the published property is a point
verdict: correctly rounded, or not. A ball operation returns an
*enclosure* of the true result over the whole input ball, so the property
splits in two. Soundness (does the ball contain the truth) is the hard,
uniform guarantee. Accuracy (how tight the sound ball is) is a separate,
graded, quality property. Conflating them is the classic interval-arithmetic
mistake; the ball's verification design already keeps them apart (ADR-0078:
containment is blocking and gets the strongest lanes, tightness is logged).
This ADR records how the published posture reflects that split, and how it
is gated against rot.

## Decision

Publish `docs/ball-enclosure-status.md`: a per-operation table of the
enclosure shape and the expected accuracy class, with the soundness-versus-
accuracy framing stated up front.

### Two layers, stated explicitly

- **Soundness** is the hard guarantee, uniform across every operation in
  the table: the output ball contains the true result over the input ball
  (the Fundamental Theorem, `src/spec.rs` Law 1), backed by the blocking
  FTIA self-consistency lane and the per-release independent Arb
  containment plus BRACKETI range-soundness lane (ADR-0078, ADR-0082). The
  posture table does not restate this per row; it holds for all of them.
- **Accuracy** is the posture column, and it is **measured, not asserted**.
  The `differential_arb` tightness lane logs the ball width against Arb's
  rigorous image per `(function, precision, magnitude)` bucket
  (`pfloat-ball/tests/arb_tightness_expected.txt`), regression-guarded per
  bucket, not on an aggregate floor (the compensating-regression rule). The
  posture column records the tightness class that *follows from each
  operation's enclosure shape*; the measured lane is the source of truth a
  regression trips.

### The per-operation classes

Two classes, keyed to the enclosure shapes in `src/elem.rs` /
`src/arith.rs`:

- **tightest** — the directed-pair arithmetic (`add`, `sub`, `mul`, `div`,
  `sqrt`, `cbrt`: correctly-rounded midpoint, radius bounding the kernel
  residual under ADR-0077) and the monotone-endpoint elementary functions
  (the `exp`/`log` family, the hyperbolics and their inverses, the inverse
  trig: the enclosure is the directed-kernel-correct endpoints of a
  monotone function, the tightest representable interval). Residual width
  beyond this is the inherent input-propagation term, not a rounding loss.
- **accurate** — the 1-Lipschitz functions (`sin`, `cos`, `hypot`) and the
  composed functions (`tan` = sin/cos, `atan2` by gradient bound, `cosh`
  through a magnitude interval). These are sound and tight in the
  near-linear or monotone region, but the shape's variation bound
  conservatively over-covers near a local extremum, where `|f'| -> 0` but
  the radius still carries the full input radius. The straddles where Arb's
  image is tightest are exactly where this radius is loosest; that
  `|f'| -> 0` effect is the one ADR-0082 records, and the concrete reason
  tightness is measured rather than asserted as uniform `tightest`.

No `valid`-only operation ships in the v1.0 surface: every operation
reaches at least `accurate`. The honest caveat is the conservatism of the
1-Lipschitz and composed families near extrema, not a soundness gap.

### The CI-checked form

The posture is authored, not generated: the enclosure shape and the class
are design declarations, not values read off a sweep, so unlike
`rounding-status-table.sh` the gate does not regenerate the doc. It gates
the doc's **completeness** instead, the role `feature-union-check.sh` plays
for the CI feature union: `scripts/ball-enclosure-status.sh --check`
asserts that the set of backtick-quoted operations in the table rows equals
the `pub fn` surface of `arith.rs` and `elem.rs`. A new operation added
without a posture row, or a row left behind after a rename, fails the
`conformance` job. `src/spec.rs` gains a closing pointer to the posture doc
so the soundness contract and the accuracy declaration cross-reference.

## Consequences

- The ball ships a published accuracy declaration that does not overclaim:
  it states `tightest` only where the midpoint is correctly rounded or the
  endpoints are monotone-exact, and `accurate` with a named conservatism
  elsewhere, with tightness measured by the lane rather than asserted.
- The completeness gate makes the declaration a maintained contract: the
  surface cannot grow a v1.x operation (or, later, a ball special function)
  without declaring its posture, the same anti-rot discipline the scalar
  conformance gates carry.
- A reader who wants the soundness obligations reads `src/spec.rs`; a reader
  who wants the accuracy expectation reads `docs/ball-enclosure-status.md`;
  the two are linked, and neither restates the other.

## Related

- ADR-0086: the v1.0 API freeze, which defers the tightness dimension here.
- ADR-0077: the directed-pair radius-soundness law behind the `tightest`
  arithmetic class.
- ADR-0078: the verification design that separates blocking containment
  from logged tightness.
- ADR-0082: the interval-input Arb bracket and the `|f'| -> 0` tightness
  finding the `accurate` class records.
- `scripts/rounding-status-table.sh` and `scripts/feature-union-check.sh`:
  the scalar-side precedents for the published table and the completeness
  gate this mirrors.
- Plan: `~/.claude/plans/plan-tower-expansion-scope-goofy-raven.md` (slice
  B0.3).
