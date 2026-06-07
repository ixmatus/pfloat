# ADR-0078: `pfloat-ball` verification design

- **Status**: accepted
- **Date**: 2026-06-07

## Context

`pfloat-ball`'s whole value is rigorous enclosure: the radius must never
under-estimate. The verification must establish that claim durably (in
repo, runs forever) and honestly (no lane over-sells what it checks). Two
constraints shape the design. First, the soundness-critical primitive
(`Mag` round-up) is `Vec`-free and so Kani-dischargeable, while the
`BigFloat`-backed ball surface is not (CBMC is hostile to heap
allocation, ADR-0062). Second, a self-consistency lane that uses the same
kernel to compute the midpoint and to supply the oracle verifies that the
radius covers the kernel's residual, **not** that the kernel is correct —
that distinction must be stated, not blurred.

## Decision

A tiered verification, with the every-push lanes committed now and the
heavier independent and formal lanes documented with a stated path.

**Tier 1 — types (every push, zero cost).** `Mag`'s `{≥0, +∞}` shape
makes a negative or NaN radius unrepresentable (ADR-0074); `Ball`'s
fallible constructor enforces the finite-midpoint invariant (ADR-0076);
`RealScalar`'s seal makes the midpoint-is-correctly-rounded fact
unbreakable from the crate's own surface (ADR-0075).

**Tier 3 — property tests (every push, pure Rust, blocking).**
`tests/property_ftia.rs` is the blocking FTIA self-consistency lane:
random balls and operations, witnesses reconstructed *exactly* inside the
denoted interval `[mid − rad, mid + rad]` (not the outward-rounded
endpoints), the true `f(witness)` bracketed by the directed scalar kernel
pair at 400 bits, asserted inside the result ball. **Bounded claim,
stated in the file:** because the midpoint and the witness oracle are the
same kernel, this verifies *the radius covers the kernel's own residual*,
not kernel correctness (Phase 1's exhaustive sweep does that). It is the
blocking self-consistency lane, not the independent soundness backstop.
The same file carries the edge cases: degenerate (exact) ball reduces to
the scalar kernel, zero-straddling product contains zero, conversion
preserves containment both directions, entire inputs never panic.

**Tier 2 — proofs (manual Kani lane, advisory).**
`src/kani_harness.rs` discharges the `Vec`-free `Mag` invariants:
round-up monotonicity (`a.add(b) ≥ a, b`), commutativity, the canonical
top-bit-set form on every result, the additive identity / multiplicative
annihilator (the Mag-level exact-in-exact-out), `+∞` absorption (incl.
`0·∞ = +∞`), and the total order across variants. These are the highest
-confidence results the phase can deliver; the universal-mantissa
containment claim over the `BigFloat`-backed ball is blocked on the
shared `[u64; N]` re-implementation ADR-0062 scopes, and is recorded with
ADR-0039's structural-analogy language rather than asserted.

**libFuzzer (60 s smoke per push).** `fuzz/fuzz_targets/ball_parse.rs`
drives the attacker-controlled decimal parser: never panics, never hangs
(its DoS bounds reject pathological literals before any bignum work,
slice 9), and produces a well-formed (`lower ≤ upper`, printable) ball.

## Deferred, with a stated path

These are the independent and formal lanes the v1.0 verification names but
does not yet wire; the one-time Phase 4 adversarial review already
produced independent evidence for each (an exact-rational Python oracle
for `Mag`, an f64 differential for the adjacency primitives, ~660k
high-precision FTIA point-containment checks for the arithmetic and
elementary surface), so the committed lanes are regression guards over an
already-reviewed surface.

- **Independent Arb containment backstop (per-release). LANDED (pf-fe5f,
  this ADR's follow-up).** The self-consistency lane is not independent;
  the independent soundness check is the `differential-arb` lane
  (`pfloat-ball/tests/differential_arb.rs`), which brackets each ball op's
  witnesses with Arb's rigorous interval (the worker's new `BRACKET` verb,
  which emits the exact rational enclosure rather than the rounded f32 that
  `CERTIFY` collapses to) and asserts the result ball is not provably
  disjoint from it. That is the same overlap predicate the self-consistency
  lane uses, so the two test the identical FTIA claim with different
  oracles. The sound direction is the ball admitting Arb's enclosure of the
  true value, never the reverse: a check that the ball lies *inside* Arb's
  interval would pass a too-narrow (unsound) ball, a false backstop. A
  negative control (a quarter-radius ball) confirms the check has teeth.
  **Correction to the original deferral rationale:** the lane reaches Arb
  purely through pfloat's python-flint subprocess and so pulls NO
  `rug`/gmp-mpfr-sys; FLINT/Arb (LGPL) never enter the link graph at rest
  or under test. It is per-release / env-gated (it needs the worker venv).
  The remaining piece is tightness (per-function, per-bucket
  regression-guarded), deferred to pf-fe5f.7: a meaningful slack needs
  Arb's enclosure of f over the input INTERVAL (the witness-image span is
  rounding-dominated for near-exact inputs), a `BRACKET` interval-input
  extension.
- **Lefèvre–Muller hardest-to-round seeding.** Seed the property
  generators with the hard-to-round corpus (hard-to-round scalars produce
  hard-to-enclose balls).
- **Creusot spike on the containment-composition lemma**, strict-revert
  stop-loss; the ADR entry is the deliverable regardless of outcome. Do
  not re-attempt the BoundedBigFloat Kani shim (ADR-0062).

## Consequences

- Every push runs the blocking FTIA self-consistency property and the
  parser fuzz smoke; the manual Kani lane discharges the `Mag`
  invariants. The soundness-critical radius round-up has the most
  verification weight, across types, proofs, properties, and (one-time)
  an independent exact oracle.
- The honest framing is load-bearing: the every-push lane is
  self-consistency, the independent Arb backstop is now wired as a
  per-release lane (pf-fe5f), and the Kani claim is scoped to the
  `Vec`-free primitive. A reader is told exactly what each lane does and
  does not establish.
- `pfloat-ball` is itself the rigorous self-oracle for Phase 1 (pfloat
  verifies pfloat); the independent Arb backstop now keeps that loop from
  being purely self-referential, checking ball enclosures against a second
  library rather than against pfloat alone.

## Related

- Plan: `~/.claude/plans/melodic-knitting-ripple.md` (slice 10); design
  reference `~/.claude/plans/pfloat-phase4-scoping.md` (the three-tier
  verification and the bounded self-consistency claim).
- Beads: `pf-icgj.10` (under epic `pf-icgj`).
- Other ADRs: the `Vec`-free Kani rationale is ADR-0062; the
  structural-analogy language is ADR-0039; the verified primitives are
  ADR-0074 (`Mag`), ADR-0075 (`RealScalar`), ADR-0076 (`Ball`), ADR-0077
  (radius soundness).
