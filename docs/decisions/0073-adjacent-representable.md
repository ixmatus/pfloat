# ADR-0073: Adjacent-representable primitives (`next_up`/`next_down`/`ulp`)

- **Status**: accepted
- **Date**: 2026-06-06

## Context

Phase 4 (`pfloat-ball`) needs the adjacent-representable value and the
unit in the last place in three places: `ulp(mid)` in the
`from_interval` soundness bound, the radius construction, and (later)
the IEEE 1788 `accurate`-tier `nextOut` and the tightest-endpoint
conversion. These are also broadly useful (interval stepping, error
analysis, `Float`-trait conformance), so they belong in pfloat core,
not in the ball crate. The v1.0 surface (frozen under ADR-0054) has no
`nextUp`/`nextDown`/`ulp`.

The semantics are IEEE 754-2019 §5.3.1 (`nextUp`/`nextDown`) plus the
common `ulp`, but pfloat's representation forces the boundary behaviour
to be derived afresh rather than copied from a binary64 implementation.
pfloat has **no subnormals and no `emin`/`emax`**: a finite value is
`±m·2^(e−p+1)` with `m` a `p`-bit top-bit-set integer and `e` an `i64`
that *saturates*. The exponent floor and ceiling are where the adjacent
operations meet `±0` and `±∞`, and those transitions differ from a
bounded IEEE format:

- The smallest positive value is `MinPos = 2^(i64::MIN)` (top-bit-only
  mantissa, `e = i64::MIN`), reached in a single step from `±0` (no
  gradual-underflow ladder of subnormals).
- The largest finite magnitude is `MaxFinite` (all-ones mantissa,
  `e = i64::MAX`); its successor is `+∞` because the exponent cannot
  increment past the ceiling.
- The true `ulp` of a value near the floor (`e − p + 1 < i64::MIN`) is
  below `MinPos` and therefore not representable.

For a rigorous-enclosure consumer the soundness direction of every
boundary choice is load-bearing: a radius built from an `ulp` that
*under*-estimates the real gap is unsound.

## Decision

Add `next_up`, `next_down`, `ulp` (each `(&self) -> (Self, Status)`,
with `_with_flags` siblings) to `BigFloat` in `src/ops/adjacent.rs`, and
the same signatures to `FixedFloat<PREC>` (delegating through
`to_big`). Purely additive; the v1.0 freeze is preserved. The result
always carries `self`'s precision.

Semantics, derived from §5.3.1 adapted to the saturating-exponent model:

- **`next_up(x)`** = least representable value greater than `x`:
  `+∞ → +∞`; `−∞ → −MaxFinite`; `±0 → +MinPos`; `−MinPos → −0` (the
  zero's sign is negative, fixed by §5.3.1); `MaxFinite → +∞`; a
  positive finite adds one ulp (a mantissa carry crosses a power of two
  upward, incrementing `e`); a negative finite subtracts one ulp toward
  zero (a borrow crosses a power of two downward, decrementing `e`).
- **`next_down(x)` = `−next_up(−x)`**, implemented by that identity.
  Negation is exact and preserves NaN payload/signaling, so the
  negative-side boundaries fall out of the positive ones and any
  `INVALID` is raised exactly once.
- **`ulp(x)`** = the positive `2^(e−p+1)`: `ulp(±0) = MinPos`,
  `ulp(±∞) = +∞`, `ulp(NaN) = NaN`, always sign-positive. When
  `e − p + 1 < i64::MIN` the true ulp is below `MinPos`, so `ulp`
  **saturates upward** to `MinPos` and raises `UNDERFLOW`. `e − p + 1`
  cannot exceed `i64::MAX`, so there is no upward overflow.
- A **signaling NaN** raises `INVALID` and returns a quiet NaN for all
  three (`nextUp`/`nextDown` are general-computational; `ulp` follows
  the same rule). A quiet NaN passes through with `OK`.

## Consequences

- The ball's `from_interval` bound and radius construction get exact
  adjacency and a sound `ulp`. The one inexact case (`ulp` near the
  exponent floor) saturates **upward**, the sound direction for a
  radius: an over-estimated ulp keeps the enclosure valid, and the
  `UNDERFLOW` flag lets a caller detect the measure-zero saturation
  rather than silently trusting a clamped value.
- `next_up`/`next_down` are total and panic-free across the whole
  representation, including the `±∞`/`±0` transitions that a bounded
  IEEE format reaches via `emax`/subnormals. The `MaxFinite → +∞` and
  `−MinPos → −0` choices are the saturating-model analogues of the IEEE
  boundary rules, written down here rather than left implicit.
- No raw-parts constructor is introduced: the boundary values
  (`MinPos`, `MaxFinite`) are assembled internally from validated
  canonical mantissas, so the top-bit-set normalization and storage
  shape stay type-checked invariants.
- `next_down(x) = −next_up(−x)` keeps the two directions in lockstep by
  construction; they cannot drift, and only one code path carries the
  boundary logic.

## Related

- Plan: `~/.claude/plans/melodic-knitting-ripple.md` (slice 2); design
  reference `~/.claude/plans/pfloat-phase4-scoping.md`.
- Beads: `pf-icgj.2` (under epic `pf-icgj`).
- Other ADRs: `scale_by_pow2` (ADR-0072) is the sibling slice-1
  primitive and the means to construct `2^k` for the `ulp`
  cross-check; the no-`emax` saturating exponent originates in
  `ops::mul`; `parts()`-without-a-converse-constructor is ADR-0016.
- Verification: an independent f64 differential oracle (pfloat at
  precision 53 vs Rust `f64::next_up`/`next_down`/`ulp` over the
  binary64 normal range) plus a five-lens adversarial review of the
  saturating-exponent boundaries.
