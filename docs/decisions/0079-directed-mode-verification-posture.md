# ADR-0079: directed-mode rounding verification posture

- **Status**: accepted
- **Date**: 2026-06-07

## Context

pfloat claims correct rounding under all five IEEE 754-2019 rounding
modes: `NearestEven`, `NearestAway`, `TowardZero`, `TowardPositive`,
`TowardNegative`. Four of those are not the default. They also carry more
weight than their standalone use suggests, because `pfloat-ball` builds
every rigorous enclosure from a directed pair: the lower endpoint is the
operation rounded `TowardNegative`, the upper endpoint the same operation
rounded `TowardPositive` (ADR-0077). A ball is therefore only as sound as
pfloat's directed-mode rounding, and a directed rounding that lands one
ULP on the wrong side of the true value is an unsound enclosure, not a
cosmetic error.

The exhaustive 2^32 binary32 sweep (ADR-0033, ADR-0058) established
correct rounding bit by bit, but only for `NearestEven`: its decimal
oracle bridge rounds to nearest, so the four other modes were never
certified by that lane. The Phase 4 ball review then found six real
shipped directed-mode defects (the `asin`/`atan`/`Si`/`atan2` irrational
constant special cases, the negate-without-mirror class, fixed in 1.2.0),
which confirmed the gap was live rather than theoretical.

A ground-truth pass over the surface before this phase found the position
better than the gap suggested, but still incomplete:

- The `BigFloat` to f32/f64 conversion under a directed mode already
  ships as `to_f32_round` / `to_f64_round` and is cross-checked against
  MPFR over a random off-grid sweep, though the f64 `NearestAway` arm is
  only covered indirectly and there is no boundary-complete corpus.
- The MPFR differential lanes already exercise all five modes at
  `BigFloat` granularity for most of the kernel surface, with the
  `NearestAway` oracle synthesized correctly (MPFR has no roundTiesToAway;
  `MPFR_RNDA` is directed round-away, a different function).
- Every oracle status row records `sampled(65536)`, not exhaustive
  coverage. Twelve functions still carry `unswept` directed rows.
- The composed kernels `log2` and `log10` reach the target precision by a
  fixed guard then a single directed round, with no Ziv interval test
  certifying the directed result.

## Decision

The directed-mode correctness claim rests on three committed corpora, run
every push under the `differential-mpfr` feature, and explicitly not on an
exhaustive directed sweep.

1. **Sampled oracle (per function).** The f32 oracle sweep runs all five
   modes on a `sampled(65536)` input grid, comparing pfloat against an
   independent MPFR (and, for the Arb-primary functions, Arb) certified
   bracket. The directed modes are covered by the mode-aware
   `certified_round_bf_to_f32` bridge, not the nearest-only decimal
   bridge.

2. **All-five-mode differential (per kernel).** Each `differential_*`
   lane compares pfloat against MPFR at `BigFloat` granularity across
   `BIT_EXACT_ROUNDING_MODES`, routing `NearestAway` through the
   synthesized roundTiesToAway oracle. This is the finer-grained lane: it
   checks the full working-precision mantissa, not the f32 projection.

3. **Boundary-complete conversion corpus.** `to_f32_round` /
   `to_f64_round` are checked against MPFR on a corpus that includes every
   structural boundary the conversion can round at: grid neighbours, exact
   midpoints, the subnormal floor, and the overflow threshold, for both
   widths and a direct `NearestAway` oracle.

Correct rounding under directed modes also requires that each kernel be
Ziv-certified per input rather than relying on a fixed working-precision
guard. `log2` and `log10` move onto the `ziv_round` interval test like the
rest of the transcendental surface; the kernels that compute an irrational
constant on a special-case input route through the mode-aware constant
helpers (`signed_constant_at_round`, ADR for 1.2.0) rather than rounding
to target and then negating.

**Exhaustive directed sweep is out of scope, deliberately.** A 2^32
directed sweep across all five modes is the only path to a literal
per-input clean bill, but the sampled oracle, the `BigFloat`-granularity
differential, and the boundary corpus together cover the surface densely
enough that the marginal assurance does not justify the cloud cost. The
decision is revisitable if a future defect is found that all three local
corpora missed.

## Consequences

- The honest claim is "correctly rounded under directed modes as
  established by the sampled oracle, the all-five-mode differential, and
  the boundary conversion corpus, with the residual measure-zero caveat
  below", not "proven correctly rounded on every binary32 input under
  every mode". A reader is told what each lane does and does not
  establish.
- **The residual caveat is real and named.** The Ziv driver caps its
  guard-doubling at `ZIV_MAX_ITERS` (`src/math/ziv.rs`), so on the
  measure-zero inputs whose true value lies within the final guard band of
  a rounding boundary a directed result may be one ULP off. This is the
  same caveat MPFR documents. It is exactly the failure mode the project's
  development disclosure names for numerical code (a directed rounding at
  a boundary that lands a hair on the wrong side), and the sampled corpora
  cannot exclude it.
- Three differential lanes stay deliberately `NearestEven`-only and are
  documented as such, not as gaps: `beta` (a loose two-ULP oracle, so
  bit-exact five-mode agreement is not even its contract), `parse`, and
  `zeta` at p=1024 (cost). They are named here so the exception reads as a
  decision.
- Frugality: once `log2`/`log10` are Ziv-certified and the audited kernels
  are confirmed, the directed-mode claim holds for every future ball
  operation and every directed-rounding caller without re-derivation. The
  verification is entropy reduction over the foundation both the ball's
  soundness and pfloat's own five-mode claim rest on.

## Related

- Plan: `~/.claude/plans/crystalline-seeking-pancake.md`.
- Beads: epic `pf-3rtr` (this phase); `pf-3rtr.1` (this ADR).
- Other ADRs: ADR-0022 (the Ziv interval test), ADR-0038 (the
  all-five-mode differential widening), ADR-0058 (the `NearestEven`
  exhaustive sweep this lane complements), ADR-0077 (the ball directed
  pair whose soundness this protects), ADR-0078 (the ball verification
  tiers). The `log2`/`log10` Ziv decision is recorded separately as the
  next ADR in this phase.
