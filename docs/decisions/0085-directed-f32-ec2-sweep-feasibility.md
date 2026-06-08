# ADR-0085: directed f32 EC2 sweep for zeta and Bessel J — feasibility and go/no-go

- **Status**: accepted
- **Date**: 2026-06-07

## Context

The directed-mode rounding verification (pf-3rtr, ADRs 0079 through 0081)
closed the scalar surface to all five IEEE modes, but `pf-3rtr.6` deferred the
f32-grid directed sweep for four hard transcendentals to a future EC2 run:
`zeta`, `J0`, `J1`, `Jn_5` (the status TOMLs carry the deferral note). Their
directed verdicts currently rest on the differential lanes
(`differential_zeta`, `differential_jn`), which compare pfloat to MPFR bit-for-bit
at BigFloat granularity across all five modes at `p <= 256`. The deferred piece
is the f32-grid coverage at exhaustive scale, the directed analogue of the
NE-only `pf-lm3` sweep.

This ADR records the feasibility scoping (no spend authorized): the shard set,
the harness, the per-input cost measured locally, and the go/no-go. (The plan's
"p320 Bessel" phrase referred to the modified Bessel `I`/`K` internal working
precision; the functions actually deferred by `pf-3rtr.6` are the `J` family
plus `zeta`.)

## Decision

**No-go on the exhaustive 2^32 directed f32 EC2 sweep this cycle. Recommend a
targeted local lane instead.** Three findings drive this.

**The harness already exists and runs.** A local release-timed probe confirmed
the pipeline: pfloat `zeta(s)` evaluated at an f32 working precision then
`to_f32_round(mode)` produces directed f32 results that bracket the true value
(spot-checked against an MPFR `Round::Down`/`Round::Up` bracket via `rug`). The
EC2 sweep would reuse the `tests/oracle` MPFR-primary bracket and the
`certified_round_f32` client-side rounding, which already supports all five
modes; NearestAway is synthesized from the bracket, since MPFR has no
`roundTiesToAway`. So no new oracle worker is needed; the gap is purely
coverage, not machinery.

**The cost is steeply input-dependent and the expensive tier is a guard-kill
risk.** Measured per-input cost of `zeta -> to_f32_round` at an 80-bit working
precision:

- Large `|s|` (`>~ 2^10`): about 0.3 us per input (the kernel short-circuits to
  1, or overflows).
- The expensive band `|s|` in roughly `[1, 1000]`: 100 us to 15 ms per input.
  `s` near the pole at 1 costs about 5 ms, and the negative reflection (via the
  functional equation, gamma times sine times the reflected zeta) costs about
  15 ms.

The expensive band is about eleven f32 exponents, about `9.2e7` inputs per
shard. A conservative extrapolation (500 us average over the band) is about 13
core-hours per shard and about 53 core-hours across the four shards: modest
compute (tens of dollars on a fleet), but the near-pole `s` in `[1, 2]` subtier
at 5 to 15 ms per input is exactly the `pf-hzup` failure mode, where a
brute-force shard does not finish inside the run guard and forces a
targeted-tail re-sweep.

**The marginal value is low.** The directed modes are already verified
bit-for-bit against MPFR at BigFloat granularity (all five modes, `p <= 256`).
What the f32 sweep adds is only the `to_f32_round` boundary behaviour at f32
scale, the cases where `zeta(x)` lands within sub-ULP of an f32 rounding
boundary. Those are far better hit by a targeted corpus (f32-rounding-boundary
inputs, the near-pole and functional-equation regions, plus a sampled grid)
than by enumerating `4.3e9` mostly-trivial large-`|s|` inputs.

**Recommendation.** Run a TARGETED directed f32 lane LOCALLY rather than an
exhaustive sweep on EC2: boundary corpus, near-pole and reflection sampling,
f32-rounding-boundary inputs, and a sampled grid, in all five modes, bit-exact
against the certified MPFR bracket. This captures the residual f32-specific risk
at near-zero spend and with no EC2 ceremony, and it sidesteps the near-pole
guard-kill entirely by not enumerating the trivial bulk.

## Consequences

- No spend is requested. The exhaustive EC2 run bead stays deferred; the user
  retains the explicit go decision, now with a grounded cost and a cheaper
  alternative.
- The honest coverage statement is unchanged: zeta and the Bessel `J` directed
  modes are differential-verified at BigFloat granularity; the f32-grid
  exhaustive directed claim remains open and is disclosed as such in the status
  TOMLs.
- If the targeted local lane is built, it is a self-contained piece of work (a
  new bead), reusing the `tests/oracle` bracket and `certified_round_f32`. It
  needs the `differential-mpfr` feature, not EC2.
- If exhaustive f32 directed coverage is later judged worth the spend, the path
  is the `pf-lm3` launcher with the expensive band sub-sharded by exponent (so
  the near-pole subtier gets its own short shards), on-demand not spot, with the
  Noble cloud-init ceremony.

## Related

- Plan: `plans/nested-prancing-lovelace.md` (S6).
- Beads: `pf-6r3o` (this feasibility scoping); the gated cloud-run bead stays
  deferred behind it.
- Other ADRs: ADR-0049 (the `pf-hcz4` / `pf-lm3` sweep infrastructure), ADRs
  0079 through 0081 (the directed-mode verification this completes the f32 view
  of), the `pf-hzup` finding (brute-force infeasible for saturating kernels)
  that the near-pole tier reprises.
