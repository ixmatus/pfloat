# ADR-0074: `pfloat-ball` crate and the `Mag` radius primitive

- **Status**: accepted
- **Date**: 2026-06-06

## Context

Phase 4 is pfloat's rigorous-enclosure tower over the correctly-rounded
scalar base. The first shippable cut is `pfloat-ball`: a real
midpoint-radius ball crate. A ball `[m ± r]` carries a full-precision
pfloat midpoint and a small upward-rounded radius; a ball operation
computes the midpoint with the existing correctly-rounded kernel and
then bounds the radius from the rounding error those kernels already
compute. The ecosystem cell "rigorous enclosure × pure-Rust, no_std,
arbitrary-precision, verifiable" is empty (inari is binary64-only;
rug/Arb/MPFI are C, FFI, LGPL, proof-free), and pfloat already owns the
one asset that makes the tower mostly composition.

The radius needs its own representation. Two soundness hazards must be
impossible to write, not merely avoided by convention: a negative radius
and an inward-rounded (too-small) radius. Either turns the enclosure a
ball certifies into a falsehood. The radius also wants to be cheap
(`Copy`, alloc-free) so it does not dominate the cost of a high-precision
midpoint operation, and `Vec`-free so its invariants can discharge under
Kani (ADR-0062 records that CBMC is hostile to heap allocation, so any
`Vec`-backed value is unverifiable at that level).

Arb's `mag_t` uses a 30-bit mantissa, a choice rooted in 32-bit-FMA
hardware that pfloat has no reason to inherit.

## Decision

Create `pfloat-ball` as a workspace member (mirroring the `pfloat-libm`
member: `pfloat` path dependency with default features off and the ball
crate's own features forwarding the exact pfloat features it needs; crate
attrs `feature(generic_const_exprs)`, `forbid(unsafe_code)`,
`cfg_attr(not(std), no_std)`; a CI job building/linting/testing the
member). The default profile is `std` + `big`; a bare
`--no-default-features` build exposes only `Mag` (alloc-free, the minimal
embedded surface and the Kani target).

Define **`Mag`**, the radius: an unsigned single-limb binary float
`m · 2^(e − 63)` with a `u64` mantissa (top bit set) and an `i64`
exponent, plus `0` and `+∞`. It is an enum `{ Zero, Finite, Infinity }`,
`Copy`, with no `Vec` on its path. The exponent width matches pfloat's,
so a radius and a midpoint exponent compose without a second overflow
regime.

- **Soundness as a type fact.** `{ finite ≥ 0, +∞ }` has no sign field
  and no NaN variant, so a negative or not-a-number radius is
  unrepresentable. Every operation (`add`, `mul`) and every conversion
  into `Mag` (`from_bigfloat_ceil`) rounds the result *up* to the
  single-limb mantissa, so an inward-rounded radius is unrepresentable
  too. `to_bigfloat` is exact (a 64-bit mantissa is representable at
  precision 64) and is the bridge for the exact `lower`/`upper`
  endpoints.
- **64-bit mantissa over Arb's 30-bit.** The single 64-bit limb
  dominates on both axes. Tightness: a wider radius mantissa yields
  tighter enclosures, strengthening the `tightest`/`accurate` story;
  64 bits give a markedly tighter radius than 30. Verifiability: the
  `Copy`, alloc-free, `Vec`-free shape is exactly what lets the round-up
  invariants discharge under Kani. This is the permacomputing-horizon
  lens applied to a representation choice.
- **The `2^-64` resolution cap is a documented, sound choice.** At
  midpoint precision above 64 bits the radius has a relative-resolution
  floor near `2^-64`. This is sound because the radius is only ever an
  upward-rounded upper bound, never an equality: the midpoint carries
  the precision, the radius carries the certified slack, so the radius
  is not the accuracy bottleneck. The crate doc states this so the cap
  is a stated choice, not a latent surprise.

The narrowing strategy (decision from the scoping doc): a ball op forms
the radius bound at the midpoint's precision (the directed spread) and
narrows it upward to `Mag` once at the end, so the only `Mag`-precision
loss is a single final upward rounding.

## Consequences

- Every later slice builds on a radius that cannot be negative or
  inward-rounded by construction, with the round-up arithmetic verified
  to never under-estimate (a BigFloat-oracle grid now, an independent
  exact-rational oracle and an eventual Kani harness in the verification
  slice). An under-estimating radius is the one defect that silently
  breaks every enclosure, so it gets the most verification weight.
- The bare `Mag` surface is no_std, alloc-free, and `thumbv6m`-clean,
  so the embedded and Kani targets are available from slice one.
- `Ball<FixedFloat<PREC>>` is alloc-free only to the degree `FixedFloat`
  is, and `fixed` implies `big` today, so a fully alloc-free embedded
  ball is blocked on the same `[u64; N]` kernel re-implementation
  ADR-0062 scopes. `pfloat-ball` ships `Ball<BigFloat>` as the headline
  and architects nothing that precludes that path later.
- The Arb containment backstop (per-release, `rug`/gmp-mpfr-sys) is not
  wired in this slice to keep the member build dependency-light; it
  lands with the verification slice.

## Related

- Plan: `~/.claude/plans/melodic-knitting-ripple.md` (slice 3); design
  reference `~/.claude/plans/pfloat-phase4-scoping.md` (decisions 2, 6,
  and the "narrow upward once" radius-precision call).
- Beads: `pf-icgj.3` (under epic `pf-icgj`).
- Other ADRs: `scale_by_pow2` (ADR-0072) is used by `to_bigfloat`; the
  `Vec`-free Kani rationale is ADR-0062; the `RealScalar` trait
  (ADR-0075) and `Ball<T>` (ADR-0076) build on this crate.
