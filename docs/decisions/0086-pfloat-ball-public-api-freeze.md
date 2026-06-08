# ADR-0086: pfloat-ball 1.0 public API freeze

- **Status**: accepted
- **Date**: 2026-06-08

## Context

Slice B0.1 through B0.4 take pfloat-ball from 0.1.0 to 1.0.0. The surface
is the arithmetic-plus-elementary real-ball cut (ADR-0074 through ADR-0078)
plus the verification deepening of ADR-0082 and ADR-0083; no new
construction is part of the 1.0 cut. This ADR is the API stability review
that precedes the tag, the ball analog of ADR-0054 for the scalar crate.
Once 1.0 ships, the public surface is under semver: additive changes only
until 2.0.

The review found the surface already consistent and freeze-ready. Every
ball operation follows one return convention, the construction surface is
uniformly fallible through one error type, and the radius type makes the
load-bearing soundness invariant a type fact rather than a runtime check.
The one open API question carried from Phase 3 (the crate-wide
`NonZeroU32`-versus-`u32` precision call) resolves to "mirror the scalar
crate," for the reason recorded below. The review therefore recommends no
pre-tag reshaping. The changes the ceremony still makes are the version
bump `0.1.0` -> `1.0.0`, the crate-root and README prose flip from the 0.x
disclaimer to the 1.0 commitments (slice B0.4), and the published
enclosure-accuracy posture (slice B0.3, ADR-0087).

## Decision

Freeze the public surface as it stands at this commit. The version bump to
1.0.0 lands at slice B0.4.

### The frozen surface, by feature

- **Always (no_std, alloc-free):** `Mag`. A bare `--no-default-features`
  build exposes `Mag` alone; it is the embedded radius primitive and the
  Kani target.
- **`big` (default):** `Ball<T>`, `BallError`, the sealed `RealScalar`
  trait, `refine_to_accuracy`, `BallParseError`, and the parser bounds
  `MAX_INPUT_BYTES` / `MAX_ABS_EXPONENT`.
- **`fixed`:** `Ball<FixedFloat<PREC>>`, through the `RealScalar` impl for
  `FixedFloat`. `fixed` implies `big` (the embedded alloc-free ball is
  blocked on the same `[u64; N]` kernel reimplementation ADR-0062 scopes;
  the constraint is stated in the crate docs, not worked around here).
- **`exp-log` / `trig`:** the matching elementary methods on `Ball`
  (`exp`-family and hyperbolics under `exp-log`; forward, inverse, and
  reciprocal trig plus `atan2`/`hypot` under `trig`), each gated to the
  pfloat kernel family it wraps.
- **`serde`:** `Serialize` / `Deserialize` for `Mag` and `Ball`.

The arithmetic, accuracy, decimal-I/O, conversion-boundary, and elementary
surface is exposed as inherent methods on `Ball<T>`; the `arith`, `elem`,
`accuracy`, and `io` modules are public only as the home of those impls,
and carry no other public items. `spec.rs` is the in-tree enclosure
contract (prose plus the law statements), the durable specification a 2040
maintainer verifies conformance against without reading Arb's source.

### Precision is `u32`, mirroring the scalar crate

The scalar crate froze precision as `u32` everywhere (ADR-0054:
`try_from_i64_exact(.., precision: u32)`, `precision() -> u32`), with no
`NonZeroU32`. pfloat-ball mirrors this: `Ball::precision() -> u32`,
`refine_to_accuracy(.., start_precision: u32, max_precision: u32, ..)`,
`parse_decimal(.., precision: u32)`, and `to_decimal_interval(digits: u32)`
all take `u32`, with a clamp to at least 1 at the entry points that would
divide by it. A `NonZeroU32` ball over a `u32` scalar would be incoherent:
the midpoint is a pfloat value whose own precision API is `u32`, so the
ball cannot enforce a stronger invariant than the type it wraps. The
crate-wide `NonZeroU32` question therefore stays a single, deferred,
2.0-gated decision across *both* crates, not a divergence the ball
introduces. This closes the open call recorded in the Phase 3 notes.

### Method conventions (locked)

- Every ball operation returns `(Self, Status)`, the ball analog of the
  scalar `(value, Status)` convention (ADR-0007). `Status` is the single
  IEEE 754-2019 sticky-flag set, composed through the existing OR-monoid.
  There is no per-component `Status`: a real ball has one value, and the
  per-component channel is a pfloat-complex concern (a separate crate).
- `INEXACT` on a ball is not a failure signal. A small positive radius is
  the normal, correct outcome of a sound operation, so the radius (read
  through `rel_accuracy_bits` / `abs_error` / `rel_error_bound`) is the
  primary accuracy channel and `Status` is the secondary IEEE-flag channel.
  This is stated so a caller does not treat `INEXACT` as an error the way a
  scalar caller might.
- Construction is uniformly fallible through `BallError`: `new`, `point`,
  and `from_interval` each return `Result<Self, BallError>`. The
  conversion boundary is `lower(&self) -> T` and `upper(&self) -> T` (the
  exact outward-rounded endpoints) and `from_interval(lo, hi)` (the sound,
  never-assume-centered construction). "Lossless" is reserved for the
  ball-to-endpoints direction only; `from_interval` is sound and inflating,
  per ADR-0076.

### Soundness as a type fact (locked)

- `Mag` is `{ finite >= 0, +inf }` with no sign and no NaN, and every `Mag`
  operation rounds toward `+inf` by the type's contract. A negative radius
  and an inward-rounded radius are therefore unrepresentable, not
  runtime-rejected (ADR-0074). `Mag` derives `Clone`, `Copy`, `Debug`,
  `PartialEq`, `Eq`, `PartialOrd`, `Ord`, `Hash` (and serde under the
  feature); the derived order is value order on the canonical form.
- `Ball<T>` carries private fields reached only through accessors
  (`midpoint`, `radius`, `precision`, `is_exact`, `is_entire`). The
  fallible constructors reject a non-finite midpoint and a reversed or
  unordered interval, so a constructed `Ball` always denotes a real
  interval.

### Enums stay exhaustive

`Mag`, `BallError`, and `BallParseError` are exhaustive (no
`#[non_exhaustive]`), mirroring ADR-0054's rationale: the sets are small
and stable, exhaustiveness lets callers `match` without a catch-all, and
the cost of being wrong is a single 2.0-gated breaking change. Revisiting
`#[non_exhaustive]` for the error enums, should a new failure mode appear,
is recorded here as a deliberate post-1.0 question, not an oversight.

### One optional call, ratified: `Mag` carries `Debug`, not `Display`

`Mag` derives `Debug` but not `Display`. Adding a `Display` impl later is a
non-breaking additive change, so the freeze ships without it; a radius is
surfaced to users mainly through the ball's own decimal interval printer
(`to_decimal_interval`) and the accuracy accessors, not by printing the
bare `Mag`. Recorded so the absence is a choice, not an omission.

### Outside the 1.0 stability guarantee

The `differential-arb` feature (the per-release Arb containment plus
BRACKETI lane, reached through the python-flint subprocess and never in the
shipped link graph) and the `kani` feature (the `Mag` proof harnesses) are
verification hooks, off in production builds, and are explicitly **not**
covered by the semver guarantee. The nightly toolchain pin (ADR-0011)
applies as it does to the scalar crate: the MSRV is the pinned nightly in
`rust-toolchain.toml`, because the ball builds pfloat, which needs
`generic_const_exprs`.

## Consequences

- The surface is committed. From 1.0, changes to the items above are
  semver-significant: additive until 2.0. A user moving between the scalar
  and the ball sees one method vocabulary (`(value, Status)`, `u32`
  precision, fallible construction), which is the point of mirroring the
  scalar freeze rather than reshaping.
- The soundness claim the surface makes (every operation returns an
  enclosure that contains the true result over the input ball) is backed by
  the in-tree enclosure spec (`spec.rs`), the blocking FTIA self-consistency
  property lane, and the per-release independent Arb containment plus
  BRACKETI range-soundness lane (ADR-0078, ADR-0082). Tightness is the
  separate, weaker property: measured per bucket and logged, not asserted,
  recorded honestly in the enclosure-accuracy posture (ADR-0087).
- The `NonZeroU32` precision question is closed for 1.0 as "mirror the
  scalar crate (`u32`)"; revisiting it is a cross-crate 2.0 decision.
- `differential-arb` and `kani` consumers take no stability promise, by
  design.

## Related

- ADR-0054: the scalar crate's v1.0 public API freeze, the template this
  mirrors.
- ADR-0007: the `(value, Status)` return convention.
- ADR-0074 through ADR-0078: the `Mag`, `RealScalar`, `Ball`, arithmetic,
  and verification design this freezes.
- ADR-0082: the interval-input Arb bracket (range soundness and tightness)
  the soundness claim leans on.
- ADR-0087: the published enclosure-accuracy posture (slice B0.3).
- ADR-0011: the nightly toolchain pin.
- ADR-0062: the `[u64; N]` reimplementation that the alloc-free embedded
  ball is blocked on.
- Plan: `~/.claude/plans/plan-tower-expansion-scope-goofy-raven.md` (the
  ball 1.0 ceremony, slice B0.2).
