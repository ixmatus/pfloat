# ADR-0068: Optional serde serialization behind a feature

- **Status**: accepted
- **Date**: 2026-06-06

## Context

pfloat froze its public API at v1.0 (ADR-0054) with no serialization
support. Roadmap Phase 3 ("pfloat adoption polish") calls for serde
interop so downstream consumers can persist, cache, and transmit
values (configuration, snapshots, inter-process messages) without
hand-rolling a codec.

Two constraints shape the design. First, the crate ships zero runtime
dependencies, and the default build must keep that posture; any serde
support is therefore opt-in. Second, deserialization runs on bytes the
caller did not produce, a named vulnerability surface, so it is a trust
boundary: a malformed encoding must not be able to construct a
`BigFloat` that violates the internal canonical-form invariants the
arithmetic kernels rely on.

## Decision

Add an optional `serde` dependency (`default-features = false`,
`features = ["derive", "alloc"]`) behind a new `serde` feature
(`serde = ["dep:serde"]`). The default build is unchanged and stays
no_std clean.

The simple public types (`Sign`, `Status`, `RoundingMode`, `IeeeClass`,
`BuildError`, `ParseError`) carry a `#[cfg_attr(feature = "serde",
derive(...))]` at their definition. `BigFloat` and `FixedFloat<PREC>`
are hand-written in `src/serde_impls.rs`.

The wire form is a `{ precision, class }` struct where `class` is an
externally-tagged enum mirroring the internal `Class` variants (`Zero`,
`Infinity`, `Nan`, `Normal`) carrying the raw sign, exponent, and
mantissa or payload limbs. It is exact by construction and does not
branch on `is_human_readable`, so a round trip recovers the value to
the bit under any serde format, human readable or compact.
Serialization borrows the limbs (no clone); deserialization owns them.

Deserialization revalidates the canonical form rather than trusting the
input: precision at least one bit, the `Normal` mantissa limb count
equal to `limbs_for(precision)`, the top bit of the most significant
limb set, the storage-padding bits below the precision clear, and the
NaN payload limb count equal to `limbs_for(precision)`. Malformed input
is rejected with a serde error, never silently coerced.
`FixedFloat<PREC>` serializes through its `BigFloat` view and
deserializes through `try_from_big_exact`, which rejects any precision
other than `PREC`.

The `serde` feature is added to the CI test-matrix union and the clippy
union (kept byte-identical per `feature-union-check.sh`), so the impls
are linted and the round-trip and rejection tests run on every push.

## Consequences

Downstream crates gain serde interop. The round trip is exact across
every `Class` variant and both precision profiles, verified in
`tests/serde_roundtrip.rs` (eleven tests, six of them confirming the
deserialize trust boundary rejects malformed encodings). The trust
boundary is explicit and tested rather than implicit.

The default build keeps its zero runtime dependencies; `serde` enters
the graph only when a caller opts in, and the `alloc` feature keeps it
no_std clean. `serde_json` is a new dev dependency, used only by the
test.

The wire form mirrors the internal `Class` shape, so a future change to
the representation would be a serde-format break. The form is
documented in `src/serde_impls.rs`; a versioned envelope can be layered
on later if the representation ever needs to evolve without breaking
stored data. Because the impls are format agnostic, one format (JSON in
the test) validates them and compact codecs (postcard, bincode) inherit
the same structure for free.

`serde` and `serde_derive` are the canonical ecosystem serialization
crates; the dependency is justified by the interop the feature exists to
provide, and it is optional, so the foundation-crate posture holds for
every caller that does not need it.

## Related

- Plan: `~/.claude/plans/melodic-knitting-ripple.md` (Phase 3, slice C3b)
- Issue: `pf-c9ok`
- Other ADRs: builds on ADR-0054 (v1.0 public API freeze), ADR-0001 and
  ADR-0002 (mantissa canonical form and bit-level precision), ADR-0005
  (`Class` tagged union), ADR-0006 (`i64` exponent)
