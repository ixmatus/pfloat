# ADR-0028: Allocation profiling, and `BigFloat` inline storage deferred to 1.x with data

- **Status**: accepted; the inline-storage prediction is contradicted by ADR-0037
- **Date**: 2026-05-18

**Outcome recorded by ADR-0037 (slice pf-cvs, 2026-05-24).** The
measurement and scheduling this ADR specifies both stand; the
allocation pattern numbers below remain the data of record. What
ADR-0037 contradicts is the implicit prediction that swapping the
mantissa container to `SmallVec<[u64; 4]>` at inline-cap 4 would
deliver a significant per-call allocation reduction. pf-cvs ran the
swap, measured 1.25x / 1.003x / 1.15x against the 2x land bar, and
reverted: `SmallVec::from_vec` does not relocate, and the wider
workspace conversion does not meet the bar either because the
intermediates exceed the inline cap at every `p >= 256`. Future
inline-storage work needs either a much larger inline cap (at
real stack-pressure cost), a workspace pool / bump allocator, or
a structural rewrite of the kernel loops to avoid creating
intermediate `BigFloat`s at every step. See ADR-0037 for the
measurement table and the lesson.

## Context

ADR-0004 chose `Vec<u64>` for the `BigFloat` mantissa and deferred
smallvec-style inline storage with an explicit trigger: "Once Phase 7
lands, if a hot path shows the allocation cost dominating, this ADR
gets revisited." No allocation profiling was done in Phase 6, so the
trigger was unmeasured: slice 7f could neither proceed nor be honestly
deferred. The Phase 7 plan resolved this by measuring first, then
deciding the disposition with the data.

### The unsafe-free invariant constrained the measurement design

A zero-dependency allocation profiler needs a custom
`#[global_allocator]`, which is necessarily `unsafe impl GlobalAlloc`.
`pfloat`'s `Cargo.toml` sets `[lints.rust] unsafe_code = "forbid"`
package-wide, and that applies to examples and benches too. `forbid`
cannot be locally overridden by `#[allow(unsafe_code)]` (that is
itself an error under `forbid`). So neither the plan's original sketch
(a `cfg`-gated counting allocator in the library) nor the obvious
refinement (a pfloat example binary) is possible without weakening a
load-bearing security invariant.

Four options were considered:

1. A standalone in-tree tool crate (`tools/alloc-profile/`, its own
   workspace root, path dependency on pfloat, its own lint scope).
   Preserves pfloat's `forbid(unsafe_code)` fully; adds no dependency
   to pfloat; the harness stays committed in-tree.
2. Relax pfloat's package lint from `forbid` to `deny` and
   `#[allow(unsafe_code)]` in the harness. Weakens the unsafe-free
   invariant for every pfloat target, not just the harness.
3. A `dhat` dev-dependency. The plan explicitly chose zero-dependency
   over `dhat` for frugality.
4. Drop the measurement and defer 7f with the trigger still
   unmeasured. Leaves the data gap ADR-0004 named.

Option 1 was chosen (user decision at slice time). It is the only
option that both measures and keeps the security posture intact, at
the cost of one small standalone crate. The counting allocator
forwards every call verbatim to `System` with a stated `SAFETY`
justification confined to that crate.

## Decision

### The measurement

`tools/alloc-profile/` counts heap allocations per operation on three
representative kernels (release, `aarch64-apple-darwin`), warmed up to
exclude one-time lazy constant initialisation:

| kernel       | allocs/op | bytes/op | mantissa |
|--------------|----------:|---------:|---------:|
| `mul` p=256  | ~5        | 224 B    | 4 limbs  |
| `exp` p=256  | ~900      | ~48 KB   | 4 limbs  |
| `gamma` p=113| ~7820     | ~352 KB  | 2 limbs  |

Arithmetic is allocation-light. The composing transcendental and
special kernels allocate a fresh `Vec` for essentially every
intermediate `BigFloat` in their Taylor and Stirling loops, hundreds
to thousands per call, while each mantissa is tiny (2 to 4 limbs, 16
to 32 bytes). This is precisely the many-small-allocations pattern
inline storage targets.

### The disposition

The ADR-0004 trigger ("a hot path shows the allocation cost
dominating") is acknowledged **met** for the transcendental and
special kernel class. The inline-cap data now exists.

The storage change itself (slice 7f.1) is **deferred to 1.x as
concrete, data-backed work**, not landed against the v1.0 timeline.
Rationale, a Parnell judgment call parallel to the 7e deferral
(ADR-0029):

- 7f.1 is a crate-wide change to `Class::Normal`'s mantissa container
  touching every kernel and every helper that constructs a
  `BigFloat`. The library is verified (Kani, differential vs MPFR,
  property tests); a storage refactor risks correctness regressions
  across the whole surface right before the v1.0 tag.
- v1.0's value is a stable, correct, complete surface. A measured
  performance opportunity that needs an invasive refactor is exactly
  the kind of work the roadmap defers past the tag (the ADR-0010
  posture for Toom-Cook and FFT).
- The deferral is now **data-backed, not hand-waved**. ADR-0004 asked
  for the measurement before revisiting; the measurement is done,
  recorded, and turns the deferral from "unmeasured, unknown" into
  "measured, scheduled". That is the outcome ADR-0004 wanted.

### Inline-cap guidance for 7f.1 (1.x)

The measured mantissas are 2 to 4 limbs. An inline cap of 4 `u64`
(256 bits) keeps binary64, binary128, and the common p=256 working
precision off the heap, covering the exp and gamma allocation classes
profiled here, while the heap spill handles arbitrary precision. 7f.1
should re-run `tools/alloc-profile/` before and after to confirm the
allocation count drops, and re-run the full verification matrix and
the 7d `mul_thresholds` bench (inline storage shifts the Karatsuba
constant factor, so ADR-0027's threshold may need re-calibration).

## Consequences

- `tools/alloc-profile/` is a permanent in-tree measurement artifact.
  The next person re-runs it rather than re-deriving the question; it
  is the 7f.1 before/after harness.
- pfloat's package-wide `unsafe_code = "forbid"` is preserved exactly.
  The cost is one standalone crate outside the package; it is not
  built by pfloat CI (a manual tool, like the 7d bench).
- 7f.1 is scheduled 1.x work with a concrete scope and an inline-cap
  starting point, tracked in the issue graph (`discovered-from` the
  7f decision). It is not relitigated from zero in 1.x.
- ADR-0004's deferral clause is resolved: its status now points here.
- Phase 8 slice 8b cites this ADR in the conformance evidence block as
  a measured, scheduled deferral, alongside ADR-0027 (7d) and ADR-0029
  (7e).

## Related

- Plan: `~/.claude/plans/quizzical-prancing-lighthouse.md` (Phase 7).
- Commits: `757e31b` (7f.0 tool crate).
- Other ADRs: ADR-0004 (revisited here; the `Vec<u64>` choice and its
  deferral clause), ADR-0003 (the dual API the storage choice serves),
  ADR-0027 (7d; the Karatsuba constant 7f.1 would shift), ADR-0029
  (7e; the parallel data-or-risk-backed v1.0 deferral), ADR-0010 (the
  defer-invasive-perf-past-v1.0 posture).
- Tool: `tools/alloc-profile/`.
