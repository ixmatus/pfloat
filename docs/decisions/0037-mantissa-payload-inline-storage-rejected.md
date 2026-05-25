# ADR-0037: `SmallVec<[u64; 4]>` for `Class::Normal::mantissa` and `Class::Nan::payload`, rejected

- **Status**: rejected (see "Decision" below for the precise scope of
  what was rejected; the data-backed deferral itself is retained)
- **Date**: 2026-05-24

## Context

ADR-0028 measured the allocation baseline on `aarch64-apple-darwin`
and recommended an inline cap of 4 `u64` limbs (256 bits) for
`Class::Normal::mantissa` and `Class::Nan::payload`, targeting the
transcendental and special kernel allocation pattern (gamma p=113
≈ 7820 allocs/op, exp p=256 ≈ 900, mul p=256 ≈ 5). The slice 7f.1
(1.x) work item, tracked as bead `pf-cvs`, was the implementation of
that recommendation.

This slice ran the implementation under CLAUDE.md measurement-led
discipline (read the asm at the call site, strict revert stop-loss on
neutral measurements, ADR is the deliverable either way). The
measurement disproves the ADR-0028 prediction. The slice reverts.

## Decision

The field-only swap with `.into()` wrapping at the rounding-pipeline
boundary, the design the slice plan specified (design choice (a) in
the plan; helpers and workspaces keep `Vec<u64>`), measures
**effectively neutral**. The land criteria require a ≥ 2× allocation
reduction on at least 2 of 3 measured kernels (mul p=256, exp p=256,
gamma p=113). The measurement was:

| kernel       | baseline | field swap | + 3 workspace sites | ratio | meets 2×? |
|--------------|---------:|-----------:|--------------------:|------:|:----------|
| mul p=256    | 5.00     | 5.00       | 4.00                | 1.25× | no        |
| exp p=256    | 915      | 915        | 912                 | 1.003×| no        |
| gamma p=113  | 7880     | 7875       | 6847                | 1.15× | no        |

Two distinct attempts measured:

1. **Field-only swap with `SmallVec::from_vec(storage)`** at the
   destination (the plan's design choice (a), 7 files / +54/-41 lines,
   commit at `78dfb0c` on the slice branch before revert). Saves
   essentially nothing; the alloc count is unchanged within noise.
2. **Field swap plus the three rounding-pipeline workspace allocation
   sites** (`src/rounding.rs:177`, `:243`, `src/big.rs:191` converted
   from `vec![0u64; n]` to `smallvec![0u64; n]`, uncommitted on the
   slice branch). mul drops by exactly 1 alloc/op (the final mantissa
   allocation), gamma drops by ~1000 (rounding-pipeline allocations
   for the p=113 = 2-limb workspaces), exp stays flat (its allocations
   are dominated by something other than the rounding pipeline). Best
   case still misses 2× on every kernel.

The root cause is mechanical, not a measurement artifact. `SmallVec`'s
`From<Vec<T>>` impl is a heap-ownership transfer, not a relocation:
the `Vec` was already heap-allocated by `vec![0u64; n]` upstream, and
`storage.into()` at the destination field merely hands the heap buffer
to the SmallVec rather than copying into inline storage and freeing
the heap. To actually save the heap allocation, the **allocation site
itself** must be `smallvec![0u64; n]` from the start, so that
`n <= 4` uses the inline buffer and `n > 4` spills.

Converting the rounding-pipeline workspace allocation sites (probe 2
above) confirms this: each rounded result that previously allocated
its `target_limbs`-sized `Vec` workspace now allocates inline when
`target_limbs ≤ 4` (p ≤ 256). The mul reduction of exactly 1 alloc/op
matches the count of rounding-pipeline workspaces per multiply.
The gamma 13% reduction (~1000 allocs saved) matches the count of
rounding-pipeline workspaces per gamma call (gamma's Stirling and
constant-composition loop rounds many intermediates).

The remaining allocations live outside the rounding pipeline:

- `multiply_limbs`, `divmod_limbs`, `extract_as_integer`, `isqrt_limbs`
  in `src/ops/limbs.rs` each return `Vec<u64>` results, sized to
  `a.len() + b.len()` or `limbs_for(precision)`.
- Workspace allocations in `src/ops/{fma,div,addsub,mul,sqrt}.rs`
  (~10 sites) sized to operands or to intermediate precision.

The original slice plan (`~/.claude/plans/prompt-for-fresh-session-delightful-wadler.md`)
explicitly excluded these from the slice's scope (design choices (a)
and (b) selected "limbs helpers and workspaces stay Vec"), with the
rationale that at p ≥ 256 these buffers routinely exceed 4 limbs and
would spill anyway, paying a 32-byte inline-header tax on the stack
with no allocation savings. That rationale is correct in isolation;
what the plan missed is that **without those wider conversions, the
field-only swap saves nothing at all**, because the only allocations
that would inline (the small mantissa destination Vecs) are precisely
the ones that `SmallVec::from_vec` cannot relocate.

The hypothesis the slice tested ("the field-only swap with destination
wrapping delivers the ADR-0028 predicted alloc reduction") is
disproved. **The slice reverts per CLAUDE.md strict revert stop-loss.**

What is **not** rejected is the broader ADR-0028 measurement and the
finding that the transcendental and special kernels are allocation-
heavy on a per-call basis. That data is still durable; ADR-0028 is
unchanged. What is rejected is the specific implementation strategy
the pf-cvs slice plan specified, and the implicit assumption that the
field swap is the right granularity for an inline-storage win at the
inline cap of 4.

## Consequences

- **The smallvec dependency is not added to pfloat.** The slice
  branch's C1 (dep add) and C2 (field swap) commits are reset, not
  merged. pfloat retains its zero-runtime-deps posture.
- **`pf-cvs` bead closes with this ADR as the deliverable.** The 7f.1
  (1.x) work item is closed as "attempted, measurement disproved the
  field-level granularity at the inline cap of 4". A future 1.x slice
  could attempt a wider workspace conversion if a fresh measurement
  justifies it, but the slice plan's analysis suggests the wider
  surgery would also not meet a 2× criterion because the limbs-helper
  return values routinely exceed 4 limbs at p ≥ 256 and would spill.
- **Slice 8c (the v1.0 tag) unblocks** without an inline-storage win
  in 1.0. The transcendental kernel allocation pattern remains a known
  measured cost; future perf work on allocation reduction would need a
  different approach (a much larger inline cap accepting the stack
  pressure, a workspace pool or bump allocator, or a structural
  rewrite of the kernel loops to avoid creating intermediate
  `BigFloat`s at every step).
- **ADR-0028 stands.** Its measurement is still the data-of-record
  for the transcendental allocation pattern. Its prediction
  ("inline-cap 4 will deliver a significant drop") is the part this
  ADR contradicts; ADR-0028 should be read alongside this one when the
  inline-storage question comes up again.
- **Lesson for the next slice that touches storage layout.** The
  `.into()` wrap pattern is a no-op for allocation savings; only the
  allocation site matters. Any future inline-storage experiment must
  modify the `vec![...; n]` or `Vec::with_capacity(n)` call sites
  themselves, not the field assignment downstream. The asm spot-check
  at the conversion boundary the slice plan specified (per CLAUDE.md)
  would not have detected this; the alloc-profile measurement did.

## Related

- Plan: `~/.claude/plans/prompt-for-fresh-session-delightful-wadler.md`
  (the pf-cvs slice plan, retained as the historical record of what
  was attempted).
- Commits: none merged (the C1 dep add at `bf6741d` and the C2 field
  swap at `78dfb0c` were reset on the slice branch).
- Other ADRs: ADR-0004 (the original `Vec<u64>` storage choice, with
  the deferral clause this slice attempted to resolve), ADR-0028 (the
  measurement and prediction this slice contradicts), ADR-0029 (the
  parallel data-or-risk-backed v1.0 deferral; same shape but for a
  different opportunity), ADR-0010 (the defer-invasive-perf-past-v1.0
  posture).
- Beads: `pf-cvs` (closes with this ADR as the deliverable).
- Tool: `tools/alloc-profile/` (used for the C0 baseline and C4
  post-measurement; remains in tree for future attempts).
