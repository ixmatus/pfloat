# ADR-0004: Mantissa storage, `Vec<u64>` for `BigFloat`, `[u64; N]` for `FixedFloat`

- **Status**: accepted
- **Date**: 2026-05-10

## Context

Given the dual API in ADR-0003, each profile picks a mantissa
container.

Candidates considered for `BigFloat`:

- `Vec<u64>` from `alloc`. Standard, predictable, requires `alloc`.
- `SmallVec<[u64; N]>` (or equivalent inline-storage type) from a
  third-party crate or built in-house. Inlines short mantissas
  (≤ N limbs) on the stack, spills to the heap beyond. Saves
  allocations in the common case at the cost of extra type
  complexity and a larger struct.
- A custom thin pointer plus length, hand-rolled. Saves the
  capacity field that `Vec` carries. Marginal storage win, real
  unsafety added.

Candidates considered for `FixedFloat<const PREC: u32>`:

- `[u64; ceil(PREC / 64)]`. Stack-allocated by definition.
- `Box<[u64; ceil(PREC / 64)]>`. Heap-allocated; defeats the point
  of the fixed profile.
- Same `Vec<u64>` as `BigFloat`. Defeats the point.

## Decision

`BigFloat` uses `Vec<u64>` from `alloc`. The `Class::Normal` variant
holds the mantissa directly:

```rust
enum Class {
    Zero { sign: Sign },
    Infinity { sign: Sign },
    Nan { quiet: bool, sign: Sign, payload: NanPayload },
    Normal { sign: Sign, exponent: i64, mantissa: Vec<u64> },
}
```

`FixedFloat<const PREC: u32>` uses `[u64; ceil(PREC / 64)]`,
expressed as a const expression on `PREC`:

```rust
const fn limbs_for(prec: u32) -> usize {
    ((prec as usize) + 63) / 64
}

enum ClassFixed<const PREC: u32> {
    Zero { sign: Sign },
    Infinity { sign: Sign },
    Nan { quiet: bool, sign: Sign, payload: NanPayload },
    Normal { sign: Sign, exponent: i64, mantissa: [u64; limbs_for(PREC)] },
}
```

(Exact spelling depends on which subset of `generic_const_exprs` is
stable on the MSRV; if needed, a wrapper type with a manually
written `Mantissa` impl plays the same role.)

Smallvec-style inline storage for `BigFloat` is deferred. The
profiling data needed to pick the inline cap honestly does not
exist yet. Once Phase 7 lands, if a hot path shows the allocation
cost dominating, this ADR gets revisited.

**Revisited by ADR-0028 (slice 7f.0, 2026-05-18).** Phase 7
measured it: the composing transcendental and special kernels
allocate heavily per call (exp ~900, gamma ~7820 allocs/op), so the
trigger is acknowledged met for that kernel class. The inline-cap
data now exists. The storage change itself (slice 7f.1) is a
crate-wide `Class::Normal` refactor with correctness risk across
every verified kernel, so it is scheduled as data-backed 1.x work
rather than landed against the v1.0 timeline. See ADR-0028.

**Resolved by ADR-0037 (slice pf-cvs, 2026-05-24).** The 1.x work
ran: the `SmallVec<[u64; 4]>` swap at the inline cap ADR-0028
recommended measures effectively neutral (1.25x / 1.003x / 1.15x
allocation reduction on mul / exp / gamma against the 2x land bar),
because `SmallVec::from_vec` is heap-ownership transfer not
relocation; only changing the allocation site (`vec![v; n]` to
`smallvec![v; n]`) saves anything, and even the wider workspace
conversion misses the bar because the intermediates exceed the
inline cap at every `p >= 256`. `Vec<u64>` therefore remains the
chosen mantissa container for both `BigFloat` and the `BigFloat`
form `FixedFloat` converts through.

## Consequences

**Wins:**

- `BigFloat`'s storage is the most predictable shape in Rust:
  `Vec<u64>` clones with `clone`, drops without surprises, walks via
  the standard slice API, and interoperates cleanly with the
  `Mantissa` trait that the kernels read through.
- `FixedFloat`'s storage is on the stack at every precision the
  caller declares. No allocation, no `Drop` overhead, no cache
  misses on the mantissa pointer.
- The `Mantissa` trait abstracts `&[u64]` for `BigFloat::Normal` and
  `&[u64; N]` for `FixedFloat::Normal`; the kernels do not
  distinguish.

**Costs:**

- `BigFloat` pays one allocation per construction. For
  multiplication-heavy workloads at moderate precision, this could
  matter; the smallvec path is the future-work answer if it does.
- `FixedFloat`'s mantissa size is fixed at compile time. Increasing
  precision means a new type instantiation, not a runtime change.
  This is the design intent (see ADR-0003) but is a real constraint
  for callers who think they want fixed-precision and discover they
  do not.
- The const expression `((PREC + 63) / 64) as usize` requires
  `generic_const_exprs` features at a level the MSRV may not yet
  carry. If the workaround is awkward, a sentinel-precision
  enumeration (`P53`, `P113`, `P256`, ...) covers the common cases
  while const generics catch up. The MSRV sets the floor.

## Related

- ADR-0001 (limb representation)
- ADR-0002 (bit-level precision)
- ADR-0003 (dual API)
- DESIGN.md, "Type architecture" section.
