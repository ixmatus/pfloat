# ADR-0011: MSRV moves to nightly to use `generic_const_exprs`

- **Status**: accepted
- **Date**: 2026-05-10

## Context

ADR-0003 commits to a dual API: `BigFloat` with runtime precision and
`FixedFloat<const PREC: u32>` with compile-time precision. ADR-0004
specifies `FixedFloat<PREC>`'s mantissa storage as
`[u64; ((PREC + 63) / 64) as usize]` so the limb count derives from the
bit-level precision parameter. ADR-0002 makes that bit-level precision
non-negotiable for both types.

Stable Rust 1.84 (the original MSRV in ADR-0003) does not stabilize
const-expression evaluation in generic positions. The spelling
`[u64; ((PREC as usize + 63) / 64)]` and the equivalent
`where [(); ((PREC as usize + 63) / 64)]:,` both require
`feature(generic_const_exprs)`, which is unstable. Compiling pfloat on
stable would fail at the FixedFloat declaration site.

The alternatives surveyed during Phase 1 planning were:

1. Keep stable. Replace `FixedFloat<const PREC: u32>` with a sealed
   trait pattern that enumerates the supported precisions
   (`FixedFloat<P53>`, `FixedFloat<P113>`, etc.). Loses bit-level
   precision flexibility for FixedFloat; contradicts ADR-0002 in
   spirit.
2. Keep stable. Change `FixedFloat<const PREC: u32>` to
   `FixedFloat<const LIMBS: usize>` and carry a runtime
   `precision: u32` field. Defeats the point of the type-level
   precision; contradicts ADR-0002 directly.
3. Move to nightly and use `feature(generic_const_exprs)`. Pays the
   nightly cost and the incomplete-feature warning; preserves the
   intended bit-level precision shape.
4. Defer FixedFloat past 1.0. Contradicts ADR-0003's commitment to
   shipping both types in 1.0.

## Decision

Move pfloat's required Rust toolchain from stable 1.84 to nightly,
date-pinned in `rust-toolchain.toml`, and enable
`feature(generic_const_exprs)` at the crate root. Drop `rust-version`
from `Cargo.toml`. Drop the MSRV check job from CI.

The pin format is `nightly-YYYY-MM-DD`, set initially to
`nightly-2026-05-10`. Bumps are deliberate: the toolchain pin moves
when we choose, not when CI refreshes its cache.

This ADR supersedes the MSRV-on-stable stance recorded in ADR-0003's
"Costs" section. ADR-0003's design (dual API, conversion rules) is
otherwise unchanged.

## Consequences

**Wins:**

- The mantissa storage spelling `[u64; ((PREC as usize + 63) / 64)]`
  works directly. No sealed-trait enumeration of precisions, no
  runtime precision field on FixedFloat.
- Bit-level precision (ADR-0002) holds for both types; the user's
  `FixedFloat<53>` is binary64 with correct rounding under any mode,
  exactly as designed.
- Const-generic precision arithmetic in kernels (loop bounds, mask
  computations) is straightforward; the optimizer sees each
  precision instantiation as a distinct concrete type.
- The differential lane (Phase 5+) compares pfloat against MPFR at
  the natural precision granularity. No precision-rounding-up cost.

**Costs:**

- pfloat consumers need a nightly toolchain. This is a real adoption
  tax: CI pipelines that pin to stable cannot use pfloat without
  changes; users who do not control their Rust toolchain (some
  enterprise environments) cannot adopt pfloat at all. Documented
  in the README's "Status" section.
- The `incomplete_features` lint fires at the crate root because
  `generic_const_exprs` is incomplete. The crate root allows this
  lint explicitly; documented at the allow site.
- `feature(generic_const_exprs)` may evolve before stabilization.
  Periodic toolchain bumps may surface compilation errors as the
  feature gates redraw their boundaries. The pin lets us test the
  bump in isolation before merging.
- Some Rust ecosystem tooling (Kani, MIRI, some cargo extensions)
  has historically lagged on the bleeding-edge nightly. We pin to
  a specific nightly date, not "latest nightly," so the ecosystem
  has time to catch up.

**Trigger to revisit:**

- When `generic_const_exprs` (or a stable-mvp subset that covers the
  pfloat use case, currently tracked as `min_generic_const_args` /
  rust-lang/rust#76560) stabilizes on stable Rust, this ADR is
  superseded by a follow-up ADR that moves pfloat back to stable.
- If the toolchain treadmill costs more than the precision flexibility
  buys (monitored across the next four toolchain bumps), this ADR is
  reconsidered.

## Implementation

The following files change in slice 1a:

- `rust-toolchain.toml`: channel `stable` → `nightly-2026-05-10`.
- `Cargo.toml`: remove `rust-version = "1.84"`; add a comment
  pointing at this ADR.
- `.github/workflows/ci.yml`: replace `dtolnay/rust-toolchain@stable`
  with `dtolnay/rust-toolchain@master` plus
  `toolchain: ${{ env.PFLOAT_NIGHTLY }}`; drop the MSRV job.
- `src/lib.rs`: add
  `#![feature(generic_const_exprs)]` and
  `#![allow(incomplete_features)]` at the crate root.

## Related

- Supersedes the MSRV stance in [ADR-0003](0003-dual-api.md).
- Enables [ADR-0002](0002-bit-level-precision.md) and
  [ADR-0004](0004-mantissa-storage.md) at the const-generic site.
- Tracking issue for stabilization: rust-lang/rust#76560.
