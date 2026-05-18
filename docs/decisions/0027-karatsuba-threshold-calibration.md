# ADR-0027: Karatsuba threshold calibrated to 48 limbs

- **Status**: accepted
- **Date**: 2026-05-18

## Context

Phase 7 is pure performance. Slice 7d is its only unconditional slice:
calibrate `KARATSUBA_THRESHOLD` in `src/ops/limbs.rs`, the operand limb
count at or below which the multiplication dispatcher (and the Karatsuba
recursion base case) uses schoolbook instead of Karatsuba. The constant
shipped at 30, carrying the comment "match MPFR's default ballpark":
an asymptotic guess, never measured against this implementation.

ADR-0010 fixed the algorithm hierarchy for v1.0 at schoolbook plus
Karatsuba, with Toom-Cook and Schönhage-Strassen deferred to 1.x. The
threshold between the two is the one performance knob in that hierarchy
v1.0 ships, so it is worth setting from a measurement.

The CLAUDE.md performance discipline governs the slice: read the asm at
the call site before patching, because LLVM under thin-LTO and one
codegen unit already inlines small helpers, constant-folds, and unrolls
tight loops, so a patch aimed at work the compiler has done measures
neutral; apply a strict revert stop-loss, where a neutral measurement
reverts the patch and the ADR entry is the deliverable regardless of
outcome.

### Asm read (before any patch)

Disassembly on `aarch64-apple-darwin` at `opt-level=3` (the host the
bench runs on):

- `multiply_limbs_schoolbook`'s inner limb loop is an irreducibly
  serial scalar multiply-accumulate: `mul` plus `umulh` for the
  64×64→128 product, then an `adds`/`cset`/`adds`/`adc` carry chain,
  then a store. The loop carries a dependency through `carry`/`prod`,
  so it is not vectorizable and LLVM has not collapsed it. Schoolbook
  cost is genuine O(n·m) scalar work, not compiler-removable. The bench
  therefore measures a real algorithmic crossover, not a micro-op the
  compiler already performed.
- `multiply_limbs_karatsuba`'s body issues several heap operations per
  recursion node (the `split_at_or_zero` `to_vec` copies, `add_owned`,
  the `vec![0u64; …]` accumulators, the result buffer): on the order of
  17 `__rust_dealloc` plus `__rust_alloc_zeroed` sites per invocation.
  Karatsuba's constant factor here is allocation bound, which pushes
  the empirical crossover for this implementation higher than a
  textbook allocation-free Karatsuba's.

The asm read predicts what the bench then confirmed: 30 is too low for
this implementation.

### Calibration methodology

`benches/mul_thresholds.rs` (criterion, `harness = false`) sweeps
`BigFloat::mul` over 8 to 512 limbs, dense around 30. Operands are
built dense and non-zero so the schoolbook zero-skip fast path is never
taken, and each is self-validated via `Parts::Normal` to occupy exactly
the target limb count. The equal-size group isolates the algorithmic
crossover (the dispatcher's `a.len().min(b.len())` reduces to the
size); the skewed group (n by n/4) exercises the dispatch `min()` and
the unbalanced split.

The dispatcher uses the threshold for both the top-level choice and the
recursion base case, so the crossover was found by measuring the pure
schoolbook curve (threshold transiently forced to 1_000_000, the edit
never committed) against the Karatsuba-enabled curves at thresholds 30
and 48. Settings: warm-up 0.3 s, measurement 1.0 s, 20 samples, plus an
independent confirm run at the decisive sizes. Reported figures are
medians; the cited differences have non-overlapping confidence
intervals.

Equal-size, median time, `aarch64-apple-darwin`:

| limbs | schoolbook | Karatsuba @30 | Karatsuba @48 |
|------:|-----------:|--------------:|--------------:|
| 32    | 986 ns     | 1238 ns       | 992 ns        |
| 40    | 1540 ns    | 1570 ns       | 1558 ns       |
| 48    | 2199 ns    | 2200 ns       | 2318 ns       |
| 56    | 2969 ns    | 2680 ns       | 2623 ns       |
| 64    | 3908 ns    | 3909 ns       | 3359 ns       |
| 128   | 15.7 µs    | 13.0 µs       | 11.0 µs       |
| 256   | 62.5 µs    | 40.2 µs       | 35.2 µs       |
| 512   | 248 µs     | 127 µs        | 108 µs        |

At threshold 30, multiplications in the 32 to 48 limb band dispatch to
Karatsuba and run slower than schoolbook would (n=32: 1238 ns against
986 ns, ~20% slower), the allocation constant the asm read predicted.
Schoolbook stays competitive through ~48 limbs; Karatsuba wins clearly
from ~56 up. The skewed group showed the same shape one rung down via
the `min()`: at threshold 30 the narrower operand in the 31 to 48 band
forced Karatsuba (n=128 skewed: 5.95 µs against 3.64 µs at 48, the
narrower operand 32 limbs).

## Decision

Set `KARATSUBA_THRESHOLD = 48`.

This is a reproducible win against the targeted bench, not a neutral
measurement, so it lands rather than reverting under the stop-loss.
Forty-eight keeps the 32 to 48 limb band on schoolbook (removing the
~20% regression threshold 30 caused there) and still wins the larger
sizes through a better recursion base case (n=512 ~15% faster than at
30, n=256 ~13%, n=64 ~14%).

The value is host- and arch-dependent: it is a measured point for
`aarch64-apple-darwin`, not a portable constant. The doc comment on the
constant states this. Re-calibrating per target is future work, gated
on it mattering; a single honest measured value beats the unmeasured
ballpark it replaces.

No public API change. Toom-Cook and Schönhage-Strassen remain out of
scope per ADR-0010.

## Consequences

- Multiplication in the 32 to 48 limb band is ~20% faster (no longer
  mis-dispatched to allocation-heavy Karatsuba). Larger operands are
  10 to 15% faster via the higher recursion base case. Skewed products
  whose narrower operand falls in the 31 to 48 band gain up to ~39%.
- One honest counter-point: skewed n=384 (narrower operand 96 limbs,
  both thresholds Karatsuba) is ~4% slower at 48 than at 30, an
  isolated loss against many larger wins and a uniformly better
  equal-size curve. It is recorded here rather than hidden; a
  per-target re-tune or a size-aware split could revisit it, deferred
  as not worth the complexity for v1.0.
- The threshold is now a measured artifact with a reproducible bench
  behind it, not a comment citing another library's default. The bench
  stays in tree so the next person (or the 1.x Toom-Cook work) re-runs
  the calibration rather than re-deriving it.
- criterion enters the dev-dependency graph (slice 7d.2). It is
  dev-only and absent from the shipped graph, so the no_std and
  permacomputing-horizon posture is unaffected; `default-features =
  false` drops the plotters and rayon HTML-report tree.
- The 1.x Toom-Cook slice inherits this bench and methodology; its
  schoolbook/Karatsuba/Toom thresholds get calibrated the same way.

## Related

- Plan: `~/.claude/plans/quizzical-prancing-lighthouse.md` (Phase 7).
- Commits: `3732741` (7d.2 bench + criterion dev-dep), `7a4800a`
  (7d.3 threshold 30→48).
- Other ADRs: ADR-0010 (the algorithm hierarchy; Toom-Cook and FFT
  deferred to 1.x, this calibrates the one knob v1.0 ships); ADR-0004
  (the `Vec<u64>` mantissa storage whose per-node allocation is the
  constant factor the asm read found, revisited by slice 7f).
- Bench: `benches/mul_thresholds.rs`.
