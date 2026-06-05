# ADR-0062: The BoundedBigFloat Kani discharge of Ziv soundness is blocked at the Vec allocation level, not the arithmetic

Status: Accepted (2026-06-04)

## Context

The Ziv interval-test soundness theorem (`src/verify/ziv_soundness.rs`,
pf-hdh8) is discharged under Kani only over the canonical eight-constant
operand set (qNaN, sNaN, the infinities, signed zero, plus or minus one).
ADR-0039 deferred the lift to true universal quantification over the
mantissa domain, sketching a `BoundedBigFloat<80>` fixed-array shadow as
the encoding that would make CBMC's symbolic execution tractable: a
`[u64; 80]` mantissa has a statically known length where pfloat's
`Vec<u64>` mantissa does not, so every limb loop CBMC unwinds would be
bounded. pf-25zw set out to build that shadow and lift the discharge.

This ADR records what the investigation measured. The deferral was
qualitative ("Vec is hostile to CBMC"); it is now quantitative, and the
diagnosis is sharper than expected: the wall is not the arithmetic, it is
the `Vec` allocation itself.

## Decision

**The fixed-array-shim-into-the-real-ops approach cannot discharge. Do
not land the `BoundedBigFloat` shadow type.** The genuine path is a fully
`Vec`-free re-implementation of the operations, which this ADR scopes as
the open follow-up.

### What was measured (Kani 0.67.0, CBMC, aarch64-apple-darwin)

1. **The existing real-op harness unwinds without bound.** Running the
   simplest harness (`ziv_interval_test_zero_half_width_is_trivially_sound`,
   one eight-constant operand, zero half-width) with no unwind bound,
   CBMC unwinds the limb comparison loop in `cmp::limb_cmp_aligned`
   (`for offset in 1..=len`, `len` from `Vec::len()`) past 1600 iterations
   and never terminates. The `Vec` length is symbolic to CBMC, so the loop
   is unbounded.

2. **Unwind bounds make it finite but not tractable.** With
   `--default-unwind 8` the limb loops bound correctly (CBMC reports
   "not unwinding" at the cap), so the problem is finite. It still did not
   discharge within a six minute budget on that simplest harness: the
   bit-precise CBMC modelling of the multi-limb `Vec` arithmetic
   (`add` / `sub` / `round_to_precision`) is too large.

3. **Even a copy-and-compare round-trip through `Vec` fails.** A
   `BoundedBigFloat<N>` shadow plus conversion shims (`to_bigfloat` /
   `from_bigfloat`) and a round-trip lemma `from(to(x)) == x` was built.
   The lemma performs no arithmetic, only a fixed-array copy into a `Vec`
   and back. It still reported a verification failure, accompanied by
   "dereference failure: dead object" and `NonNull` dereference checks
   left undetermined: CBMC's model of the `Vec` heap allocation breaks. A
   native concrete-playback of CBMC's own counterexample does not panic,
   confirming the shim is correct and the failure is a CBMC modelling
   artifact of the `Vec`, not a bug in the round-trip.

### Why the shim approach cannot work

A `BoundedBigFloat<80>` operand fixes the *input* length, but the Ziv
soundness theorem is a property of the rounding and comparison logic, and
evaluating it runs the real `add` / `sub` / `round_to_precision` /
`partial_cmp`. Those operate on `Vec<u64>`. Converting the bounded operand
to a `BigFloat` re-introduces the `Vec` at the first step, and measurement
(3) shows CBMC fails on the `Vec` at the allocation level, before any
arithmetic. So no amount of operand bounding or unwind tuning rescues a
shim that calls the real, `Vec`-backed operations.

### The genuine path, deferred

A discharge requires the four operations the theorem needs
(`add`, `sub`, `round_to_precision`, `partial_cmp`) re-implemented on
`[u64; N]` fixed arrays with no `Vec`, so CBMC never allocates a heap
object. Fidelity to the real kernel would then be established separately,
by differential testing the fixed-array operations against the real ones
over a large random corpus (CBMC cannot check the fidelity itself, since
that comparison would run the real `Vec` ops). This is a substantial
effort (the four operations plus the rounding logic, re-derived on fixed
arrays, plus the fidelity corpus), and it verifies a faithful shadow
rather than the kernel directly. It is scoped here as the open follow-up;
it was not pursued in pf-25zw.

## Consequences

- ADR-0039's deferral stands, now on a reproduced, quantified boundary
  rather than an assumption. The eight-constant discharge plus the pf-tqzz
  per-release sweep cross-check remain the stand-in for the universal
  claim, exactly as ADR-0039 records.
- ADR-0012's qualitative "Vec is hostile to CBMC's symbolic execution" is
  upgraded: the hostility is at the heap-allocation level (a copy-and-
  compare round-trip fails), not only in the data-dependent loop bounds.
  This is the load-bearing fact for any future verification work on
  pfloat: a tractable Kani discharge of a `Vec`-backed kernel needs a
  `Vec`-free shadow of every operation on the path, not merely bounded
  operands.
- No code lands. The `ziv_soundness` module doc is updated to point here
  so the next person does not re-attempt the shim.

## Related

- pf-25zw (this investigation); pf-hdh8 (the scaffolded theorem).
- ADR-0039 (the deferral this measures), ADR-0012 (the CBMC-on-Vec lesson
  this sharpens).
