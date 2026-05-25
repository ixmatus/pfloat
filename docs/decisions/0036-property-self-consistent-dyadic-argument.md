# ADR-0036: `property_jn::self_consistent` argument constrained to dyadic rationals (pf-ok9 lesson)

- **Status**: accepted
- **Date**: 2026-05-24

## Context

Slice pf-cvs C0 baseline measurement on `bdcf361` ran the full feature
union integration suite locally and surfaced
`tests/property_jn::self_consistent` failing deterministically at
`n=2, num=16, den=3` (`J_2(16/3) ≈ J_2(5.33)`). The persisted proptest
seed reproduced the failure on every re-run; the input lies near
`J_2`'s first positive zero at `x ≈ 5.136`. The initial diagnosis
attributed the failure to a kernel-precision defect in the Miller
backward recurrence near zeros, modeled after slice p1.4's
`BesselJ1` Ziv-envelope fix.

The diagnosis was wrong. Isolation against the assumption found three
inconsistent pieces of evidence:

1. **Per-input bit-exact self-consistency.** Computing `J_2(x_p96)`
   directly at `p=96` returned the same bit pattern as computing
   `J_2(x_p96)` at `p=500` and rounding to `p=96`. The kernel is
   correctly rounded for the input it sees.
2. **Oracle independent certification.**
   `tests/oracle/status/Jn_5.toml` records the MPFR oracle's
   `correctly-rounded` status across `sampled(65536)` f32 inputs with
   `worst_ulp = 0` and `mismatch_count = 0`. A 13-ULP kernel defect
   would have shown up in that sweep.
3. **Parallel proptest consistency.** Three of four `property_jn`
   proptests passed (`boundary`, `parity`, `recurrence`), including
   the DLMF 10.6.1 three-term recurrence cross-tie binding three
   independently descended orders. A real Miller-regime error would
   have broken `recurrence`.

The failure decomposes cleanly when traced. The test does

```rust
let x_lo = rat(num, den, p);
let x_hi = rat(num, den, p + 96);
```

and treats `x_lo` and `x_hi` as the same value at different
precisions. They are not: `rat(num, den, p)` rounds `num/den` to `p`
bits under nearest-even, and for non-dyadic `num/den` the rounded
value depends on `p`. `rat(16, 3, 96) = 16/3 + (1/3)·2^-93` (rounded
up because the fractional ULP position is at 2/3); `rat(16, 3, 192) =
16/3 + (1/3)·2^-189`. The two differ in real magnitude by
`≈ (1/3)·2^-93 ≈ 3.35e-29`. With `|J_2'(16/3)| ≈ 0.322`, the legitimate
output difference is `≈ 1.08e-29 ≈ 6.8 ULPs at p=96`, exceeding the
test's 4-ULP tolerance (`p-2 = 94` bits below the value's binary
exponent).

The amplification scales with `|F'/F|`, which diverges near zeros of
`F`. Any `self_consistent`-shaped test that uses non-dyadic
denominators will fail spuriously near zeros of the function under
test for inputs where the random fractional ULP position lands at the
worst (≈ 1/3 or 2/3) of the grid.

The `pf-ok9 lesson`, already encoded in `tests/property_yn.rs`,
`tests/property_ik.rs`, and `tests/property_zeta.rs`, is to constrain
the argument's denominator to powers of two so `rat(num, den, p)` is
exact at every `p ≥ 1` and both precisions evaluate the same real
point. `property_jn` was the only file in this family that had
`den in 1i64..=4`, which admits `den = 3` (the non-dyadic case the
shrinking proptest found).

## Decision

Constrain `property_jn::self_consistent`'s `den` to
`prop_oneof![Just(1i64), Just(2), Just(4)]`, matching
`property_yn::self_consistent` (the sibling Bessel-Y proptest). Update
the doc comment to explicitly cite the pf-ok9 lesson and the
`|Jn'/Jn|` amplification mechanism, so the next person who reads it
sees the constraint's reason. No production code changes; the kernel
was always correct.

The `recurrence` and `parity` proptests in the same file keep
`den in 1i64..=4`. They evaluate all values at the same precision in
the same invocation, so the input-rounding-dependence cancels: no
`x_lo` vs `x_hi` comparison, no need for the dyadic constraint.

## Consequences

- `property_jn::self_consistent` now passes deterministically, with no
  loss of meaningful coverage. The four-precision-corpora dyadic
  argument set still spans all three Bessel J regimes (tiny `|x|<1`
  via `num=1,den=2,4`; moderate via `num=1..=20,den=1,2,4`; large via
  `num=20,den=1`) and still lands near zeros of `J_n` for various `n`
  at dyadic-representable arguments.
- The `tests/property_jn.proptest-regressions` file created by the
  failing run is deleted, since the persisted seed `n=2, num=16, den=3`
  is now out of range and would never fire. The discovery itself is
  recorded here in the ADR for the next maintainer to find.
- The oracle layer is reaffirmed as the load-bearing kernel-correctness
  verification: `Jn_5.toml` would have caught a real defect. The
  property tests cover precision-self-consistency and DLMF identities,
  not bit-exact correctness; the oracle does the latter.
- **The CI gap that let this slip is identified but not fixed in this
  slice.** CI's `.github/workflows/ci.yml` test matrix runs six
  feature combos per OS, not the full feature union; the bessel
  property tests run only under combos that enable `bessel`, and
  proptest's stochastic draws (12 cases per property, no fixed seed)
  happened to miss the `n=2, num=16, den=3` neighbourhood on every
  per-push run since slice 6o.7 (`fb1b916`) added the test. The full
  feature union integration run belongs in CI or in a documented
  pre-merge gate. That is a separate slice (file as a follow-up bead
  if not already filed; not in pf-jn1y's scope).
- **The general lesson for the proptest template across the codebase**:
  any `self_consistent`-shaped test that compares `F(rat(num,den,p))`
  to `F(rat(num,den,p+k))` must constrain `den` to dyadic. `property_yn`,
  `property_ik`, and `property_zeta` already do. `property_jn` now
  does. If a new `property_*.rs` file adds a `self_consistent`
  proptest in the future, it must inherit this constraint; an
  `ast-grep` rule encoding it is the right shape for a follow-up.

## Related

- Plan: pf-cvs slice plan at
  `~/.claude/plans/prompt-for-fresh-session-delightful-wadler.md`
  (the C0 baseline that surfaced this); pf-jn1y bead description for
  the diagnostic record.
- Commits: TBD on the slice branch.
- Other ADRs: ADR-0022 (Ziv driver, exonerated here as not at fault),
  ADR-0023 (Bessel J kernel design, exonerated here as not at fault),
  ADR-0034 (oracle layer, the kernel-correctness gate that actually
  caught the case across 65536 inputs).
- Beads: pf-jn1y (this slice; closes with this ADR as the deliverable).
  pf-ok9 (the original lesson, encoded in yn/ik/zeta).
