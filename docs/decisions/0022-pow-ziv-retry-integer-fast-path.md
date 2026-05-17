# ADR-0022: `pow` Ziv interval-test retry and integer-exponent fast path

- **Status**: accepted
- **Date**: 2026-05-16

## Context

`pow` shipped in slice 3c as a single evaluation of `exp(y · ln(x))`
at a fixed `target + 64` guard, with no correct-rounding retry. The
composition accumulates rounding from `ln`, the multiply, and `exp`,
so the result is not correctly rounded: `tests/differential_pow.rs`
carried a 2 ULP tolerance and ran NearestEven only, and ADR-0014's
slice-6h status update recorded this as structural limitation #3,
naming the fix as "either a Ziv-strategy retry pass or an
integer-exponent fast path." The roadmap's slice 7c does both. It is
sequenced before the Bessel slices because their recurrence
machinery wants a working integer `pow` primitive.

This is the first Ziv loop in the codebase. `exp`, `ln`, `si`, and
`erf` all use the fixed 64-bit guard. DESIGN.md §"Ziv's strategy"
already states the contract: compute at `target + guard`, detect the
rounding-boundary case, double the guard and retry, cap the
iteration count and document it. pfloat has no error-bound plumbing
through `ln`/`mul`/`exp`, and building it is out of scope and
against the frugality posture, so the realization had to work from
the values the existing kernels already return.

## Decision

**Two evaluation paths feed one Ziv driver.** `integer_exponent`
reconstructs an exact integer exponent (sibling of the existing
`integer_parity`, same `scale = exponent − precision + 1`
decomposition, value built from the lowest to the highest set bit so
a small integer survives the long trailing-zero mantissa it carries
at high precision). When the exponent is an exact in-range integer,
`pow_int` forms `x^|n|` by square-and-multiply and the driver
reciprocates for `n < 0`; otherwise the driver evaluates
`exp(y · ln(x))`. The negative-base branch already delegated to
`pow_positive(|x|, …)`, so the fast path covers negative integer
bases with no extra wiring. A result-exponent feasibility guard
(`POW_INT_RESULT_EXPONENT_CAP`, applied where the base is in scope
since `integer_exponent` does not see it) and an
overflow/underflow check defer the extreme cases to the `exp·ln`
path, which carries the correct OVERFLOW/UNDERFLOW status.

**The Ziv test is the interval test, not recompute-and-compare.**
The first cut compared the target-rounded value at two adjacent
guards and returned on agreement. That is a heuristic, not the Ziv
criterion: on a hard-to-round input both insufficient guards agree
on the same wrong value and the loop converges falsely. The
bit-exact differential lane caught this on the first run
(`pow(63, -3)` at p=53). The shipped `pow_ziv` instead bounds
`eval`'s error at the working precision by the half-width
`|y| · 2^-(working − ZIV_ERROR_GUARD)` (`ZIV_ERROR_GUARD = 24`,
comfortably above the accumulated NearestEven rounding of `pow_int`'s
at most 64 multiplies or the `exp·ln` path); if both ends of that
uncertainty interval round to the same target value under the
caller's mode, the true value rounds there too and the result is
correctly rounded. Otherwise a rounding boundary lies inside the
interval, the guard doubles (schedule 64, 128, 256, 512, 1024), and
the loop retries, bounded by `ZIV_MAX_ITERS = 5`. The half-width is
formed by the exact power-of-two exponent decrement already used in
`math::mod::pi_over_2_at`.

**The iteration cap is the honest caveat.** On the measure-zero
exact-tie inputs that exhaust `ZIV_MAX_ITERS` the result may be 1
ULP off in directed modes, the same caveat MPFR documents. It is
stated in the `pow_round` doc comment per the DESIGN.md contract.

**Differential lane tightened to bit-exact across all five modes.**
`tests/differential_pow.rs` now asserts exact equality against MPFR.
The four MPFR-mappable modes use `mpfr_round_of` directly.
NearestAway has no MPFR oracle: `MPFR_RNDA` is directed round-away
(the farther neighbour of any inexact value) and `MPFR_RNDN` is
ties-to-even, while an integer base to a small integer power is
frequently an exact tie at the target precision (`99^8` sits exactly
between two p=53 values), so neither is roundTiesToAway. The lane
synthesizes roundTiesToAway from a `p + 128` MPFR result and rounds
it itself. The latent `mpfr_round_of` NearestAway mapping (unused
elsewhere, since every other lane is NearestEven only) is filed as
separate cleanup, not widened here.

**No new Kani harness.** `src/verify/pow.rs` (the IEEE 754-2019
§9.2.1 special-case table) is unchanged. Ziv correct-rounding is a
precision property over unbounded-precision loops, not amenable to
bounded model checking; it is pinned by the bit-exact differential
lane plus the new property tests. **No new fuzz target**:
`fuzz/fuzz_targets/exp_log_family.rs` already exercises
`x.pow(&y, mode)`.

## Consequences

`pow` is correctly rounded under every IEEE rounding mode for
positive base and finite exponent, subject to the documented Ziv
cap, and is the first transcendental off the NearestEven-only
differential tier. ADR-0014 limitation #3 is closed (limitations #1
and #2 were closed earlier by slices 7a/ADR-0016 and 7b/ADR-0017).
Exact integer powers settle bit-exactly at the first guard, matching
MPFR's integer fast-path parity.

Costs and honest tradeoffs. The Ziv driver recomputes when the
interval straddles a boundary; the common path is one evaluation
plus a cheap interval check, so there is no routine doubling, but a
genuinely hard input pays up to five evaluations. No performance
machinery was added on top of this without a bench, per the
project's measurement-first rule. The error model assumes the
accumulated rounding of `eval` stays under `2^24` ULP at the working
precision; that holds for the domain this kernel serves (bounded
multiply counts, the bounded `exp·ln` composition) and is the same
class of assumption MPFR's own guard heuristics make. The
`pow_ziv` driver and the half-width helper are written so a later
slice can reuse them for other transcendentals when those move off
the NearestEven-only tier.

## Related

- Plan: `plans/polymorphic-growing-nygaard.md`
- Commits: `b84ea1b` `a34e57b` `db8aa29` `9bea08e` `8ed6e89`
  `e7ee70d`
- Other ADRs: closes limitation #3 of ADR-0014; follows ADR-0016
  (slice 7a bit-exact converter) and ADR-0017 (slice 7b AGM
  constants), which closed limitations #1 and #2.
