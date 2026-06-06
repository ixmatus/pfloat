# ADR-0067: Public constants-at-precision API

- **Status**: accepted
- **Date**: 2026-06-06

## Context

pfloat computes the transcendental constants its kernels need (`π`,
`π/2`, `2/π`, `ln 2`, `ln 10`, `ln 2π`, `2/√π`, the Euler–Mascheroni
constant `γ`) through `pub(crate)` accessors in `src/math/mod.rs`. Each
accessor reads a 1024-bit pinned table for the common precisions and
falls back to the arithmetic-geometric mean (Brent–Salamin for `π`, an
atanh series for the logarithms, Brent–McMillan for `γ`) above the
table ceiling (ADR-0017, ADR-0018). None of this is reachable from
outside the crate.

Roadmap Phase 3 ("pfloat adoption polish") calls for a constants at
precision API so a caller can request a constant at any precision
without learning the `agm` module directly. Two gaps stood in the way.
First, the accessors except `pi_at_round` and `pi_over_2_at_round`
round to nearest only; a public API must honor the caller's IEEE
rounding mode. Second, the one mode aware helper that existed,
`pi_at_round`, rounds a fixed 128 bits above the target and then
rounds down under the mode. That fixed guard was validated only at
`f32` precision by the exhaustive sweep (ADR-0038); it is a heuristic,
not a soundness argument, and is unfit for a public arbitrary precision
surface that claims correct rounding under every mode.

## Decision

Add a public free-function module `pfloat::constants`, gated behind
`exp-log`, exposing `pi`, `pi_over_2`, `pi_over_4`, `two_over_pi`,
`ln_2`, `ln_10`, `euler_gamma`, `two_over_sqrt_pi`, and `ln_2pi`. Each
function has the signature
`fn name(precision: u32, mode: RoundingMode) -> (BigFloat, Status)`,
matching the crate's `*_round` convention (tuple return,
`debug_assert!(precision >= 1)`, `#[must_use]`).

Per-function availability mirrors the cluster feature that already
carries the underlying kernel: `pi`, `pi_over_2`, `pi_over_4`, and
`two_over_pi` need `trig`; `ln_2`, `ln_10`, and `euler_gamma` need
`exp-log`; `two_over_sqrt_pi` and `ln_2pi` need `specials`. No new
feature is introduced.

Each function routes its round-to-nearest accessor through the Ziv
interval-test driver (`crate::math::ziv::ziv_round`) with
`error_guard = DEFAULT_ERROR_GUARD`. The accessor evaluated at a
working precision `w` returns the constant correctly rounded to `w`
bits, whose error against the true value is at most about one unit in
the last place of `w`, which is the premise the interval test assumes.
The driver then certifies the rounding to the target precision under
the caller's mode and escalates `w` until the certification succeeds,
crossing automatically from the 1024-bit table into the AGM path when
the boost pushes past the ceiling. This is the same soundness argument
the rest of the Ziv-driven kernel surface rests on (ADR-0022).

`pi_over_4_at` is added to `src/math/mod.rs` mirroring `pi_over_2_at`:
`π/4` shares `π`'s mantissa, so it is an exact exponent shift of the
table or AGM value.

Scope is the constants that already have kernels. Euler's number `e`
(reachable today as `BigFloat::try_from_i64_exact(1, p)?.exp(mode)`),
Catalan's constant, and Apéry's constant `ζ(3)` are deferred; they
would need new kernels, not just a facade.

## Consequences

The constants are genuinely correctly rounded under all five IEEE
754-2019 rounding modes, not merely round to nearest. `tests/
differential_constants.rs` cross-checks every constant against MPFR
across the four transcendental precisions and all five modes (the
`NearestAway` oracle is synthesized because MPFR lacks it); the unit
tests add a rug-independent check against parsed reference decimals so
the directed-mode rounding is pinned even in feature combinations that
do not build the MPFR lane.

The choice of the Ziv driver over a wider fixed guard buys soundness:
the interval test certifies that the whole evaluation-error interval
rounds to one value, where a fixed guard only makes a double rounding
unlikely. For the irrational constants the first guard (64 bits) wins
in every tested case, so the cost over the fixed-guard path is nil in
practice.

Costs: the public name surface grows by nine functions. Constants
requested above the pinned table recompute through the AGM, memoized
per `(kind, precision)` in a thread-local under `std`. As elsewhere on
the Ziv surface, soundness rests on the interval test terminating;
`ZIV_MAX_ITERS` caps the doubling schedule, the measure-zero caveat
MPFR also documents. A zero precision request is a `debug_assert` and
an `expect` in release, matching every other `*_round` entry point.

No new dependency. The default build stays zero runtime dependencies;
the module is absent unless `exp-log` (or a feature implying it) is on.

## Related

- Plan: `~/.claude/plans/melodic-knitting-ripple.md` (Phase 3, slice C3a)
- Issue: `pf-ebjq`
- Other ADRs: builds on ADR-0017 (AGM transcendental constants),
  ADR-0018 (`γ` via Brent–McMillan), ADR-0022 (Ziv interval-test
  driver), ADR-0038 (`pi_at_round` mode-aware special-case pattern),
  ADR-0039 (per-kernel calibrated error guards)
