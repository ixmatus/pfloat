# ADR-0012: Kani harness architecture and CI gating

- **Status**: accepted
- **Date**: 2026-05-10

## Context

Phase 6 lands the verification surface that DESIGN.md scopes:

> Phase 6 lands the harness layout copy-pasted from ferrodec, then
> adapted. Initial properties:
>
> - No panic on bounded-precision inputs for `+`, `−`, `×`, `÷`,
>   `sqrt`, `fma`.
> - Rounding direction is correct under each mode for fixed small
>   precisions.
> - Sign-of-zero correctness across all operations.
> - NaN propagation matches IEEE 754-2019 §6.2.

ferrodec's verification scaffolding (ADR-0009) is the structural
template. ferrodec ships 74 Kani harnesses across `src/verify/<op>.rs`,
organized by op, with a 10-constant operand set (qNaN, sNaN, ±∞, ±0,
±1, ±MAX, ±MIN_POSITIVE) plus `kani::assume()` constraints to keep the
SAT problem tractable.

Two adaptations matter for pfloat:

1. **Operand bounding.** pfloat has no ±MAX / ±MIN_POSITIVE constants:
   it is arbitrary precision with an `i64` exponent. The natural
   substitute is an **eight-constant set** (qNaN, sNaN, ±∞, ±0, ±1,
   neg one) plus two **bounded-normal generators** parameterized over
   exponent range. The generators return non-deterministic finite
   `BigFloat` values constrained by `kani::assume(exp >= lo && exp <= hi)`
   to keep the SAT problem tractable.

2. **CI gating.** ferrodec's CI runs Kani as a blocking job; the
   `feedback_kani_ci_timeout_ok.md` engineering memory records that the
   job times out on every recent GitHub Actions run, and ferrodec
   merges anyway when every other job is green. Importing that
   experience: pfloat's Kani lane is **advisory, not blocking** from
   day one. Failures appear in PR status; they do not block merge.
   Deep verification runs locally on demand.

The non-blocking decision is load-bearing. Treating Kani as advisory
is not "lower standards"; it is matching the gating posture to the
infrastructure's actual reliability. A blocking job that times out
forces every PR author to either rerun (waste minutes) or merge
through a red CI (waste evidence). Advisory output gives Kani's
authority where it earns it (proofs that complete) and routes around
the cost where it does not (proofs that time out).

## Decision

Harnesses live under `src/verify/<op>.rs`, one file per operation,
cfg-gated under `#[cfg(kani)]` at the module level. `src/lib.rs`
declares `#[cfg(kani)] pub mod verify;` after the existing module
list. The module is internal; it has no public surface.

Per-op file shape (mirroring ferrodec's `src/verify/addsub.rs`,
`mul.rs`, `div.rs`, `sqrt.rs`):

```rust
#[cfg(kani)]
mod tests {
    use crate::*;
    use super::super::helpers::*;

    #[kani::proof]
    fn add_nan_propagates() {
        let a = nondet_bigfloat_constant(53);
        let b = nondet_bigfloat_constant(53);
        let (r, status) = a.add(&b, RoundingMode::NearestEven);
        if a.is_nan() || b.is_nan() {
            assert!(r.is_nan());
            if a.is_signaling_nan() || b.is_signaling_nan() {
                assert!(status.contains(Status::INVALID));
            }
        }
    }
}
```

`src/verify/helpers.rs` provides the operand-bounding API:

- `nondet_bigfloat_constant(p) -> BigFloat` — non-deterministic
  choice across the eight-constant set, parameterized by precision.
- `nondet_normal_small(p) -> BigFloat`,
  `nondet_normal_large(p) -> BigFloat` — bounded-normal generators
  for small-exponent and large-exponent ranges.

Rounding-direction harnesses at small fixed precision (PREC=4) land
in slice 6b alongside the arithmetic core. Those use
`#[kani::unwind(N)]` directives with N derived from the precision and
documented inline.

CI invocation, replacing the existing stub in `.github/workflows/ci.yml`:

```yaml
kani:
  name: kani harnesses
  runs-on: ubuntu-latest
  continue-on-error: true
  steps:
    - uses: actions/checkout@v4
    - uses: dtolnay/rust-toolchain@master
      with:
        toolchain: ${{ env.PFLOAT_NIGHTLY }}
    - uses: model-checking/kani-github-action@v1.1
    - run: cargo kani --features=kani
```

The `continue-on-error: true` line is the load-bearing change. The job
runs; its result appears in the dashboard; failures and timeouts do
not block merge.

Local deep verification runbook: `cargo kani --features=kani` runs
the same harnesses without time pressure. Run before significant
arithmetic-kernel changes.

## Consequences

**Wins:**

- Every public op gets harness coverage on its dispatch tree:
  NaN propagation, infinity arithmetic, sign-of-zero, signaling-NaN
  status. Adds a verification layer that proptest cannot deliver
  (proptest finds counterexamples; Kani proves universal properties).
- The advisory posture removes a coupling between Kani infrastructure
  health and pfloat's release cadence. Each release records the
  Kani lane's state; readers can see at a glance whether the proofs
  were green.
- The harness layout copies cleanly across plant-flag projects.
  Future projects can use this ADR as the template for their own
  verification scaffolding.

**Costs:**

- A regressed proof might land if the author does not run Kani
  locally before merge. Mitigation: the slice cadence
  (`feedback_pfloat_slice_cadence.md`) already includes a clippy +
  fmt + thumbv6m sweep; a future revision can add `cargo kani` to
  that list once the verification habits are established.
- The advisory posture is a softer guarantee than ferrodec's
  blocking-but-flaky one. We are explicit about that softness in
  this ADR rather than masking it with a blocking job that
  intermittently fails.
- The eight-constant set is smaller than ferrodec's ten-constant
  set. Two harnesses-worth of coverage (the ±MAX and ±MIN_POSITIVE
  edge cases) do not transfer. The bounded-normal generators are
  intended to cover the same ground; if they fall short, future
  ADRs will revise.

## Trigger to revisit

- A Kani release that reliably completes inside GitHub Actions
  free-tier time limits on pfloat's harness set. At that point the
  `continue-on-error: true` flag can be dropped and the lane
  upgraded to blocking.
- Anthropic-funded paid GitHub Actions runners with longer time
  budgets, which would also enable the upgrade.
- A future plant-flag project whose verification surface is small
  enough that blocking-Kani is feasible from day one; that project's
  ADR may override this one.

## Related

- ADR-0008 (differential testing oracle)
- ADR-0009 (verification scaffolding, copy-paste from ferrodec)
- DESIGN.md, "Verification" section
- `feedback_kani_ci_timeout_ok.md` (engineering memory, load-bearing
  for the non-blocking decision)
- ferrodec/src/verify/* (the template for harness shape)
