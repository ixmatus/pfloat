# ADR-0084: Creusot containment-composition lemma feasibility spike

- **Status**: accepted
- **Date**: 2026-06-07

## Context

ADR-0078 deferred a Creusot feasibility spike on the containment-composition
lemma: that ball operations preserve the FTIA enclosure because the radius
round-up never under-estimates. Creusot is a deductive verifier (it can prove
functional properties, not just absence of panics), so it could in principle
discharge the soundness-critical `Mag` round-up monotonicity lemma as a total
specification rather than the bounded model-checking Kani provides. This was
filed as a strict-revert, time-boxed probe (pf-6upt): the deliverable is this
ADR, feasible or not, with the precise blocker named. The spike was explicitly
forbidden from re-attempting the `BoundedBigFloat` Kani shim (ADR-0062).

## Decision

**Not feasible to wire in this session. Recorded, not pursued further now.** The
spike's first gate (can Creusot install cleanly without diverging from the
project's pinned toolchain?) fails on three independent counts, established by
probing rather than recall:

1. **Creusot is a separate rustc-driver verifier with its own pinned nightly.**
   Like Kani, Prusti, and MIRAI, Creusot ships a custom `rustc` driver tied to a
   specific nightly revision named in its own `rust-toolchain`. It does not run
   on an arbitrary project-pinned nightly. pfloat pins `nightly-2026-05-10` for
   `generic_const_exprs` (ADR-0011); running Creusot therefore diverges from the
   project toolchain by construction. (This is the same posture Kani already
   has, run as a separate lane.)

2. **Creusot is not packaged in nixpkgs and is not installed here.** A nix
   search surfaces no Creusot package, and neither `creusot`, `cargo-creusot`,
   nor `creusot-rustc` is on the path. Installation means cloning the Creusot
   repository and bringing up its full backend: the Why3 platform plus external
   SMT solvers (Z3, CVC, Alt-Ergo) that discharge the verification conditions.
   That is a substantial external toolchain, not a single `cargo install`, and
   it is exactly the kind of off-orbit dependency chain the project treats with
   skepticism.

3. **`pfloat-ball` enables `generic_const_exprs` crate-wide** (`src/lib.rs`, for
   the `FixedFloat<PREC>` ball path). The soundness-critical `Mag` is a single
   `u64` limb and uses no const generics, but the feature gate is on the whole
   crate, so a Creusot translation of `pfloat-ball` must contend with
   `generic_const_exprs`, an incomplete feature that analysis backends do not
   reliably support. Even a `Mag`-only target would need `Mag` factored into a
   crate that does not enable the feature.

Beyond the install gate, the marginal value is low. The `Mag` round-up
monotonicity lemma is `Vec`-free and is already discharged by the Kani lane
(ADR-0078, Tier 2). The genuinely unproven property, the containment-composition
lemma over the `BigFloat`-backed ball surface, runs the real `Vec`-backed
arithmetic and so hits the same heap wall that ADR-0062 documents for CBMC; a
deductive tool reaches that surface through the same `Vec`, so Creusot is not a
way around the wall ADR-0062 named.

## Consequences

- No code was written; there is nothing to revert. The strict-revert stop-loss
  resolved at the first gate, as intended for a feasibility probe.
- The honest verification posture is unchanged: `Mag` invariants by Kani
  (proofs), the ball surface by the property and independent-Arb lanes
  (ADR-0078, ADR-0082). Creusot is not currently a fit.
- The future path, if Creusot is pursued, is named: factor `Mag` (and the
  radius round-up lemma) into a `generic_const_exprs`-free crate, then run
  Creusot against that crate under Creusot's own toolchain as a separate lane,
  the way Kani already runs. This is a self-contained piece of work, not a
  precondition for any shipped behavior.

## Related

- Plan: `plans/nested-prancing-lovelace.md` (S5).
- Beads: `pf-6upt` (discovered from `pf-fe5f`).
- Other ADRs: ADR-0078 (the deferral and the Kani `Mag` lane), ADR-0062 (the
  `Vec` heap wall a deductive tool does not avoid), ADR-0011 (the
  `generic_const_exprs` toolchain pin).
- References: the Creusot project (rustc-driver verifier translating to Coma /
  Why3, nightly-pinned), <https://github.com/creusot-rs/creusot>.
