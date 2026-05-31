# ADR-0053: CI gate widening for v1.0 (oracle smoke coverage and a feature-union drift guard)

- **Status**: accepted
- **Date**: 2026-05-30

## Context

Bead pf-xyaq was filed to "widen CI to gate the full-feature-union
integration test run," on the stated rationale that the per-push matrix
ran a 6-combo feature subset and not the union, and that this gap "let
pf-jn1y and pf-1axr slip in." Closing the bead before slice 8c required
first establishing what the gap actually is. The git history does not
support the original framing.

The test matrix in `.github/workflows/ci.yml` has, since the kernel
slices, carried a seventh entry that is the full library-feature union:

```yaml
- "--features=std,fmt,big,fixed,ops,exp-log,trig,specials,agm,integrals,airy,bessel,zeta"
```

`cargo test` on that entry compiles and runs every `tests/property_*.rs`
file, because each is gated only by a file-top `#![cfg(...)]` on the
features it needs, all of which the union enables. The union entry was
extended kernel by kernel as each landed: `bessel` was added to it at
slice 6o.6, before `property_jn` shipped at slice 6o.7, and the entry
predates both pf-jn1y (2026-05-24) and pf-1axr (2026-05-27). So
`property_jn` and its siblings have run per-push the whole time.

Neither implicated defect was a feature-coverage gap:

- **pf-jn1y** (ADR-0036) was a latent proptest flake. `property_jn`'s
  `self_consistent` case treated `rat(num, den, p)` and
  `rat(num, den, p+96)` as the same value, which holds for dyadic `den`
  but not otherwise; the failing input was only sometimes sampled. The
  test ran in CI; the flake hid in the random draw, not behind a
  disabled feature.
- **pf-1axr** (ADR-0042) was a too-narrow input sweep. `differential_yn`
  swept random `|x|` in 1..40, never reaching the `|x| >= 128`
  asymptotic regime where the `bessel_y` range-cap pre-check returned
  `NaN`. The lane ran per-push; the boundary input was simply never
  generated until slice 2b.2.a added it to the dyadic table.

The genuine residual gap is narrower and real. The MPFR-backed oracle
smoke tests (`oracle_types_smoke`, `oracle_mpfr_smoke`,
`oracle_verify_smoke`, `oracle_driver_smoke`, `oracle_smoke_gate`,
`oracle_certified_round_bf_to_f32`) require `differential-mpfr` to
compile, so the test matrix never builds them, and the dedicated
`differential` job runs `cargo test ... --test 'differential_*'`, whose
glob does not match `oracle_*`. No per-push job ran them at all; they
were exercised only by a local full-feature run.

A second, deeper hazard sits behind the union entry: it is a hardcoded
list, duplicated in the `clippy` job. A new kernel feature added to
`Cargo.toml` but forgotten in those two lists would silently stop being
tested and linted per-push, reopening exactly the kind of hole the
union entry was meant to prevent. That hazard had been a convention
held in reviewer memory; a convention is the weakest tier of
verification.

## Decision

Two changes, both per-push, both cheap.

1. **Oracle smoke coverage, via the full MPFR feature union.** Drop the
   `--test 'differential_*'` filter from the existing `differential`
   job's run step, leaving:

   ```yaml
   - run: cargo test --release --features=differential-mpfr --verbose
   ```

   A `--test 'oracle_*'` glob does not work: cargo treats a globbed test
   selection as explicit, so it errors on `oracle_arb_midpoint_smoke`
   and `oracle_cross_check_smoke` (which `require differential-arb`)
   rather than skipping them. With no `--test` filter at all, cargo runs
   every target whose `required-features` are met and silently skips the
   rest. That is the matrix's own behaviour, and it is exactly what is
   wanted: the differential lanes, the MPFR oracle smokes, and the
   property/unit suite all run once under the MPFR feature union, while
   the Arb- and `ziv-instrumented`-gated oracle targets stay skipped for
   their EC2 cross-check tier (ADR-0049, with its Python and Arb system
   dependencies). The differential sweep is unchanged and still runs
   once; the added property/oracle/unit time is small against it, well
   inside the job's 90-minute budget. The job is renamed from "MPFR
   differential" to "MPFR full-union integration" to match.

2. **Feature-union drift guard.** Add `scripts/feature-union-check.sh`
   (POSIX `sh`, the `conformance-counts.sh` pattern) and wire it into
   the cheap `conformance` job. It asserts two facts: the test-matrix
   union entry and the `clippy` job's `--features=` list are
   byte-identical, and every `[features]` key in `Cargo.toml` either
   appears in that union or is on an explicit, commented exclusion
   allowlist (`default`, `alloc`, `kani`, `ziv-instrumented`,
   `differential-mpfr`, `differential-arb`). A new kernel feature then
   fails CI until it is added to the union or deliberately excluded.

The slice also corrects the pf-xyaq paragraph in `docs/ROADMAP.md` to
record the accurate gap rather than the original misdiagnosis.

## Consequences

- The six MPFR oracle smoke tests now run on every push, on the Linux
  differential runner. They are smoke-scale, so the added wall-clock is
  negligible against the differential sweep already in that job.
- The drift guard moves a convention from reviewer memory into a
  runnable check, one tier below a type and above documentation. This
  is the same per-bucket-not-aggregate discipline `conformance-counts.sh`
  embodies, applied to feature coverage instead of harness counts.
- The guard is intentionally strict: adding a feature is now a
  three-line edit (Cargo.toml, the matrix entry, the clippy line) or a
  deliberate allowlist entry. The friction is the feature; it is what
  makes silent drift impossible.
- The Arb / cross-check oracle tier remains gated to the EC2 ceremony
  (ADR-0049), not per-push, because of its system dependencies.

## Related

- ADR-0014: MPFR differential CI gating (the `differential` job this
  extends).
- ADR-0036: pf-jn1y, the proptest flake the original framing
  misattributed.
- ADR-0042: pf-1axr, the narrow-sweep defect the original framing
  misattributed.
- ADR-0049: the EC2 cross-check sweep tier the Arb oracle tests belong
  to.
- `scripts/conformance-counts.sh`: the per-bucket gate this guard is
  modelled on.
