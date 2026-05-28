# ADR-0049: pf-hcz4 full pf-tqzz cross-check sweep — first execution and v1.0 baseline

- **Status**: in-progress (results pending sweep run)
- **Date**: 2026-05-28

## Context

ADR-0039 closed Phase 1g with the pf-tqzz cross-check assertion shipped
as a per-push smoke (245 assertions at `tests/oracle_cross_check_smoke.rs`)
and the full per-release sweep deferred. The bound the smoke validates is

```
|eval(working) − midpoint| ≤ 2^(error_guard − working) · |midpoint|
```

evaluated in `rug::Float` at `oracle_prec = working + 64`, with
`error_guard` read from per-kernel calibrated constants in
`src/math/ziv_calibration.rs`. All 38 distinct `<KERNEL>_ERROR_GUARD`
constants currently sit at `DEFAULT_ERROR_GUARD = 24` (algebraic-primary
calibration per pf-yupm); the cross-check sweep is the empirical
secondary that confirms the algebraic claim on every `f32` input.

Phase 2 perf work (ADR-0040 through ADR-0048) closed across 2026-05-27
and 2026-05-28; the full sweep is the next v1.0-chain gate per
[[project_perf_before_full_sweep]]. The bead `pf-hcz4` records the
acceptance criteria: extend or land a full-sweep harness, run it,
handle violations, produce a durable per-kernel pass/fail summary in
`tests/oracle/status/`.

The sweep scale: 65 536 `f32` inputs × 5 IEEE rounding modes ×
47 v1.0 FnIds (63 (FnId, order) shards once parametric Bessel orders
unroll) = **15.4M assertions**. Single-threaded wall-clock projection
post-Phase-2 perf: 3–6 hours per slowest shard (Bessel small-arg at
`verification_precision = 320`). Per-FnId sharding spreads the cost.

## Decision

Four pieces shipped under this ADR:

### 1. Shared cross-check module

Factor `cross_check_one` + `error_guard_for` + `midpoint_for` +
`CheckOutcome` + the new `ViolationRecord` type out of the smoke
harness into `tests/oracle/cross_check.rs`. Gated on
`#[cfg(all(feature = "differential-arb", feature = "ziv-instrumented"))]`.
The smoke (`tests/oracle_cross_check_smoke.rs`) and the new sweep
example (`examples/pf_tqzz_sweep.rs`) both consume this module via
the existing `mod oracle;` / `#[path = "../tests/oracle/mod.rs"]` pattern.

The shared `cross_check_one` returns `CheckOutcome::Violation(v)` instead
of panicking; the smoke's caller wraps the variant in `panic!` to
preserve its pre-refactor structured-panic UX, and the sweep's caller
appends each `v` to a sidecar `Vec<ViolationRecord>`.

### 2. Full-sweep CLI binary

`examples/pf_tqzz_sweep.rs` models on `examples/oracle_sweep.rs`:
`--fn-id <name>` (canonical `FnId::name()` or parametric `Jn:5` syntax),
`--modes all|RNE,...`, `--sample <N>` (default 65 536), `--output <PATH>`
(default `/tmp/pf_tqzz_<fn>.json`), `--skip-lm-seeds`,
`--instance-type <S>` (informational tag). Default input grid: the
Lefèvre-Muller hard-to-round corpus for the FnId (when one exists) plus
`(0u32..sample)` `f32` bit patterns, deduplicated.

The binary emits a single JSON sidecar per shard. Schema documented in
the file header; aggregator consumes it.

### 3. EC2-sharded execution ceremony

Per the ferrodec one-unit-per-instance template
([[project_fd_ykr_campaign_closed]]):

- `scripts/pf-hcz4-launch.sh` — idempotent prereqs (bucket, default
  VPC, latest Ubuntu Noble ARM64 AMI via SSM Parameter Store, quota
  check); discovers 63 shards from `ls tests/oracle/status/*.toml`;
  per-shard `aws ec2 run-instances` with cloud-init baking `FN_ID`,
  `RUN_ID`, `S3_BUCKET`, `GIT_SHA`; tags `RunId` + `Shard` + `Project`;
  schedules `shutdown -h +180` SSM fallback after launch.
- `scripts/pf-hcz4-status.sh` — tag-filtered instance discovery;
  bulk `aws ssm send-command` for `tail -n 1 /var/log/sweep.log`;
  checks `_DONE` + `_FAILED_PREFLIGHT` sentinels; auto-`s3 sync` on
  all-done.
- `scripts/pf-hcz4-relaunch.sh` — single-shard recovery with
  `MARKET=on-demand` default for spot reclaim cases.
- `scripts/pf-hcz4-aggregate.py` — reads per-shard `result.json`;
  emits global summary + per-FnId `[cross_check]` table appended to
  `tests/oracle/status/<fn>.toml`.

Default instance: `c8g.large` spot (~$0.03/h), 2 vCPU Graviton4.
Projected cost: 63 × 3h × $0.03 ≈ $5.67, well under the $10 ceiling.

Critical ferrodec-template details baked in:

- Cloud-init opens with `exec > >(tee -a /var/log/sweep.log) 2>&1`
  (not plain `>`); console-output stays useful when SSH is unavailable.
- `pip install awscli` rather than `apt install awscli` (apt awscli
  is flaky on Ubuntu Noble).
- Per-instance `--instance-initiated-shutdown-behavior terminate`
  plus an explicit `shutdown -h now` at the end of user-data.
- Defensive `shutdown -h +180` via SSM after launch, bounding billing
  if user-data hangs.
- Each cloud-init runs the per-push smoke first as a pre-flight gate;
  failure uploads `_FAILED_PREFLIGHT` instead of `_DONE` so the
  aggregator can distinguish setup vs sweep failure.

### 4. Collect-all-violations failure handling

Violations do NOT panic. The sweep runs to completion across every
input × mode for the given FnId. Each violating triple lands in the
shard's JSON with `(input_u32, mode, working_prec, eval_w, midpoint,
abs_diff, bound, gap, ratio_log2)`. The aggregator sorts all
violations across all shards by `ratio_log2` descending — the priority
signal for triage.

Triage protocol if violations land: NO auto-widening of
`*_ERROR_GUARD` constants. Per fault class, file a follow-up bead
`bd create --deps "discovered-from:pf-hcz4"` and a `feedback_*.md`
memory if the pattern generalises. The pf-hcz4 slice ships the
measurement; widening is a separate slice with its own ADR.

## Results

*Filled post-run. Placeholder structure below.*

- **Run ID**: `<RUN_ID>`
- **Tip SHA**: `<GIT_SHA>` (main at slice start: `9453215`)
- **Shards launched**: 63 (full surface) + 1 smoke pre-flight
- **Total assertions**: TBD (`63 × 65 536+lm-seeds × 5` ≈ 15.4M+)
- **Total violations**: TBD
- **Per-kernel violation counts**: TBD
- **Top violations by `ratio_log2`**: TBD
- **Wall-clock**: max-shard TBD, sum TBD
- **AWS cost**: TBD (target <$10)
- **Spot reclaims**: TBD
- **Durable artifact**: `tests/oracle/status/<fn>.toml` `[cross_check]`
  tables (63 files); per-shard `result.json` archived to
  `s3://${S3_BUCKET}/${RUN_ID}/`

### If zero violations

All 38 per-kernel `*_ERROR_GUARD` constants confirmed empirically across
the swept input surface. The v1.0 "actively guarded" claim per ADR-0039
strengthens from algebraic-only to algebraic-plus-empirical.

### If non-zero violations

Per fault class:

- Largest-`ratio_log2` cell named in the follow-up bead.
- Bead description carries top-5 violations for the kernel +
  smallest power-of-two `error_guard` that would have passed.
- Pattern crosswalks (e.g. "every cell at `working_prec ≥ 400`")
  go to `feedback_*.md`.

## Consequences

- **v1.0 baseline established.** The 63 `[cross_check]` tables in
  `tests/oracle/status/` become the durable evidence the slice 8c
  release ceremony attaches to the v1.0 announcement.
- **Per-release cadence locked in.** ADR-0039's per-release gate is
  no longer hypothetical; this is its first execution. Subsequent
  v1.x releases re-run the same sweep against `pf_tqzz_sweep`'s
  CLI surface; no harness changes per release expected.
- **EC2 ceremony reusable.** The `scripts/pf-hcz4-*.sh` set
  generalises to any future cross-check campaign (the assertion
  surface lives in `oracle::cross_check`; only the CLI dispatch
  is FnId-specific). A similar campaign against an expanded
  IEEE-binary64 grid would re-use the same launcher with
  `--sample` adjusted.
- **AWS-account-specific configuration.** `S3_BUCKET`,
  `IAM_INSTANCE_PROFILE`, `AWS_REGION` are env-var-supplied; the
  scripts do not bake account identifiers. Re-running on a different
  account requires only `S3_BUCKET` + `IAM_INSTANCE_PROFILE`
  re-provisioning (the IAM role needs `s3:PutObject` on the bucket
  and `AmazonSSMManagedInstanceCore` for the `aws ssm send-command`
  bulk-tail).
- **No source-of-truth changes.** This slice does not modify
  `src/math/ziv_calibration.rs`. Any constant-widening lives in a
  separate follow-up slice if violations surface.

## Alternatives considered

- **Local single-instance run on `aarch64-apple-darwin`.** Considered;
  rejected because the 6–12-hour run blocks the dev box and the
  cost-vs-disruption math favours the sharded EC2 fan-out (≈$6 vs
  half a day of local iteration loss). The script ceremony cost
  amortises against future per-release sweeps.
- **Macro-generated `#[test] fn sweep_<fn_id>` per FnId.** Considered;
  rejected because cargo-test's string-prefix filter does not map
  cleanly to parametric Bessel orders (`Jn:5` vs `Jn:10`), and the
  bead's acceptance criterion calls for a structured JSON artifact,
  which is more naturally produced from a CLI binary than a
  pass/fail test.
- **MPFR-only first pass, Arb in a follow-up.** Considered;
  rejected because the bead's "65 536 × 5 × 47 = 15.4M" surface
  explicitly covers both backends, and the 12 Arb-primary FnIds
  include the slowest-but-most-stress-relevant Bessel small-arg
  cells (the case where the calibration discipline matters most).

## Cross-ties

- **ADR-0039** — Phase 1g verification architecture closure.
  Defines the per-kernel `error_guard` calibration table and the
  pf-tqzz assertion form; this ADR is the per-release-gate
  execution that ADR-0039 deferred.
- **ADR-0034** — Oracle layer (MPFR + Arb dual-backend).
- **ADR-0035** — Oracle worker protocol; the MIDPOINT verb the
  sweep uses for Arb-primary FnIds.
- **ADR-0038** — Five-mode kernel completeness; the sweep covers
  every IEEE 754-2019 rounding mode per ADR-0038.
- **ADR-0042 / ADR-0043 / ADR-0047 / ADR-0048** — Phase 2 perf
  changes the sweep validates remained calibration-compatible.
- **`feedback_slice_landing_workflow`** — YubiKey-boundary protocol
  followed for this slice (three prompt points: pre-EC2-launch,
  pre-signed-merge, pre-push).
- **`feedback_precision_gated_verification_surface`** — the 2b.1
  rejection lesson; informed the sweep's "validate at the
  precisions the calibrated constants apply at" methodology.
- **ferrodec `project_fd_ykr_campaign_closed`** — the ferrodec
  Decimal32 brute-force campaign whose one-unit-per-instance
  EC2 pattern this sweep adopts.
