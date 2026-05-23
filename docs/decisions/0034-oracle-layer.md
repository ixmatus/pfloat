# ADR-0034: Oracle layer for the Phase 1 exhaustive `f32` sweep

- **Status**: accepted
- **Date**: 2026-05-23

## Context

Slice p1.2 closed the five-finding has-errors class on pfloat's
v1.0 surface via the Ziv-retry kernel upgrade across the exp/log
and specials family. The Lefèvre-Muller worst-case corpus
(twenty-four functions / 1200 cases at p=53 NE) verifies the
sampled hard-to-round inputs upstream identified. It does not
verify correct rounding across the full binary64 input space, and
its coverage is sampled rather than exhaustive even within the
hard-to-round region.

The Phase 1 plan's exhaustive `f32` sweep is the rigor escalation
that closes that gap. ADR-0033 already records the sequencing
decision: the sweep runs to completion before the v1.0 tag rather
than after, because shipping a v1.0 that claims correctly-rounded
transcendentals while carrying latent rounding defects would
permanently anchor the credibility of an immutable crates.io
version against a less-rigorous-than-stated standard.

The sweep needs an oracle. "Matches MPFR" is the existing
differential-lane statement and is necessary but not sufficient
for the v1.0 claim: MPFR's correct-rounding guarantee covers the
basic operations and a documented subset of functions; for tier-2
specials MPFR is itself the reference implementation and is
believed correctly rounded in practice without shipping the formal
guarantee. For modified Bessel I and K, MPFR has no primitive at
all. The Phase 1 plan addresses both gaps by layering two
independent oracle backends behind one interface:

- MPFR (in process, via `rug`) for the standard elementary and
  most specials surfaces.
- Arb (via `python-flint` in a long-lived subprocess) for the
  Bessel I/K family and Airy, where MPFR has no primitive.

The oracle layer's correctness lift over the existing differential
lane is the **enclosure** posture: the oracle returns a proven
bracket `[lo, hi]` of the true value, not a rounded scalar. The
verifier then asks whether both bracket endpoints round to the
same `f32` under the caller's mode; if they do, every point in the
bracket including the true value rounds there too and the
correctly-rounded `f32` is determined. If they straddle a rounding
boundary, the oracle's working precision doubles and the bracket
tightens; this is Ziv-at-oracle, the same interval test pfloat's
own `ziv_round` driver (slice p1.2, ADR-0022) uses internally.

The lift over "evaluate at one mode, compare" is the elimination
of the "did the oracle round the same way" question. The verifier
never compares pfloat to "the oracle's rounded answer." It
compares pfloat to the unique `f32` inside the oracle's certified
enclosure. The trait returns a proven bracket, not a value.

## Decision

**Enclosure type and OracleBackend trait.** A bracket-shaped
return type and a small trait wrap any backend that can produce
proven enclosures of a function value at a requested working
precision:

```rust
pub struct Enclosure {
    pub lo: OracleReal,
    pub hi: OracleReal,
}

pub trait OracleBackend {
    fn enclose(&self, f: FnId, input: F32Bits, working_prec: u32)
        -> Enclosure;
    fn name(&self) -> &'static str;
}
```

`OracleReal` is the backend's internal high-precision real type
(an `mpfr::Float` for the MPFR backend, an `arb` for the Arb
backend); the verifier never inspects it beyond the
`round_to_f32(&OracleReal, mode)` call that asks "what `f32` does
this real value round to under `mode`." The trait is `Send` so
shards can fan out across cores.

**FnId enum, single global dispatch.** A single
`enum FnId { Exp, Ln, Sin, ..., I0, I1, In(i32), ... }` identifies
each function the sweep verifies. The MPFR backend's `enclose`
matches on `FnId` and dispatches to the corresponding rug call;
the Arb backend matches on `FnId` and dispatches to a
`python-flint` worker request. The alternative (per-function
typestate trait) is more idiomatic Rust but produces boilerplate
the harness side would have to bridge through the trait object;
the enum is simpler at the harness level and is the convention the
Phase 1 plan's verify_input signature already uses.

**MPFR backend produces brackets via directed rounding.** For
each `enclose(f, x, w)` call the backend evaluates `f` at working
precision `w` twice through `rug`: once with MPFR_RNDD
(toward-negative-infinity), once with MPFR_RNDU
(toward-positive-infinity). The two values provably bracket the
true result; the bracket tightens as `w` grows. The cost is two
MPFR calls per oracle invocation rather than one, which is
acceptable: the existing differential lane already runs `rug`
single-mode evaluations at every input; the oracle lane shifts to
two-mode evaluations as the price for the enclosure guarantee.

Special cases. For `f(x) = +inf` (overflow at MPFR's working
precision), both directed roundings overflow; the enclosure is
`[+inf, +inf]` and `certified_round_f32` returns the f32 infinity
when both round there. Similarly for NaN inputs and domain-edge
NaN outputs. Subnormals are handled by reading the directed
rounding results without re-emulating the binary32 subnormal
behavior: MPFR's directed modes produce mathematically correct
brackets that the `round_to_f32` step interprets through the f32
subnormal mapping.

**Arb backend posture (next slice).** The Arb backend is recorded
here for design completeness; its implementation lands in a
follow-up slice. The shape: a long-lived `python-flint` worker
process reads input batches over a pipe and streams back `(lo, hi)`
pairs from the `arb` ball type's midpoint and radius (Arb is ball
arithmetic; the enclosure is constructive in the type). The Rust
harness owns enumeration and the boundary check; the worker only
encloses. LGPL isolation: FLINT and Arb are LGPL; the worker is a
test-time only subprocess invoked via the `differential-mpfr`
feature lane (no in-process link, no shipped-crate dependency).

**Verification core: Ziv-at-oracle.** `certified_round_f32(enc,
mode)` returns `Some(f32)` when both endpoints round to the same
value under `mode`, `None` otherwise. `verify_input(oracle, f, x,
mode)` runs a Ziv loop at the oracle: start working precision at
`START_PREC = 64` (comfortably above f32's 24 bits), double on
`None`, cap at `MAX_PREC = 1024`. An `OracleInconclusive` verdict
at the cap is rare, is not a pfloat failure, and goes to a
separate worst-case-candidate file rather than the failure corpus.

**Tiered sweep budget.** Per the Phase 1 plan and the slice p1.3
scope decision: cheap MPFR-primary functions get literal 2^32
enumeration (over many hours in the release CI lane); expensive
ones run a dense sample of ~2^20 inputs supplemented with all
boundary inputs (NaN, ±inf, ±0, subnormal min/max, normal min/max,
power-of-two boundaries). The status table records `exhaustive`
vs `sampled(N)` per function. Never claim exhaustive where
sampled.

**Two run modes.** A per-push smoke gate runs ~2^10 inputs per
MPFR-primary function under all five rounding modes on every
differential-mpfr CI build; asserts zero mismatches, zero panics,
zero oracle-inconclusive. A standalone runner
(`examples/oracle_sweep.rs`) takes `--function`, `--exhaustive` /
`--sample N`, `--mode`, `--output` flags and runs the full sweep
when invoked manually or on release CI. The runner is the
artifact a future maintainer would re-execute; the smoke gate is
the per-commit signal.

**Status table schema (machine-readable).** One row per
(function, rounding-mode set) emitted as TOML, matching the Phase
1 plan's schema: `function`, `order` (for Bessel),
`kernel_kind`, `domain_coverage`, `oracle`,
`oracle_independence`, `rounding_modes`, `rounding_status`
(`correctly-rounded` | `faithful` | `has-errors`), `worst_ulp`,
`mismatch_count`, `inconclusive_count`, `panic_count`, `vectors`
(path to the per-function regression corpus). At p1.exit the
status table publishes in the README and docs; pre-exit it lives
under `tests/oracle/status/<fn>.toml` as the in-tree machine
artifact.

**Regression corpus capture.** Every `Mismatch` and every
`faithful`-not-`correctly-rounded` input has its exact `f32` bits
(and mode) appended to `tests/vectors/<fn>_regression.bin`. Fixes
are then regression tested against the captured inputs, and the
expensive exhaustive sweep does not have to re-discover known hard
cases on every run. `OracleInconclusive` inputs go to
`tests/vectors/<fn>_inconclusive.bin` (worst-case candidates).
`Panic` inputs go to `tests/vectors/<fn>_panic.bin` and the panic
regression file runs on every CI push, not only at the next
per-release sweep.

**Location: `tests/oracle/`.** The oracle and harness modules live
under `tests/oracle/` as test infrastructure (matching the
`tests/differential/mod.rs` precedent for the existing
differential lane). Gated `#[cfg(feature = "differential-mpfr")]`
through the existing CI lane; never linked into the shipped crate.
The runner binary lives at `examples/oracle_sweep.rs` for the same
gating.

**L-M corpus as adversarial seeds.** The sweep prefixes each
function's input iteration with the existing L-M corpus entries
(from `tests/differential/lefevre_muller_data.rs`). These are the
adversarial cases least likely to be hit by random or enumerated
sampling within a finite budget; they earn the function the
strongest correctly-rounded claim the sweep can make on a sampled
budget.

## Consequences

**The v1.0 correctly-rounded claim is substantiated per function
under the table.** A user evaluating "can I depend on this for X"
gets a per-function answer from the published status table without
reading the source.

**The harness scales from per-push smoke to per-release exhaustive
without architectural change.** The same `verify_input` runs
inside both modes; the difference is the input iterator
(`(0..1024).map(F32Bits::from)` for smoke vs
`(0..u32::MAX).map(F32Bits::from)` for exhaustive). The runner
binary handles the iterator selection.

**Two-evaluation oracle cost is the durable price of the enclosure
posture.** MPFR backend doubles per-input runtime compared to the
existing single-mode differential lane. Acceptable for CI; the
exhaustive sweep's CPU budget is per-release rather than
per-push.

**Ziv-at-oracle adds a precision-doubling loop on hard inputs.**
Most inputs (especially mode-uniform overflow/underflow and
non-hard-to-round normals) certify at the first guard. Hard inputs
pay up to MAX_PREC / START_PREC = 16 retries (log2(1024 / 64) = 4
guard doublings). The honest measure-zero caveat: if an input
exhausts MAX_PREC the verdict is `OracleInconclusive` and the
input goes to the worst-case-candidate file, not a pfloat failure.

**Arb backend isolation by subprocess.** Choosing `python-flint`
in a subprocess (rather than a direct Rust FFI link) preserves
pfloat's pure-Rust shipped surface: the LGPL FLINT/Arb libraries
are CI-only, never linked into the shipped crate. The IPC cost
(input batches over a pipe, `(lo, hi)` streamed back) is
acceptable for the Arb-only specials whose enclosure is already
expensive in absolute terms.

**Sharding is not part of slice p1.3.** Single-threaded driver
ships first; the shard coordinator and rebalancer the plan
describes ("a shard coordinator detects this shard has been
running an order of magnitude longer than its peers and
rebalances") is a follow-up slice if release-CI runtime motivates
it. The single-threaded runner is sufficient for the per-release
budget when run on a modest core count without explicit
parallelism.

**Findings are discovery, not failure.** The sweep's job is to
surface latent has-errors that the L-M corpus did not hit. If a
function comes back has-errors at slice p1.3's smoke gate, the
finding is captured as a regression corpus entry and a defect
bead; the kernel fix lands in a follow-up slice unless cheap
enough to absorb in-slice.

## Alternatives considered

**Single-evaluation oracle (no enclosure).** Compare pfloat's
output to MPFR's output under the same mode at working precision
high enough to be "definitive." Rejected: this is what the
existing differential lane does, and it cannot answer "the oracle
might have rounded the same wrong way pfloat did" for functions
where MPFR's correct-rounding guarantee is implementation-belief
rather than published proof. The enclosure posture makes the
oracle's rounding the verifier's question to answer, not assumed
correct.

**Per-function trait dispatch instead of an FnId enum.**
Idiomatic Rust would define `trait OracleFn { type Input; fn
enclose_with(...); }` and one type per function. The harness side
would then have to bridge through a trait-object collection or
runtime dispatch. Rejected for slice p1.3 in favor of the simpler
FnId enum; if the per-function-trait shape proves more
maintainable as the surface grows, a follow-up slice can refactor
without changing the verifier's correctness argument.

**In-process Arb via `flint-sys` FFI.** Rejected on the
permacomputing-horizon discipline: FFI to LGPL libraries pollutes
the shipped crate's dependency tree (even behind a feature gate
the toolchain has to find the system library). The subprocess
posture isolates Arb completely.

**Per-push exhaustive smoke (no separate runner).** Rejected on
runtime grounds: a 2^32 sweep on even one function is hours of
CPU; per-push cannot accommodate it. The split into smoke (per
push) plus runner (per release / manual) is the durable shape.

## References

- Plan: `docs/decisions/plans/phase-1-correctness-sweep.md` (the
  full Phase 1 work breakdown; the Oracle layer section spelled
  out the type signatures this ADR adopts verbatim).
- ADR-0014 — the existing MPFR differential lane gating posture
  (CI-only, behind `differential-mpfr`). The oracle layer extends
  this same posture.
- ADR-0022 — the Ziv interval test for the `pow` kernel and its
  lift to the shared `ziv_round` driver at slice p1.2; the
  Ziv-at-oracle loop here uses the same interval-test correctness
  argument applied to the oracle's bracket rather than to a
  pfloat working-precision evaluation.
- ADR-0032 — libm reciprocal and root kernels stay direct primary
  (the alias-absence rationale recorded in the surface freeze).
- ADR-0033 — the Phase 1 sweep runs to completion before v1.0;
  the sequencing decision that makes this slice's harness a v1.0
  gate.
- `docs/v1.0-surface.md` — the frozen unary surface this oracle
  layer verifies.
