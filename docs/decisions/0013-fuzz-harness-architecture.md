# ADR-0013: Fuzz harness architecture

- **Status**: accepted (Phase 6 complete)
- **Date**: 2026-05-10

## Status update (slice 6g, 2026-05-10)

Phase 6 closed with seven fuzz targets covering the full op
surface: `parse`, `arith`, `exp_log_family`, `trig`, `hyperbolic`,
`specials`, and `fmt`. The cluster-level grain (one target per
feature cluster, not one per op) kept the count bounded and
manageable. The non-workspace `fuzz/` subcrate layout worked
without friction; the `path = ".."` dep on pfloat compiled
cleanly under every slice.

The OSS-Fuzz upstream PR scaffolding lives in `fuzz/oss-fuzz/`
(Dockerfile, build.sh, project.yaml, README.md). Submission is a
post-Phase-6 follow-up. The upstream PR is tracked in the slice
6g merge commit.

The no-checked-in-corpus policy held. Counterexamples that
surface in CI's 60-second smoke or in local deep runs are
expected to be promoted to `.proptest-regressions` entries
against the relevant property test, per
`feedback_proptest_regression_seeds.md`. No counterexamples
surfaced during slice cadence; the property suite stayed clean.

## Context

DESIGN.md scopes Phase 6 to include `cargo-fuzz` harnesses and an
eventual OSS-Fuzz integration. ferrodec ships six libfuzzer-sys
targets under `fuzz/fuzz_targets/`, with no checked-in corpus, and
treats fuzz as a panic-freedom + identity-invariant smoke gate rather
than an oracle-driven correctness gate (accuracy is the differential
oracle's job per ADR-0008 / ADR-0014).

pfloat's surface is bigger than ferrodec's. Phase 1–4 has shipped
~38 operations; fuzz targets organized per-op would balloon the count.
The right grain is the cluster: one target per feature cluster
(`parse`, `arith`, `exp_log_family`, `trig`, `hyperbolic`, `specials`,
`total_cmp`, `fmt`).

Three design questions matter:

1. **Where does fuzz live?** ferrodec ships a sibling cargo subcrate
   at `fuzz/`, with its own `Cargo.toml` and a `path = ".."` dep on
   the main crate. pfloat copies the pattern.
2. **What does each target check?** Panic-freedom is mandatory.
   Cheap identity invariants (`a + 0 ≡ a`, `parse(fmt(x)) == x` at
   matching precision) are good additions. Oracle-based accuracy
   checks belong to the MPFR differential lane, not fuzz.
3. **Checked-in corpus or libFuzzer-evolved?** A checked-in corpus
   commits a regression set with the repo. libFuzzer evolves its own
   corpus per run; counterexamples that survive get promoted to
   `.proptest-regressions` via the existing seed-commit convention.
   The latter has lower maintenance overhead; ferrodec uses it.

## Decision

`fuzz/` is a sibling cargo subcrate, not a workspace member of pfloat.
Layout:

```
fuzz/
├── Cargo.toml                # publish = false, libfuzzer-sys + arbitrary deps,
│                             # path = ".." dep on pfloat with all-features.
├── .gitignore                # target, corpus, artifacts
└── fuzz_targets/
    ├── parse.rs              # 6a: canonical target — parse panic-freedom.
    ├── arith.rs              # 6b: add/sub/mul/div/sqrt/fma identity invariants.
    ├── exp_log_family.rs     # 6c
    ├── trig.rs               # 6d
    ├── hyperbolic.rs         # 6d
    ├── specials.rs           # 6e
    ├── total_cmp.rs          # 6f
    └── fmt.rs                # 6f
```

Each target body follows the libfuzzer-sys `fuzz_target!` macro
shape:

```rust
#![no_main]
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if let Ok(s) = core::str::from_utf8(data) {
        let _ = pfloat::BigFloat::parse_str(s, /* prec */ 113, RoundingMode::NearestEven);
    }
});
```

Targets that exercise arithmetic (`arith`) use `arbitrary::Arbitrary`
derives on a small input struct that bundles `(op_tag, a_bits,
b_bits)`, with the operation dispatched on the tag. Identity
invariants (`a + 0 ≡ a`, etc.) appear after the dispatch, only
asserted when both inputs are finite.

**No corpus is checked into the repo.** `fuzz/corpus/` is gitignored.
libFuzzer evolves its working corpus across runs; meaningful
counterexamples get extracted as proptest regression seeds (the
`feedback_proptest_regression_seeds.md` workflow) so the regression
record lives with the proptests rather than in the fuzz corpus.

CI invocation (new job in `.github/workflows/ci.yml`):

```yaml
fuzz-smoke:
  name: fuzz smoke (60s per target)
  runs-on: ubuntu-latest
  steps:
    - uses: actions/checkout@v4
    - uses: dtolnay/rust-toolchain@master
      with:
        toolchain: ${{ env.PFLOAT_NIGHTLY }}
    - run: cargo install cargo-fuzz
    - run: cargo fuzz build
    - run: cargo fuzz run parse -- -max_total_time=60
    # Slices 6b–6f add additional 60s smoke runs per target.
```

The smoke run is blocking. 60 seconds per target is enough to surface
new panics introduced by a slice while keeping CI minutes bounded.
Deep coverage runs locally: `cargo fuzz run <target> -- -max_total_time=3600`.

OSS-Fuzz upstream PR is **deferred to slice 6g**. A
`fuzz/oss-fuzz/` directory with `Dockerfile`, `build.sh`, and
`project.yaml` lands then. The PR submission target is
`google/oss-fuzz/projects/pfloat/`.

## Consequences

**Wins:**

- Panic-freedom across the surface is detectable in CI without a
  Kani-grade infrastructure burden. libFuzzer's adaptive corpus
  finds bugs at scales proptest cannot reach.
- Counterexamples that matter end up as proptest regressions, where
  they get versioned with the code. The fuzz directory stays empty.
- Cluster-level targets keep the count bounded as pfloat grows
  through Phase 5 (Bessel, zeta, Airy) and beyond.
- The OSS-Fuzz upstream PR gives pfloat continuous fuzzing without
  pfloat paying the compute bill.

**Costs:**

- 60-second smoke is shallow. Deep coverage requires local manual
  runs. Acceptable because the MPFR differential lane handles
  accuracy at depth; fuzz is the panic-freedom gate.
- `cargo install cargo-fuzz` in CI adds ~30s to cold runs. Mitigated
  by `Swatinem/rust-cache` caching the installed binary.
- The non-workspace `fuzz/` subcrate produces a slightly unusual
  cargo experience: `cargo build` from pfloat/ does not build the
  fuzz subcrate, which is what we want, but a contributor unfamiliar
  with cargo-fuzz might initially miss it. Documented in
  `fuzz/Cargo.toml`'s top-of-file comment.

## Trigger to revisit

- OSS-Fuzz onboarding produces feedback (corpus growth rates, bug
  classes) that suggests a different harness layout would catch
  more. Update this ADR to reflect.
- Phase 5 / Phase 7 surface introduces input shapes (e.g. complex
  numbers, intervals) that the current cluster layout does not fit
  well. New target file per shape; update this ADR.
- libFuzzer is replaced by a different fuzzing engine
  (Honggfuzz, AFL++) for a specific class of bug; would warrant
  a follow-up ADR.

## Related

- ADR-0009 (verification scaffolding)
- ADR-0008, ADR-0014 (MPFR differential — the accuracy oracle that
  fuzz delegates to)
- DESIGN.md, "Verification" section, "Fuzzing" subsection
- ferrodec/fuzz/ (the template for harness shape and corpus
  policy)
- `feedback_proptest_regression_seeds.md` (the bridge from fuzz
  counterexamples to versioned regression seeds)
