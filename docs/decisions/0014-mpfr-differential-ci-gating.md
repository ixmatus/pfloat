# ADR-0014: MPFR differential CI gating and implementation choice

- **Status**: accepted
- **Date**: 2026-05-10

## Context

ADR-0008 established MPFR via `gmp-mpfr-sys` as the primary
differential oracle, with `astro-float` as a secondary pure-Rust
oracle on the default Linux lane. Phase 6 planning revisited the
secondary-oracle clause and concluded that **MPFR alone is the right
sole oracle** for pfloat:

- One oracle keeps divergence triage single-source. A bug in pfloat
  that disagrees with two oracles produces three-way debugging; a
  bug that disagrees with one produces two-way debugging.
- `astro-float` is the closest pure-Rust peer, but it has not yet
  hit the maturity (sustained version cadence, downstream production
  use, formal-verification adjacency) that would make it a useful
  independent witness against MPFR. Until that gap closes, a
  disagreement between pfloat and astro-float is almost always
  pfloat or astro-float, not MPFR.
- ferrodec uses astro-float because ferrodec is decimal (radix 10);
  MPFR's decimal coverage is partial. pfloat is binary; MPFR is the
  canonical authority. The radix difference makes ferrodec's choice
  inapplicable to pfloat.

ADR-0008's literal text declared `gmp-mpfr-sys = { version = "1",
optional = true }` as the dev-dependency. `gmp-mpfr-sys` is a
low-level raw-bindings crate; it requires `unsafe` to call. pfloat
sets `unsafe_code = "forbid"` at the crate root, which applies to
pfloat's own code but not to dev-dependencies. Still, the natural
safe wrapper for `gmp-mpfr-sys` is **`rug`**, which provides a
high-level type-safe API over the same MPFR binary. `rug` is the
right surface for pfloat's differential tests.

A second question Phase 6 settles: the CI sweep size. DESIGN.md
states "10⁶ random inputs per op per rounding mode" as the Phase 1
exit criterion for arithmetic. Running 10⁶ inputs per
op × precision × rounding mode × cluster in every CI run is
expensive. The natural split is **10⁴ in CI, 10⁶ local**, gated
behind a `PFLOAT_DEEP=1` environment variable. CI catches
regressions cheaply; the deep sweep runs before each Phase merge.

## Decision

This ADR refines ADR-0008 on three points and stands on its own on
two. ADR-0008's "secondary astro-float oracle" clause is **superseded
by this ADR**; the remainder of ADR-0008 continues to apply.

### Implementation: `rug` (which wraps `gmp-mpfr-sys`)

`Cargo.toml`:

```toml
[features]
differential-mpfr = []

[target.'cfg(unix)'.dev-dependencies]
rug = { version = "1", default-features = false, features = ["float"] }

[[test]]
name = "differential_add"
path = "tests/differential_add.rs"
required-features = ["differential-mpfr"]
```

`rug` is Unix-only because `gmp-mpfr-sys`'s C build does not produce
a working binary on stock Windows toolchains. Windows differential
coverage is not available; users rely on cross-validation via Linux
CI.

Cargo does not allow dev-dependencies to be `optional`, so the
feature itself activates no dependency. The differential test
binaries declare `required-features = ["differential-mpfr"]` in
their `[[test]]` entries; `cargo test` without the feature skips
them entirely, and the test binaries do not link `rug`. The `rug`
crate itself still compiles on every Unix `cargo test` invocation,
but the C compile happens once per cache lifetime; warm builds are
seconds.

### Test layout

```
tests/
├── differential/
│   ├── mod.rs              # Shared helpers:
│   │                       # - bigfloat_to_rug(&BigFloat) -> rug::Float
│   │                       # - rug_to_bigfloat(&rug::Float, prec) -> BigFloat
│   │                       # - mpfr_round_of(RoundingMode) -> rug::float::Round
│   │                       # - deterministic ChaCha8Rng generators
│   └── ...
├── differential_add.rs     # Canonical example. 10⁴ inputs / 4 precisions
│                           # / 5 rounding modes; PFLOAT_DEEP=1 escalates to 10⁶.
└── (further differential_*.rs files added per slice 6b–6f)
```

Each `differential_<op>.rs` file is gated by
`#![cfg(all(feature = "differential-mpfr", unix))]`.

### Comparison strategy

Convert pfloat's result back through `rug::Float` and compare in
MPFR's normalized form:

```rust
let pfloat_as_rug = bigfloat_to_rug(&pfloat_r);
let mpfr_r = rug::Float::with_val_round(p, &mpfr_a + &mpfr_b,
                                       mpfr_round_of(rm)).0;
assert_eq!(pfloat_as_rug, mpfr_r);
```

Binary radix means one canonical normalized form per finite value —
no cohorts. Bit-for-bit equality is the right test.

### Determinism

All input generation runs through a `ChaCha8Rng` seeded from a
`u64` that `proptest` produces. Counterexamples are reproducible
from the seed; failing seeds get committed as
`.proptest-regressions` entries (the existing convention from
`feedback_proptest_regression_seeds.md`).

### CI gating

New job in `.github/workflows/ci.yml`:

```yaml
differential:
  name: MPFR differential (linux)
  runs-on: ubuntu-latest
  steps:
    - uses: actions/checkout@v4
    - uses: dtolnay/rust-toolchain@master
      with:
        toolchain: ${{ env.PFLOAT_NIGHTLY }}
    - uses: Swatinem/rust-cache@v2
    - run: sudo apt-get update && sudo apt-get install -y libmpfr-dev libgmp-dev
    - run: cargo test --features=differential-mpfr --test 'differential_*'
```

Blocking. macOS and Windows runners skip this job; the matrix `os`
key does not include them for this lane.

### Sweep size

- **CI default**: 10⁴ inputs per (op × precision × rounding mode).
- **Local deep (`PFLOAT_DEEP=1`)**: 10⁶ inputs per
  (op × precision × rounding mode), run before each Phase merge.
- Precisions exercised in CI: 53, 113, 256, 1024 bits.
- All five rounding modes per op.

## Consequences

**Wins:**

- Bit-for-bit MPFR agreement is the strongest possible authority
  for pfloat's correctness. Every op pfloat ships matches MPFR at
  every tested precision and rounding mode.
- Single-oracle triage: a divergence is unambiguously a pfloat bug
  (or a documented spec interpretation choice).
- `rug`'s safe API keeps the test code clean. No `unsafe` blocks in
  pfloat's test directory; `gmp-mpfr-sys`'s bindings stay hidden
  beneath `rug`.
- The 10⁴ / 10⁶ split balances CI minutes against coverage. CI
  catches regressions cheaply; the per-Phase deep sweep gives the
  10⁶ coverage that DESIGN.md asks for at exit time.

**Costs:**

- `rug` adds `gmp-mpfr-sys` as a transitive dev-dependency. Cold CI
  builds spend ~2 min compiling GMP + MPFR + MPC from source.
  Mitigated by `Swatinem/rust-cache`; warm builds are seconds.
- Windows differential coverage is not available. Documented in
  README.md as a CI-coverage caveat. Pure-Rust pfloat itself runs
  on Windows; only the differential lane does not.
- The deep sweep is a manual ritual. The slice cadence
  (`feedback_pfloat_slice_cadence.md`) gets a Phase-exit step
  in the next revision: run `PFLOAT_DEEP=1 cargo test --features=
  differential-mpfr --test 'differential_*'` before signing the
  Phase merge.
- This ADR drops astro-float as a secondary oracle. We lose the
  cross-oracle disagreement detector (pfloat vs MPFR vs
  astro-float). The trade is acceptable: MPFR is canonical, and a
  pfloat-vs-MPFR disagreement is already a publishable finding.

## Trigger to revisit

- `gmp-mpfr-sys` gains a Windows toolchain that works in stock
  GitHub Actions runners. At that point the Windows lane can join
  the differential matrix.
- `astro-float` reaches the maturity threshold (1.0 release with
  formal verification adjacency or sustained downstream production
  use) that would make it a useful independent witness. Re-introduce
  it as a secondary oracle.
- A faster differential oracle (LLM-based symbolic execution, formal
  verification via Creusot proving correctness directly) makes the
  10⁴/10⁶ sweep redundant for some ops.

## Related

- ADR-0008 (refined by this ADR; the secondary-astro-float clause
  is superseded here)
- ADR-0009 (verification scaffolding)
- ADR-0012 (Kani harness architecture)
- ADR-0013 (fuzz harness architecture)
- DESIGN.md, "Verification" section, "Differential testing"
  subsection
- `feedback_proptest_regression_seeds.md` (regression seed
  convention)
- `feedback_pfloat_slice_cadence.md` (slice workflow that gets the
  per-Phase deep sweep step)
