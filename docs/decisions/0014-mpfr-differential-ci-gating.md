# ADR-0014: MPFR differential CI gating and implementation choice

- **Status**: accepted (Phase 5 complete)
- **Date**: 2026-05-10

## Status update (slice ci-green, 2026-05-29): the lane runs `--release`

The `differential` CI job had been red on every push since the Phase 1f merge
(2026-05-26), each run hitting GitHub's 6-hour default job timeout. The cause was
neither a hang nor an inherent cost explosion: the job ran
`cargo test --features=differential-mpfr --test 'differential_*'` in the **debug**
profile. pfloat's differential kernels are bignum limb loops the optimiser
collapses; in debug, with overflow checks live and nothing inlined, each
multi-limb operation runs roughly an order of magnitude slower. The slice-6h note
below recorded the full sweep at "roughly 5 minutes total" precisely because it
was measured under `--release`; the CI job had simply never carried the flag.

Two changes after slice 6h compounded the debug cost until it crossed 6 hours:

- Phase 1f (ADR-0038) widened many kernels from NearestEven-only to all five
  `BIT_EXACT_ROUNDING_MODES`, an up-to-5× input multiplier on those ops.
- `TRANSCENDENTAL_PRECISIONS` regained `1024` (the slice 7b AGM-constant
  follow-up), the per-evaluation cost driver.

Both are coverage worth keeping, so neither is the fix.

### Decision

Run the lane with `--release` and cap the job with `timeout-minutes`. No input,
mode, or precision coverage is reduced; the work is run optimised rather than run
less of it (the frugal reading of CLAUDE.md: every cycle should earn its keep,
and an unoptimised bignum sweep does not).

The plan that opened this slice anticipated reducing the CI sweep size. The
measurement redirected that: the missing `--release` flag was the whole problem,
so sweep reduction was unnecessary and would have sacrificed signal for nothing.

### Measurement (receipts)

Apple-silicon dev box, under heavy concurrent load (a parallel ferrodec
all-features run plus other cargo suites pinning cores), so these are loose
upper bounds; a clean 4-vCPU runner has no such contention:

- Debug: the per-test sweep could not clear `agm`/`ai` inside a 240 s/test cap;
  the CI job had been timing out at the 6-hour ceiling.
- Release: the full 35-file suite completed ~34 of 35 tests in ~16 minutes and
  was dominated by `differential_zeta::zeta_matches_mpfr`, a single test that ran
  past 12 minutes on its own. zeta is a single non-parallelizable test function
  whose `s < 0` points at `p = 1024` compose Γ + sin + pow + Borwein; under the
  five-mode tier it is the suite's wall-clock floor.

The first green run on the real GitHub runner (run 26655573667) confirmed the
profile: the differential job took 4549 s (~76 min), of which
`differential_zeta::zeta_matches_mpfr` alone was 3184 s (~53 min) — roughly 70%
of the job in a single non-parallelizable test (the next-largest test was 124 s).
The lane passed, but at 84% of the cap, and zeta's cost is compute, not build, so
a warm cache would not have relieved it.

### Follow-up applied: zeta NearestEven-only at p = 1024 (slice ci-green-zeta)

`tests/differential_zeta.rs` now runs all five modes at p ≤ 256 and
NearestEven-only at p = 1024, restoring zeta's original user-confirmed posture
(the Phase 1f widening, ADR-0038, had added the directed modes at p = 1024). The
directed modes at p = 1024 mostly re-exercise the final-round converter already
covered by the five-mode tier at p ≤ 256, so this is a ~5× cut on the dominant
cell with no loss of meaningful signal. It brings the differential job to roughly
40 minutes, comfortably under the 90-minute tripwire.

`timeout-minutes: 90` stays a regression tripwire, not a budget: a future breach
is a signal to investigate the responsible kernel, not to raise the cap.

### Consequences

- The lane returns to minutes-scale and CI goes green, unblocking the v1.0
  ceremony's dependence on a trustworthy `main`.
- Release disables debug overflow checks. This does not affect the differential
  comparison: pfloat's kernels produce bit-identical results in both profiles
  (a user's release build must, or the library would be broken), and the lane
  asserts bit-exact agreement with MPFR either way.
- Cold release builds of the test binaries cost more compile time than debug,
  amortised by `Swatinem/rust-cache`; the run-time win dwarfs it.

## Status update (slice 7c, 2026-05-16)

Limitation #3 below is closed. Slice 7c (ADR-0022) routes `pow`
through a Ziv interval-test driver plus a square-and-multiply
integer-exponent fast path, so `pow` is correctly rounded under
every IEEE rounding mode (subject to the documented Ziv iteration
cap). `tests/differential_pow.rs` now asserts bit-exact equality
against MPFR across all five rounding modes, replacing the 2 ULP /
NearestEven-only posture, and `pow` is the first transcendental off
the NearestEven-only differential tier. Limitations #1 and #2 were
already resolved by slices 7a (ADR-0016, the bit-exact converter)
and 7b (ADR-0017, AGM-based constants); `TRANSCENDENTAL_PRECISIONS`
now includes 1024. A latent `mpfr_round_of` mapping that sends
NearestAway to MPFR's directed round-away (never exercised, since
every other lane is NearestEven only) was found while tightening the
`pow` lane and is tracked as separate cleanup.

## Status update (slice 6h, 2026-05-10)

First full local MPFR sweep surfaced three structural limitations
on the differential lane that the original ADR text did not
anticipate. All three are documented in the test code and
either patched (constants in `tests/differential/mod.rs`) or
tracked as Phase 6 / Phase 7 follow-ups:

1. **Rounding mode coverage is NearestEven only**, not all five.
   The `bigfloat_to_rug` helper goes via `BigFloat::Display` and
   `rug::Float::parse`; both default to NearestEven. Values
   produced under non-NearestEven rounding lose up to 1 ULP
   through the conversion. Concrete divergences this masks:
   `div(-966132233652331, 1233101814760529)` at p=53/NearestAway,
   `sqrt(2473446)` at p=53/NearestAway, several `fma` cases.
   pfloat under NearestEven matches MPFR bit-for-bit on every
   such case; pfloat under non-NE *also* matches MPFR but the
   conversion can't witness it. Fix is a bit-exact converter
   built on a future `pub raw_parts()` accessor on `BigFloat`
   or a hex/binary radix Display.

2. **Transcendental precisions are capped at 256 bits** via the
   new `TRANSCENDENTAL_PRECISIONS` constant. pfloat's elementary
   transcendentals (exp, ln, pow, sin, cos, tan, atan2, sinh,
   cosh, asinh, erf) all reach into hardcoded 1024-bit
   constants (`ln(2)`, `2/π`, `2/sqrt(π)`, etc.) for argument
   reduction or the leading coefficient. With a 64-bit guard,
   target precisions above ~960 bits exceed those constants'
   reach. At p=1024 the bit-exact MPFR agreement breaks for
   reasons unrelated to the arithmetic — the table simply runs
   out of bits. Phase 7 work: extend constants to 4096 bits or
   compute them on the fly.

3. **`pow` uses 2 ULP tolerance**, not bit-exact. pfloat's slice
   3c pow ships `exp(y · ln(x))` at working precision, which
   accumulates rounding from ln, the multiplication, and exp.
   MPFR's `mpfr_pow` has a fast path for integer exponents that
   avoids the exp/ln composition entirely. The 1 ULP difference
   on, e.g., `pow(63, 9)` at p=53 is the design gap, not a bug.
   The follow-up is either a Ziv-strategy retry pass or an
   integer-exponent fast path in pfloat's `pow`.

**Sweep results** (10⁴ inputs × NearestEven × 3 or 4 precisions
per op, depending on whether transcendentals are involved):
all 22 differential test files pass on macOS local under
`cargo test --features=differential-mpfr,fixed,ops --release`
in roughly 5 minutes total. 768 cargo tests pass overall
(744 existing plus 24 differential test functions).

## Status update (slice 6g, 2026-05-10)

Phase 5 closed with 22 differential test files under `tests/`,
one per op (or per cluster where the underlying kernel is
shared): `differential_{add, sub, mul, div, sqrt, fma, exp, ln,
pow, sin, cos, tan, atan2, sinh, cosh, asinh, erf, gamma,
lgamma, digamma, beta, parse}.rs`. Each is gated by
`required-features = ["differential-mpfr"]` so default builds
stay pure Rust. The `rug` safe wrapper covered the full surface
without an `unsafe` block in test code.

The 10⁴-in-CI / 10⁶-local split via `PFLOAT_DEEP=1` is in place;
the actual deep-sweep ritual is the Phase 5 exit pass and
becomes a regular slice-cadence step once Phase 6 (tier-2
specials) lands.

The macOS-runner cost of building `gmp-mpfr-sys` from source is
real (~2 minutes cold). pfloat's CI test matrix avoids this by
not enabling `differential-mpfr` in the cross-OS matrix row —
the dedicated `differential` job (Linux-only) covers it with
`apt-get install libmpfr-dev libgmp-dev` in ~30 seconds. The
gating is documented in `.github/workflows/ci.yml`.

The astro-float-removal clause from this ADR's body remains in
effect; ADR-0008's secondary-oracle text is permanently
superseded.

## Context

ADR-0008 established MPFR via `gmp-mpfr-sys` as the primary
differential oracle, with `astro-float` as a secondary pure-Rust
oracle on the default Linux lane. Phase 5 planning revisited the
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

A second question Phase 5 settles: the CI sweep size. DESIGN.md
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
