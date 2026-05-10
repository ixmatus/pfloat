# ADR-0008: Differential testing oracle, `gmp-mpfr-sys` on a feature-gated CI lane

- **Status**: proposed
- **Date**: 2026-05-10

## Context

Trust evidence for an arbitrary-precision arithmetic library has
to come from somewhere. The candidates:

- **MPFR**, via `gmp-mpfr-sys`. The reference implementation; the
  one every other library (Python's `mpmath`, Mathematica's
  arbitrary-precision module, GNU Octave's high-precision package)
  is calibrated against. Authority is unambiguous. Cost: pulls in
  a C build dependency at test time.
- **astro-float**, the closest pure-Rust peer. Authority is
  weaker (a peer can have the same bug pfloat does), but the
  build is clean.
- **Both**, treated as independent oracles. A bug in pfloat is
  flagged when it disagrees with either; a bug in MPFR or
  astro-float surfaces as a cross-oracle disagreement that pfloat
  is uniquely positioned to publish.
- **Self-test only**, via internal property checks (round-trip
  identities, monotonicity, symmetry). Necessary but insufficient;
  a uniformly-wrong implementation passes self-tests.

The C build dependency is real but limited. `gmp-mpfr-sys` builds
GMP and MPFR from source on Linux and macOS using the system
toolchain. Windows requires more work. The dependency is dev-only:
production users of pfloat never see it.

The principles in `~/.claude/CLAUDE.md` and the project CLAUDE.md
forbid FFI as a runtime dependency for performance reasons. They
say nothing against FFI as a dev-dependency for trust evidence.
The split is principled: production builds stay pure Rust; testing
uses the most authoritative oracle available.

## Decision

`gmp-mpfr-sys` is the primary differential oracle. It runs as a
feature-gated dev-dependency:

```toml
[dependencies]
# (none from gmp-mpfr-sys; it is dev-only)

[dev-dependencies]
gmp-mpfr-sys = { version = "1", optional = true }

[features]
differential-mpfr = ["dep:gmp-mpfr-sys"]
```

CI runs the differential tests on a Linux-only lane. The default
lanes (Linux, macOS, embedded cross-compile) build without the
feature and stay pure Rust.

`astro-float` ships as a secondary dev-dependency on the default
Linux lane. It is faster to run, builds without a C toolchain, and
catches a class of "we both have the same f64-style bug" issue
that MPFR alone would miss. Its conformance is partial relative to
MPFR; the secondary status is appropriate.

For each operation and each function, the differential test:

1. Generates random inputs across the precision and value range
   relevant to the operation.
2. Computes the result in pfloat.
3. Computes the result in MPFR (or astro-float) at the matching
   precision and rounding mode.
4. Compares bit patterns of the mantissa and exponent at the
   target precision.
5. Compares exception flag state.

A divergence either fails the test (CI red) or is documented in
`docs/decisions/known-divergences.md` as a deliberate spec
interpretation choice. The latter is rare; pfloat is a literal
reading of IEEE 754-2019 with MPFR's conventions for ambiguous
clauses.

## Consequences

**Wins:**

- The differential lane gives the strongest possible authority for
  pfloat's correctness. Every operation pfloat ships is bit-equal
  to MPFR at every tested precision and rounding mode.
- The default lane stays pure Rust, preserving the embedded story
  and the WebAssembly story for users who run their tests on
  those targets.
- A cross-oracle disagreement (pfloat vs MPFR vs astro-float)
  isolates the bug to exactly one of the three. Publishable
  finding either way.

**Costs:**

- One extra CI runner. GitHub Actions Linux runners are cheap; not
  blocking.
- The differential lane requires a `libgmp` and `libmpfr` install,
  done by `gmp-mpfr-sys`'s build script. Slightly slower CI build
  step (~minutes); a one-time cost per cache invalidation.
- Windows differential testing is not available in 1.0. Documented
  as a CI-coverage caveat in the README; pfloat itself runs on
  Windows via the pure-Rust default lane.
- Translation between pfloat's flag enum and MPFR's flag bits lives
  in a small module, exercised by every differential test. The
  translation is mechanical but adds a small surface to maintain.

## Related

- ADR-0007 (rounding mode and flags, both compared against MPFR)
- ADR-0009 (verification scaffolding, of which this lane is part)
- DESIGN.md, "Verification" section, "Differential testing"
  subsection.
