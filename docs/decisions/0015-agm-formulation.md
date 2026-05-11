# ADR-0015: AGM kernel uses Gauss's iteration with an independent `agm` feature flag

Status: accepted (slice 6l shipped)

## Context

Phase 6 opens with the arithmetic-geometric mean (AGM) as the smallest
tier-2 special. AGM is the iterated common limit of the arithmetic and
geometric means, and is a building block for two downstream items: the
Phase 7 plan to compute transcendental constants on-the-fly past the
1024-bit hardcoded tables (Brent–Salamin for `π` and Salamin's variant
for `ln(2)` both layer on top of AGM), and the elliptic-integral
surface that some Bessel asymptotics use post-1.0.

Two questions needed answers before any code landed:

- **Algorithmic form.** Gauss's iteration `(a, b) → ((a + b) / 2,
  sqrt(a · b))` is the canonical AGM. Brent–Salamin's variant tracks
  an auxiliary `c_n = (a_n − b_n) / 2` and a running sum of
  `2^n · c_n^2` so that the same iteration produces both AGM and a
  byproduct that converges to `π`. The byproduct only matters when
  computing `π` specifically; the standalone AGM gains nothing.
- **Feature gating.** pfloat already has `exp-log`, `trig`, and
  `specials` (with `specials = ["trig"]`) as cluster features in
  `Cargo.toml`. AGM uses only the Phase 1 arithmetic primitives
  (`add`, `sub`, `mul`, `div`, `sqrt`), so binding it to `specials`
  would force embedded users to pull in trig and the gamma family
  for an unrelated capability. The alternative is a dedicated `agm`
  feature that depends only on `big`.

## Decision

1. Implement the Gauss iteration directly. The kernel computes at
   working precision `target_precision + 64` and stops once the
   binary exponent of `|a_n − b_n|` falls below `−p_work − 4`. The
   iteration cap is `64`, which covers any precision pfloat admits
   (quadratic convergence doubles the bit-agreement count each step,
   so the cap is a safety belt, not a tuned bound). The
   Brent–Salamin auxiliary is omitted; when slice 7b lands the
   on-the-fly `π` computation, the additional bookkeeping lives in
   that kernel, not in `agm`.

2. Add a new top-level Cargo feature `agm = ["big"]`. The `mod math`
   declaration in `src/lib.rs` widens its gate to
   `any(feature = "exp-log", feature = "agm")` so the `math` parent
   module compiles when `agm` is enabled standalone. Inside
   `src/math/mod.rs`, `mod agm` is gated by `feature = "agm"`. The
   existing `specials = ["trig"]` chain stays unchanged.
   `differential-mpfr` gains `agm` so the new
   `tests/differential_agm.rs` test crate builds against the same
   feature set as the rest of the differential lane.

## Consequences

- An embedded caller wanting AGM without the rest of the math
  surface compiles with `--no-default-features --features=big,agm`.
  The compiled artifact contains no exp/log, no trig, no tier-1
  specials.
- Slice 7b (on-the-fly transcendental constants) will depend on
  both `agm` and `exp-log`; its feature gate becomes
  `any(all(exp-log, agm), exp-log)` or, more cleanly,
  `exp-log = ["big", "agm"]`. The dependency widens cleanly when
  the time comes.
- Future tier-2 work that uses AGM internally (slice 7b for
  constants, post-1.0 elliptic-integral families) inherits the
  `agm` dependency naturally.
- The CI feature matrix in `.github/workflows/ci.yml` gains two
  combinations: `--no-default-features --features=agm` (AGM alone)
  and the everything-on combo grows by `agm`.
- The kernel does not implement Ziv-strategy retry. The 64-bit
  guard above target precision is generous enough that the
  three-operation rounding error per iteration, compounded over
  the `O(log p)` iterations the loop actually runs, stays well
  below 1 ULP at the target precision. The
  `tests/differential_agm.rs` lane will catch any case where this
  guarantee breaks.

## References

- Plan: `let-s-review-the-backlog-vast-harbor.md` — slice 6l.
- ADR-0010 — multiplication algorithm thresholds (the AGM
  iteration calls `mul`, so AGM's running cost inherits the
  schoolbook / Karatsuba threshold tuning that lands in slice 7d).
- ADR-0014 — MPFR differential gating; the new
  `tests/differential_agm.rs` follows the same pattern as the 22
  existing differential test files.
- Brent, R. P. and Zimmermann, P. *Modern Computer Arithmetic*,
  §4.8 (AGM and its applications), Cambridge University Press,
  2010. — The Gauss iteration's convergence analysis and the
  Brent–Salamin variant's derivation both live here.
