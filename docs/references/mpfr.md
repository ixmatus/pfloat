---
slug: mpfr
category: software
citation: GNU MPFR, multiple precision floating point with correct rounding. Fousse, L., Hanrot, G., Lefevre, V., Pelissier, P., Zimmermann, P. "MPFR: A Multiple-Precision Binary Floating-Point Library with Correct Rounding". ACM TOMS 33(2), 2007.
edition: consumed via the rug crate (in process bindings)
canonical_url: https://www.mpfr.org/
document_number: none
doi: 10.1145/1236463.1236468
archived_url: http://web.archive.org/web/20260530152157/https://www.mpfr.org/
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: LGPL
vendor_status: pointer-only
rot_risk: community-run
provenance_class: oracle
consumers:
  - Cargo.toml
  - tests/differential_lefevre_muller.rs
  - tests/differential_zeta.rs
  - docs/decisions/0014-mpfr-differential-ci-gating.md
  - docs/decisions/0058-libm-verification-harness.md
verification: it IS the primary verification: the differential-mpfr feature gates roughly forty seven per function lanes plus conversion, parsing, and formatting oracles, in CI per push; pfloat-libm's harness is MPFR only by ADR-0058.
---

# GNU MPFR

## Why this source

The de facto semantics specification for correctly rounded arbitrary
precision floating point beyond IEEE 754's fixed formats, and pfloat's
principal behavioral oracle. Where 754 is silent (what correct rounding
means at precision p for a transcendental), MPFR's documented semantics
are the reference point; the TOMS 2007 paper is the citable statement
of those semantics.

## What it grounds

The differential verification surface: every MPFR covered function runs
a bit exact comparison lane under the `differential-mpfr` feature
(LGPL stays a dev dependency only, via `rug`; nothing links it at
runtime). Oracle, not template: behavior is cross checked, internals
were never copied, per the README provenance discipline.

Where pfloat deliberately diverges from MPFR semantics the tree records
it at the divergence site; the known classes are NaN payload handling
(payloads are not part of the correct rounding claim), roundTiesToAway
(MPFR has no RNDNA lane, so NearestAway oracles are synthesized), and
the exception flag model (pfloat carries IEEE status flags through
`Status`, MPFR's global flags differ in granularity).

## Alternatives

astro-float is the pure Rust arbitrary precision crate the tooling
defaults prefer for dev dependencies; it is not correctly rounded
across this surface, so MPFR remains the oracle of record where bit
exactness is the claim under test. Arb (`arb.md`) covers the special
functions MPFR lacks.
