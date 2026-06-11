---
slug: brent-mcmillan-1980
category: paper
citation: Brent, R. P. and McMillan, E. M. "Some New Algorithms for High-Precision Computation of Euler's Constant". Mathematics of Computation 34(149), 1980, pp. 305-312.
edition: none
canonical_url: https://doi.org/10.1090/S0025-5718-1980-0551307-4
document_number: none
doi: 10.1090/S0025-5718-1980-0551307-4
archived_url: none
archive_date: none
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: AMS copyright; Mathematics of Computation back issues are free on ams.org
vendor_status: pointer-only
rot_risk: stable-publisher
provenance_class: primary
consumers:
  - src/math/agm_constants.rs
  - docs/decisions/0018-euler-gamma-via-brent-mcmillan.md
verification: tests/differential_constants.rs pins Euler's constant against MPFR and the 1024 bit three way pinned decimals.
---

# Brent and McMillan, High-Precision Computation of Euler's Constant (1980)

## Why this source

The Bessel function identity algorithm (the B1 variant) for Euler's
constant: the fastest classical method whose error analysis gives an a
priori term count for any target precision. ADR-0018 records the
selection and the derivation.

## What it grounds

The Euler constant derivation in `src/math/agm_constants.rs` (the
Bessel identity sum pair with the term count proportional to the target
precision).

## Alternatives

The integral and series definitions converge far too slowly; the
zeta product accelerations are asymptotically worse than the Bessel
identity at high precision.
