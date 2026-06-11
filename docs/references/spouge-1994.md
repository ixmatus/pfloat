---
slug: spouge-1994
category: paper
citation: Spouge, J. L. "Computation of the Gamma, Digamma, and Trigamma Functions". SIAM Journal on Numerical Analysis 31(3), 1994, pp. 931-944.
edition: none
canonical_url: https://doi.org/10.1137/0731050
document_number: none
doi: 10.1137/0731050
archived_url: none
archive_date: none
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: SIAM copyright
vendor_status: pointer-only
rot_risk: stable-publisher
provenance_class: primary
consumers:
  - src/math/gamma_stirling.rs
  - src/math/lgamma.rs
  - docs/decisions/0041-spouge-precision-pegging.md
verification: tests/differential_gamma.rs, differential_lgamma.rs, and differential_digamma.rs MPFR lanes; the precision gated dispatch above the Spouge threshold is verified by composing kernels at high precision (the sub-slice 2b.1 lesson).
---

# Spouge, Computation of the Gamma, Digamma, and Trigamma Functions (1994)

## Why this source

Spouge's approximation gives gamma family values with coefficients in
closed form (`c_k` from factorials and powers, no precomputed table)
and an error bound that is a simple function of the parameter `a`,
which makes it derivable at any working precision. That is the property
the arbitrary precision dispatch needs and the reason it was chosen
over Lanczos.

## What it grounds

`spouge_lgamma` and `spouge_a_for` in `src/math/gamma_stirling.rs`
(coefficient computation, the `a log2 a >= p` selection rule, and the
truncation error bound the Ziv envelope consumes), with Pugh's error
analysis (`pugh-2004.md`) backing the bound and ADR-0041 recording the
precision pegging.

## Alternatives

Lanczos approximation needs tabulated coefficients per precision
(rejected for the closed form property); Stirling's series with
recurrence shift handles the large argument region and is the other arm
of the dispatch, derived from DLMF chapter 5.
