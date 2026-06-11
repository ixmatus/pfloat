---
slug: payne-hanek-1983
category: paper
citation: Payne, M. H. and Hanek, R. N. "Radian Reduction for Trigonometric Functions". ACM SIGNUM Newsletter 18(1), 1983, pp. 19-24.
edition: none
canonical_url: https://doi.org/10.1145/1057600.1057602
document_number: none
doi: 10.1145/1057600.1057602
archived_url: none
archive_date: none
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: ACM copyright
vendor_status: pointer-only
rot_risk: stable-publisher
provenance_class: primary
consumers:
  - src/math/trig_reduce.rs
verification: tests/differential_sin.rs and the trig lanes at huge arguments (the RC1 large argument precision fix class) exercise the reduction; the 2 over pi table precision scales with the argument exponent.
---

# Payne and Hanek, Radian Reduction for Trigonometric Functions (1983)

## Why this source

Accurate argument reduction modulo pi over 2 for huge arguments: use
only the relevant window of the 2 over pi expansion selected by the
argument's exponent, so reduction cost and accuracy do not degrade with
magnitude. The trig kernels reduce by this scheme before quadrant
dispatch.

## What it grounds

`src/math/trig_reduce.rs`: the windowed 2 over pi table walk and the
working precision rule that scales with the argument exponent (the RC1
review found and fixed an omitted exponent term in that rule; the
reduction structure itself is Payne and Hanek's).

## Alternatives

Cody and Waite reduction is cheaper and fine for small arguments (the
small argument path uses the same machinery degenerately); naive
division by pi over 2 at fixed precision loses all accuracy once the
argument exceeds the working precision, the classic sin(1e700) failure
the RC1 review pinned.
