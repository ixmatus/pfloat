---
slug: brent-1976
category: paper
citation: Brent, R. P. "Fast Multiple-Precision Evaluation of Elementary Functions". Journal of the ACM 23(2), 1976, pp. 242-251. With Salamin's independent companion result (Salamin, "Computation of pi Using Arithmetic-Geometric Mean", Mathematics of Computation 30, 1976).
edition: none
canonical_url: https://doi.org/10.1145/321941.321944
document_number: none
doi: 10.1145/321941.321944
archived_url: none
archive_date: none
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: ACM copyright (Brent); AMS copyright (Salamin)
vendor_status: pointer-only
rot_risk: stable-publisher
provenance_class: primary
consumers:
  - src/math/agm.rs
  - src/math/agm_constants.rs
  - docs/decisions/0015-agm-formulation.md
  - docs/decisions/0017-agm-based-transcendental-constants.md
verification: tests/differential_agm.rs and differential_constants.rs; the 1024 bit constant pins (pi, ln 2, ln 10, 2 over pi, 2 over sqrt pi, ln 2pi) are three way pinned against mpmath.
---

# Brent (and Salamin), AGM evaluation of elementary functions and pi (1976)

## Why this source

The Brent and Salamin arithmetic geometric mean iteration: quadratic
convergence for pi and the elementary constants, the right asymptotic
tool for a library whose constants must be derivable at any precision
rather than tabulated. ADR-0015 records the AGM formulation; ADR-0017
the constant derivations built on it.

## What it grounds

`src/math/agm.rs` (the AGM iteration with its convergence cutoff) and
the cached constant derivations in `src/math/agm_constants.rs` (pi via
Brent and Salamin, ln 2 and ln 10 via the AGM log machinery).

## Alternatives

Machin like arctangent series converge linearly (fine at fixed small
precision, wrong asymptotics here); Chudnovsky converges faster for pi
specifically but is a single constant special case where the AGM
machinery serves the whole constant family. This entry deliberately
covers Salamin's companion paper in the citation line rather than a
separate file: the tree cites the joint algorithm name and nothing
cites Salamin's text independently.
