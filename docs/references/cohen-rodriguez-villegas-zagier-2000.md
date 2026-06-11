---
slug: cohen-rodriguez-villegas-zagier-2000
category: paper
citation: Cohen, H., Rodriguez Villegas, F., Zagier, D. "Convergence Acceleration of Alternating Series". Experimental Mathematics 9(1), 2000, pp. 3-12.
edition: none
canonical_url: https://doi.org/10.1080/10586458.2000.10504632
document_number: none
doi: 10.1080/10586458.2000.10504632
archived_url: none
archive_date: none
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: Taylor and Francis copyright (Experimental Mathematics)
vendor_status: pointer-only
rot_risk: stable-publisher
provenance_class: primary
consumers:
  - src/math/zeta.rs
  - docs/decisions/0026-zeta.md
verification: tests/differential_zeta.rs and property_zeta.rs; the acceleration's error constant feeds the zeta working precision schedule.
---

# Cohen, Rodriguez Villegas, Zagier, Convergence Acceleration of Alternating Series (2000)

## Why this source

The general alternating series acceleration whose Chebyshev polynomial
weights give geometric convergence with an explicit `(3 + sqrt 8)^(-n)`
error constant. The zeta kernel uses it alongside Borwein's
specialization (`borwein-zeta-1995.md`); ADR-0026 records how the two
compose in the dispatch.

## What it grounds

The acceleration weights and the term count selection in
`src/math/zeta.rs` for the eta series region.

## Alternatives

Euler transformation converges too slowly for the precision schedule;
Euler-Maclaurin needs Bernoulli number infrastructure (the MPFR route)
rejected in ADR-0026.
