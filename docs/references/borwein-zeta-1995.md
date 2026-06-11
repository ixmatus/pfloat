---
slug: borwein-zeta-1995
category: paper
citation: Borwein, P. "An Efficient Algorithm for the Riemann Zeta Function". 1995 preprint; published in Constructive, Experimental, and Nonlinear Analysis, CMS Conference Proceedings 27, 2000, pp. 29-34.
edition: CECM preprint P155
canonical_url: https://www.cecm.sfu.ca/personal/pborwein/PAPERS/P155.pdf
document_number: none
doi: none
archived_url: http://web.archive.org/web/20210414133849/http://www.cecm.sfu.ca/personal/pborwein/PAPERS/P155.pdf
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: author preprint on a personal academic page; CMS holds the proceedings copyright
vendor_status: pointer-only
rot_risk: died-once
provenance_class: primary
consumers:
  - src/math/zeta.rs
  - docs/decisions/0026-zeta.md
verification: tests/differential_zeta.rs MPFR lane; the high precision composing check at p=1024 (sub-slice 2b.1 lesson) exercises the Borwein sum past the dispatch threshold.
---

# Borwein, An Efficient Algorithm for the Riemann Zeta Function (1995)

## Why this source

Borwein's algorithm 2 is the alternating series acceleration the zeta
kernel implements for the convergent region: Chebyshev shaped weights
with the `2 (3 + sqrt 8)^(-n)` error bound, which makes the working
precision schedule derivable in closed form. ADR-0026 records the
selection.

## What it grounds

The Borwein sum core of `src/math/zeta.rs` (weights, term count from the
target precision, and the error bound the Ziv envelope consumes),
composed with the DLMF 25.4.2 functional equation for the reflected
region.

## Alternatives

Euler-Maclaurin summation (the MPFR route) needs Bernoulli numbers and
tail bounds that are harder to certify at arbitrary precision;
Cohen-Rodriguez-Villegas-Zagier acceleration (see
`cohen-rodriguez-villegas-zagier-2000.md`) is the sibling method the
kernel also draws on. ADR-0026 documents the comparison.

## Archive note

Classed died-once: the canonical copy lives on Peter Borwein's CECM
personal page; the author died in 2020 and the host already refuses
non-browser clients. The Wayback snapshot is the durable pointer.
