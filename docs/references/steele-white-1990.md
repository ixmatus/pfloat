---
slug: steele-white-1990
category: paper
citation: Steele, G. L. and White, J. L. "How to Print Floating-Point Numbers Accurately". PLDI 1990, ACM, pp. 112-126.
edition: none
canonical_url: https://doi.org/10.1145/93542.93559
document_number: none
doi: 10.1145/93542.93559
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
  - src/fmt.rs
  - docs/decisions/0071-dragon4-shortest-formatter.md
verification: tests/fmt_shortest.rs round trip and shortness properties; the Dragon4 high fixup defect found in adversarial review is pinned by a regression case there.
---

# Steele and White, How to Print Floating-Point Numbers Accurately (1990)

## Why this source

The free format printing problem and its correct solution (the Dragon
family): emit the shortest decimal string that reads back to exactly
the source value. The shortest decimal formatter implements the
Dragon4 scaled integer approach for arbitrary precision mantissas
(ADR-0071).

## What it grounds

`to_shortest_decimal_string` in `src/fmt.rs`: the scaled fraction
setup, the digit loop, and the low/high boundary tracking.

## Alternatives

Grisu and Ryu are faster fixed format algorithms whose precomputed
tables assume binary32/64; they do not generalize to arbitrary
precision, so the Dragon approach (with Burger and Dybvig's
simplifications, see `burger-dybvig-1996.md`) is the right shape here.
ADR-0029 records the earlier deferral and ADR-0071 the landing.
