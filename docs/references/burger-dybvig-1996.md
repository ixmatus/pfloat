---
slug: burger-dybvig-1996
category: paper
citation: Burger, R. G. and Dybvig, R. K. "Printing Floating-Point Numbers Quickly and Accurately". PLDI 1996, ACM, pp. 108-116.
edition: none
canonical_url: https://doi.org/10.1145/231379.231397
document_number: none
doi: 10.1145/231379.231397
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
verification: tests/fmt_shortest.rs (shared with steele-white-1990).
---

# Burger and Dybvig, Printing Floating-Point Numbers Quickly and Accurately (1996)

## Why this source

The streamlined restatement of Steele and White's algorithm: the
cleaner loop structure, boundary handling, and termination conditions
that practical Dragon4 implementations follow. The formatter's loop
shape follows this presentation.

## What it grounds

The digit generation loop and boundary comparisons in `src/fmt.rs`
(with `steele-white-1990.md` as the problem statement and correctness
argument).

## Alternatives

See `steele-white-1990.md`: Grisu and Ryu are fixed format only.
