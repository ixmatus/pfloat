---
slug: burnikel-ziegler-1998
category: paper
citation: Burnikel, C. and Ziegler, J. "Fast Recursive Division". Research Report MPI-I-98-1-022, Max-Planck-Institut fuer Informatik, Saarbruecken, October 1998.
edition: none
canonical_url: https://pure.mpg.de/pubman/faces/ViewItemOverviewPage.jsp?itemId=item_1819444
document_number: MPI-I-98-1-022
doi: none
archived_url: http://web.archive.org/web/20240621153224/https://pure.mpg.de/pubman/faces/ViewItemOverviewPage.jsp?itemId=item_1819444
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: Max Planck Society research report, distributed free by the institute repository
vendor_status: pointer-only
rot_risk: community-run
provenance_class: primary
consumers:
  - src/ops/limbs.rs
  - docs/decisions/0052-recursive-burnikel-ziegler-division.md
verification: tests/differential_div.rs and the formatter and parser lanes that drive wide divisions through the recursive path above the threshold.
---

# Burnikel and Ziegler, Fast Recursive Division (1998)

## Why this source

The recursive division with remainder that reduces wide division to
multiplications (2K(n) + O(n log n) over Karatsuba time), making
division track whatever speed the multiplier has. ADR-0052 records the
adoption above the schoolbook threshold.

## What it grounds

The recursive divider in `src/ops/limbs.rs` (the D2n/1n and D3n/2n
split structure) feeding division, formatting, and parsing at large
precision.

## Alternatives

Newton inversion with multiplication wins asymptotically with FFT
multiplication, which pfloat does not have (FFT measured and deferred,
ADR-0010/ADR-0040); at Karatsuba and Toom-3 sizes the Burnikel and
Ziegler method is the standard choice (it is also GMP's mid range
divider, used as a behavioral oracle only).
