---
slug: jebelean-1993
category: paper
citation: Jebelean, T. "An Algorithm for Exact Division". Journal of Symbolic Computation 15(2), 1993, pp. 169-180.
edition: none
canonical_url: https://doi.org/10.1006/jsco.1993.1012
document_number: none
doi: 10.1006/jsco.1993.1012
archived_url: none
archive_date: none
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: Academic Press/Elsevier copyright; open access on ScienceDirect under the JSC archive
vendor_status: pointer-only
rot_risk: stable-publisher
provenance_class: primary
consumers:
  - src/ops/limbs.rs
  - docs/decisions/0061-toom3-multiplication.md
verification: the Toom-3 lane of tests/property_mul.rs exercises divexact_by3 on every interpolation; an incorrect exact division would corrupt every product above the threshold.
---

# Jebelean, An Algorithm for Exact Division (1993)

## Why this source

Exact division by a known small divisor via the modular inverse, word
by word from the low end with no remainder bookkeeping. Toom-3
interpolation requires an exact division by three on one matrix row;
Jebelean's method makes it linear with a single limb multiply per word.

## What it grounds

`divexact_by3` in `src/ops/limbs.rs` (ADR-0061).

## Alternatives

General short division computes a remainder that is provably zero here,
wasted work the modular inverse method avoids; GMP's `mpn_divexact_by3`
implements the same idea and served only as a behavioral cross check.
