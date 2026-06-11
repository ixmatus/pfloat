---
slug: knuth-taocp-v2
category: book
citation: Knuth, D. E. The Art of Computer Programming, Volume 2, Seminumerical Algorithms. 3rd edition, Addison-Wesley, 1998.
edition: 3rd edition
canonical_url: none
document_number: none
doi: none
archived_url: none
archive_date: none
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: Addison-Wesley copyright; print only
vendor_status: legally-cannot
rot_risk: stable-publisher
provenance_class: primary
consumers:
  - src/ops/limbs.rs
  - docs/decisions/0001-limb-representation.md
verification: tests/property_mul.rs and tests/differential_mul.rs exercise the base case division and multiplication paths; the decimal divmod work (parse OOM hardening) is word wise Knuth Algorithm D by design.
---

# Knuth, TAOCP Volume 2, Seminumerical Algorithms

## Why this source

Section 4.3.1's Algorithm D is the canonical base case long division
for multi limb integers; the limb kernel implements it below the
recursive Burnikel and Ziegler threshold. The limb representation
itself (ADR-0001) was reasoned against chapter 4's treatment of
positional arithmetic.

## What it grounds

The base case division in `src/ops/limbs.rs` and the word wise decimal
divmod posture (bit at a time division was rejected as O(exp squared)
during the parse hardening work).

## Alternatives

Brent and Zimmermann's Modern Computer Arithmetic covers the same
ground with modern notation (`brent-zimmermann-mca-2010.md`); Algorithm
D is cited from Knuth because that is the formulation the code follows
and names.

## Vendoring note

In print and in copyright; no canonical electronic location exists to
archive. A paper copy is the durable form.
