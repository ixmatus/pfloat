---
slug: mpmath
category: software
citation: mpmath, a Python library for real and complex floating point arithmetic with arbitrary precision (Johansson et al).
edition: consumed from the shared oracle venv; also the offline generator of the corpus expected outputs at 200 bits
canonical_url: https://mpmath.org/
document_number: none
doi: none
archived_url: http://web.archive.org/web/20260609032917/https://mpmath.org/
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: BSD
vendor_status: pointer-only
rot_risk: community-run
provenance_class: oracle
consumers:
  - scripts/mpmath_oracle_worker.py
  - scripts/extract_lefevre_muller.py
  - tests/oracle/cross_check.rs
verification: tier two of the three way agreement protocol (ADR-0035); independently computes the hard to round corpus expected outputs and the pinned 1024 bit constant decimals.
---

# mpmath

## Why this source

The independent non MPFR arbitrary precision implementation: agreement
between MPFR, Arb, and mpmath is meaningful precisely because mpmath
shares no code with the other two. It also generated the corpus
expected outputs (at 200 bits, rounded to binary64) and the pinned
constant decimals, so the hard to round lane cross checks mpmath
against pfloat rather than trusting any single upstream.

## What it grounds

The cross check tier of the oracle subsystem and the offline generation
provenance of `tests/differential/lefevre_muller_data.rs` and the
constant pins.

## Alternatives

Maxima (`maxima.md`) is the slower tier beyond it; SageMath bundles the
same mpmath and adds nothing independent.
