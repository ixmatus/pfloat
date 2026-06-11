---
slug: arb
category: software
citation: Arb, a C library for rigorous real and complex arithmetic with arbitrary precision (Johansson), since 2023 part of FLINT. Consumed via python-flint in a subprocess worker.
edition: consumed via python-flint in the dedicated oracle venv
canonical_url: https://arblib.org/
document_number: none
doi: none
archived_url: http://web.archive.org/web/20260305102239/https://arblib.org/
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: LGPL 2.1 or later
vendor_status: pointer-only
rot_risk: single-maintainer
provenance_class: oracle
consumers:
  - scripts/arb_oracle_worker.py
  - scripts/setup_arb_oracle.sh
  - tests/oracle/arb.rs
  - pfloat-ball/tests/differential_arb.rs
  - pfloat-complex/tests/differential_acb.rs
  - docs/decisions/0078-pfloat-ball-verification.md
verification: the Arb lanes are themselves verification: the twelve special functions MPFR cannot cover (Si, Ci, li, Bi, Ai', Bi', Bessel I and K), the independent ball containment backstop (ADR-0078), and the acb componentwise certified rounding differential for pfloat-complex (2940 checks at the 1.0 cut).
---

# Arb (FLINT)

## Why this source

The only widely trusted oracle for rigorous arbitrary precision special
functions beyond MPFR's surface, and the only one that returns
enclosures rather than point estimates, which is what the ball
containment lane and the certified rounding bridge need. Johansson's
design paper is a separate entry (`johansson-arb-2017.md`): design
reference and oracle are different relationships.

## What it grounds

Three lanes: the Arb primary differential surface (the twelve FnIds
MPFR lacks, via the stdin/stdout worker protocol of ADR-0035), the
independent containment backstop for pfloat-ball (BRACKET verbs,
superset semantics for interval inputs), and the acb lane for
pfloat-complex. The worker runs out of process from a dedicated venv
(`scripts/setup_arb_oracle.sh`); python-flint is the binding because
rug has no MPC or Arb feature.

## Alternatives

mpmath (`mpmath.md`) covers similar functions without rigorous error
bounds (it is the tier two cross check); MPFR covers the standard
surface in process. Classed single-maintainer: Arb and FLINT are
overwhelmingly one author's work even after the merge.
