---
slug: goldberg-1991
category: paper
citation: Goldberg, D. "What Every Computer Scientist Should Know About Floating-Point Arithmetic". ACM Computing Surveys 23(1), 1991, pp. 5-48.
edition: the Sun Microsystems edited reprint (with the differences from the ACM original documented in its appendix)
canonical_url: https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
document_number: none
doi: 10.1145/103162.103163
archived_url: http://web.archive.org/web/20260608130707/https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: ACM copyright; the edited reprint is distributed free by Sun/Oracle with permission
vendor_status: pointer-only
rot_risk: died-once
provenance_class: contextual
consumers:
  - CLAUDE.md
verification: none; this entry is proxy literature, no test derives from it.
---

# Goldberg, What Every Computer Scientist Should Know About Floating-Point Arithmetic (1991)

## Why this source

IEEE 754-2019 is paywalled, so the registry convention names free proxy
literature for it. Goldberg's survey is the canonical free introduction
to 754 semantics: rounding, guard digits, the standard's rationale. The
crate level CLAUDE.md names it (with the Muller et al Handbook) as the
754 proxy pair.

## What it grounds

Nothing in the implementation derives from it; pfloat's 754 behavior
traces to the standard's clauses directly (see the clause table in
`docs/references.md`). The entry exists so a reader without IEEE access
has a verifiable path into the semantics the code cites.

## Alternatives

The Muller et al Handbook of Floating-Point Arithmetic covers the same
ground at book depth (see `muller-handbook-fp-2018.md`); Goldberg is the
free and shorter on ramp.

## Archive note

Classed died-once: the paper's free hosting has moved repeatedly over
three decades (Sun documentation servers, validlab, assorted university
mirrors). The Oracle documentation copy is the longest lived edited
reprint; the archived snapshot is the rot insurance.
