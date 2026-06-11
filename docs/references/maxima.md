---
slug: maxima
category: software
citation: Maxima, a computer algebra system (descendant of MIT Macsyma).
edition: invoked ad hoc via nix-shell for the pinned corpus tier
canonical_url: https://maxima.sourceforge.io/
document_number: none
doi: none
archived_url: http://web.archive.org/web/20260519151942/https://maxima.sourceforge.io/
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: GPL
vendor_status: pointer-only
rot_risk: community-run
provenance_class: oracle
consumers:
  - scripts/maxima_oracle_worker.py
  - scripts/maxima_oracle_worker.sh
  - tests/oracle/maxima.rs
verification: tier six sampling of the three way agreement protocol (ADR-0035): the pinned corpus entries carry a maxima triple check provenance class.
---

# Maxima

## Why this source

The third independent implementation in the agreement protocol, from a
lineage (Lisp computer algebra) that shares nothing with MPFR, Arb, or
mpmath. Slow per request startup confines it to sampling the pinned
corpus rather than running per push.

## What it grounds

The triple check provenance class on entries in `tests/oracle/pinned/`;
disagreement between the bigfloat oracles would surface here even if
both inherited a common bug.

## Alternatives

PARI/GP would serve the same third opinion role; Maxima was chosen for
its nix availability and its genuinely separate implementation
ancestry.
