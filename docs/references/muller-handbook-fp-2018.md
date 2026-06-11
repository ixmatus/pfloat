---
slug: muller-handbook-fp-2018
category: book
citation: Muller, J.-M., Brunie, N., de Dinechin, F., Jeannerod, C.-P., Joldes, M., Lefevre, V., Melquiond, G., Revol, N., Torres, S. Handbook of Floating-Point Arithmetic. 2nd edition, Birkhauser, 2018.
edition: 2nd edition
canonical_url: https://doi.org/10.1007/978-3-319-76526-6
document_number: none
doi: 10.1007/978-3-319-76526-6
archived_url: none
archive_date: none
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: Springer/Birkhauser copyright
vendor_status: legally-cannot
rot_risk: stable-publisher
provenance_class: contextual
consumers:
  - docs/lefevre-muller-corpus-provenance.md
  - DESIGN.md
verification: none derives from it; the corpus provenance document states explicitly that pfloat does not transcribe worst cases from the Handbook (CORE-MATH is the machine readable source).
---

# Muller et al, Handbook of Floating-Point Arithmetic

## Why this source

Two roles, both contextual. It is the print consolidation of the
Lefevre and Muller worst case tables (the corpus provenance document
names it in the lineage and records that pfloat does not transcribe
from it), and it is the second of the two free standing proxies the
crate level CLAUDE.md names for readers without IEEE 754 access
(chapter level treatment of formats, rounding, and exceptions).
DESIGN.md also cites it (as HFA) at chapter level for the arithmetic
kernel discussions.

## What it grounds

Nothing in the implementation directly; it is a bibliographic anchor
and reading aid. This entry, with `berkeley-testfloat.md`, is a
canonical example of the contextual class: catalogued precisely to
record a non derivation.

## Alternatives

Goldberg (`goldberg-1991.md`) is the short free proxy; the standard
itself (`ieee-754-2019.md`) is the normative source the code cites.
