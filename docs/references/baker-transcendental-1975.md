---
slug: baker-transcendental-1975
category: book
citation: Baker, A. Transcendental Number Theory. Cambridge University Press, 1975. Chapters 1 and 2.
edition: none
canonical_url: none
document_number: none
doi: none
archived_url: none
archive_date: none
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: Cambridge University Press copyright
vendor_status: legally-cannot
rot_risk: stable-publisher
provenance_class: primary
consumers:
  - docs/decisions/0060-inexact-fidelity-transcendental-kernels.md
  - docs/decisions/0063-inexact-fidelity-elementary-transcendentals.md
  - docs/decisions/0064-inexact-fidelity-special-functions-proven.md
verification: the INEXACT fidelity lanes (the forced INEXACT dispatch tests across exp, log, trig, and the proven special function families) rest on the theorems this book states and proves.
---

# Baker, Transcendental Number Theory

## Why this source

The Lindemann Weierstrass and Gelfond Schneider theorems make the
INEXACT discipline unconditional rather than empirical: a
transcendental kernel's value on a nontrivial algebraic input is
transcendental, hence irrational, hence never exactly representable on
the target grid, so INEXACT can be forced on Ziv fall through without a
per input exactness proof. Chapters 1 and 2 carry the statements and
proofs the ADRs cite.

## What it grounds

The pre Ziv exact dispatch plus force INEXACT structure of the
elementary and special function kernels (ADR-0060, ADR-0063, ADR-0064),
including the scope split between proven clean families and the named
open problems (gamma at rational points, odd zeta values) where the
force is defensive only.

## Alternatives

Any standard transcendence text states the same theorems; Baker is
cited because it is the compact canonical monograph and the ADRs quote
its chapter numbers.

## Vendoring note

In print and in copyright; no canonical electronic location to archive.
