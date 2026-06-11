---
slug: crlibm
category: software
citation: de Dinechin, Lauter, Muller and the Arenaire project. "CRlibm, a library of correctly rounded elementary functions in double-precision". ENS Lyon, 2004-2010.
edition: none
canonical_url: https://github.com/taschini/crlibm
document_number: none
doi: none
archived_url: http://web.archive.org/web/20260309194505/https://github.com/taschini/crlibm
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: LGPL
vendor_status: pointer-only
rot_risk: died-once
provenance_class: contextual
consumers:
  - pfloat-libm/README.md
  - docs/ROADMAP.md
verification: none; named as a non-template in the pfloat-libm disclosure, and its verification surface is compared (proofs, not enumeration) in the pfloat-libm README table.
---

# CRlibm

## Why this source

CRlibm is the proof carrying ancestor of correctly rounded libm work:
the first complete demonstration that correct rounding in binary64 is
practical, with Gappa assisted error proofs per function. The
pfloat-libm README names it twice: in the alternatives comparison table,
and in the provenance disclosure as a behavioral reference that is
explicitly not a code template (its LGPL license alone forecloses
adaptation into an MIT/Apache crate; see the license discipline in the
crate level CLAUDE.md).

## What it grounds

Nothing in the implementation. Context for pfloat-libm's positioning
(the README comparison table row: proven but not enumerated verification,
C not Rust, LGPL) and a named member of the class of C implementations
the agent is instructed to treat as oracles rather than templates.

## Alternatives

CORE-MATH is CRlibm's modern successor (MIT style license, active
maintenance) and is the project pfloat actually draws vectors from; MPFR
is the behavioral oracle for the arbitrary precision core.

## Archive note

Classed died-once: the original lipforge.ens-lyon.fr hosting is gone;
the canonical URL points at the longest lived community mirror. The
Arenaire project pages survive only partially. Archived at mining time.
