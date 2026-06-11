---
slug: dlmf
category: standard
citation: NIST Digital Library of Mathematical Functions. Olver, Lozier, Boisvert, Clark et al (eds). National Institute of Standards and Technology.
edition: Version 1.2.6 (current at retrieval; the DLMF is a versioned living edition of Abramowitz and Stegun's successor)
canonical_url: https://dlmf.nist.gov/
document_number: none
doi: none
archived_url: http://web.archive.org/web/20260610051014/https://dlmf.nist.gov/
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: US government work; NIST terms of use
vendor_status: pointer-only
rot_risk: stable-publisher
provenance_class: primary
consumers:
  - src/math/bessel_j.rs
  - src/math/bessel_i.rs
  - src/math/airy.rs
  - src/math/zeta.rs
  - src/math/ei.rs
  - src/math/beta.rs
  - src/math/trig_reciprocal.rs
verification: the per function differential lanes (MPFR where covered, Arb for the DLMF heavy surface: Si, Ci, li, Bi, Ai', Bi', Bessel I and K) verify the series, recurrences, and asymptotics taken from the cited chapters.
---

# NIST Digital Library of Mathematical Functions

## Why this source

The specification of record for the special function surface: series
definitions, recurrences, asymptotic expansions, normalizations, and
special values. Where IEEE 754 stops (it names gamma and erf special
cases but no methods), the DLMF is what the kernels derive from.

## What it grounds

Per chapter, as quoted in the modules: chapter 4.14 (reciprocal trig),
chapter 5 with 5.2/5.5.3/5.12.1 (gamma family and beta), chapter 6
(Ei, Si, Ci, li: convergent series 6.6.5, asymptotics 6.12), chapter 9
with 9.7.2 (Airy, including the u_k recurrence whose transcription
pitfall is documented in docs/references.md footnotes), chapter 10
(Bessel J/Y/I/K: Maclaurin, Miller recurrence directions, Hankel
asymptotics, 10.25.2/10.30/10.35.5/10.40.1), and chapter 25 with
25.2.3/25.4.2 (zeta series and functional equation). The function major
table in `docs/references.md` maps each function to its quoted
sections.

## Alternatives

Abramowitz and Stegun is the print ancestor (the DLMF supersedes it and
is maintained); Numerical Recipes and Gil/Segura/Temme are method
handbooks the kernels did not derive from and so get no entries. The
DLMF was chosen as the single special function authority precisely
because it is citable at equation granularity and versioned.
