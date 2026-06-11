---
slug: berkeley-testfloat
category: software
citation: Hauser, J. "Berkeley TestFloat", a conformance test generator for IEEE 754 binary arithmetic in fixed formats (release 3 series).
edition: release 3e
canonical_url: http://www.jhauser.us/arithmetic/TestFloat.html
document_number: none
doi: none
archived_url: http://web.archive.org/web/20260520133238/http://jhauser.us/arithmetic/TestFloat.html
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: BSD-3-Clause (release 3 series)
vendor_status: pointer-only
rot_risk: single-maintainer
provenance_class: contextual
consumers:
  - README.md
  - docs/decisions/0094-reference-registry.md
verification: none; this entry documents why no third party conformance vector set runs against pfloat.
---

# Berkeley TestFloat

## Why this source

TestFloat is the de facto conformance vector generator for IEEE 754
implementations, which makes it the first thing a reviewer would expect
in pfloat's verification story. This entry records why it is absent:
TestFloat generates vectors for the fixed interchange formats (binary32,
binary64, binary128) and their format specific operations; pfloat
computes at arbitrary precision, where each precision is its own
"format" and no enumerable third party vector set exists or can exist
for the whole surface.

## What it grounds

The structural shape of the README disclosure's named failure mode
("inputs no test happened to cover"): with no third party conformance
corpus applicable beyond the fixed formats, pfloat's verification rests
on differential oracles (MPFR, Arb, mpmath), the hard to round corpus at
binary64, the exhaustive f32 enumeration in pfloat-libm, and property
plus Kani lanes. The gap this entry names feeds
`docs/references/coverage-gaps.md`.

## Alternatives

The IBM FPgen suite targets the same fixed formats with the same
inapplicability; differential testing against MPFR (whose semantics are
the de facto arbitrary precision extension of 754) is the route pfloat
takes instead. The fixed format boundary itself is exercised where it
exists: pfloat-libm's exhaustive f32 sweep enumerates all 2^32 binary32
inputs, which is stronger than sampled vectors on that surface.
