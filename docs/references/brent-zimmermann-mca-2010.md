---
slug: brent-zimmermann-mca-2010
category: book
citation: Brent, R. P. and Zimmermann, P. Modern Computer Arithmetic. Cambridge Monographs on Applied and Computational Mathematics 18, Cambridge University Press, 2010.
edition: free electronic version 0.5.9 (the authors' draft corresponding to the CUP edition)
canonical_url: https://maths-people.anu.edu.au/~brent/pub/pub226.html
document_number: none
doi: none
archived_url: http://web.archive.org/web/20260213054438/https://members.loria.fr/PZimmermann/mca/mca-cup-0.5.9.pdf
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: 0a7b9eff1f7c02faa7f3f350545282f1e5a5e34b055f8eb366f656c2ac3bc213
vendored_path: none
license: print edition CUP copyright; electronic version 0.5.9 distributed under CC BY-NC-ND 3.0 (stated in its front matter)
vendor_status: pointer-only
rot_risk: academic-personal
provenance_class: primary
consumers:
  - src/ops/limbs.rs
  - docs/decisions/0061-toom3-multiplication.md
  - docs/decisions/0051-formatter-magnitude-cap-and-sub-quadratic-conversion.md
  - DESIGN.md
verification: tests/property_mul.rs and the conversion property lanes exercise the Toom-3 and sub-quadratic base conversion paths derived with this book as reference.
---

# Brent and Zimmermann, Modern Computer Arithmetic (2010)

## Why this source

The standard modern reference for arbitrary precision integer and
floating point algorithms, by the authors of the GMP and MPFR
foundations. DESIGN.md names it (as MCA) a chapter level reference for
the arithmetic kernels; section 1.3.3 anchored the Toom-3 interpolation
cross check (ADR-0061) and section 1.7 the sub-quadratic radix
conversion reasoning (ADR-0051).

## What it grounds

The Toom-3 interpolation cross check in `src/ops/limbs.rs`, the
formatter magnitude cap and divide and conquer conversion rationale, and
the addition, subtraction, and transcendental range reduction
discussions DESIGN.md cites at section level.

## Alternatives

Knuth TAOCP volume 2 covers the classical algorithms (Algorithm D is
cited separately, see `knuth-taocp-v2.md`) but predates the modern
sub-quadratic toolbox; the Muller et al Handbook covers floating point
formats rather than multiprecision kernel algorithms.

## Archive note

The free PDF lives on two academic personal pages: Brent's ANU page
(`pub226.html`, the page DESIGN.md cites) and Zimmermann's Loria members
page (`https://members.loria.fr/PZimmermann/mca/mca-cup-0.5.9.pdf`).
Both are rot inventory items; both were archived at mining time. The
sha256 above is a fixity anchor for the 0.5.9 PDF (1,987,150 bytes,
retrieved 2026-06-11), recorded so a future mirror can be checked
against the copy this registry verified; the PDF is not vendored
(pointer plus archive posture, CC BY-NC-ND permitting but not requiring
a copy).
