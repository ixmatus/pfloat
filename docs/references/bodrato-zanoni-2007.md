---
slug: bodrato-zanoni-2007
category: paper
citation: Bodrato, M. and Zanoni, A. "Integer and Polynomial Multiplication: Towards Optimal Toom-Cook Matrices". ISSAC 2007 (Waterloo, Ontario), ACM, pp. 17-24.
edition: preprint "What About Toom-Cook Matrices Optimality?", Centro Vito Volterra N.605, 2006
canonical_url: http://www.bodrato.it/papers/#ISSAC2007
document_number: none
doi: 10.1145/1277548.1277552
archived_url: http://web.archive.org/web/20251026202924/http://www.bodrato.it/papers/
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: ACM copyright; author self-archived preprint free on bodrato.it
vendor_status: pointer-only
rot_risk: single-maintainer
provenance_class: primary
consumers:
  - src/ops/limbs.rs
  - docs/decisions/0061-toom3-multiplication.md
verification: tests/property_mul.rs and the limb multiplier differential coverage exercise the Toom-3 path above TOOM3_THRESHOLD; the interpolation was derived from the Vandermonde inverse and cross checked against this paper, not transcribed.
---

# Bodrato and Zanoni, Towards Optimal Toom-Cook Matrices (2007)

## Why this source

The optimal evaluation and interpolation sequences for Toom-Cook
multiplication: which point sets and which exact division orders
minimize the linear work around the recursive sub-products. The Toom-3
interpolation in `src/ops/limbs.rs` was derived from the Vandermonde
inverse at `{0, 1, -1, 2, inf}` and cross checked against this paper and
Modern Computer Arithmetic section 1.3.3 (ADR-0061); GMP's
`mpn_toom33_mul` was an oracle for behavior, not a template.

## What it grounds

The Toom-3 interpolation sequence and the exact division by three
(`divexact_by3`, after Jebelean) in `src/ops/limbs.rs`.

## Alternatives

Brent and Zimmermann section 1.3.3 covers the same interpolation at
textbook level (see `brent-zimmermann-mca-2010.md`); Bodrato's solo
WAIFI 2007 paper (LNCS 4547, pp. 116-133, "Towards Optimal Toom-Cook
Multiplication for Univariate and Multivariate Polynomials in
Characteristic 2 and 0") treats the polynomial characteristic 2 case the
integer kernel does not need.

## Citation defect in tree

The in tree citations (`src/ops/limbs.rs`, ADR-0061, the old
bibliography in `docs/references.md`) name this paper as "Bodrato and
Zanoni, WAIFI 2007". That conflates two works: the author pair and the
title belong to the ISSAC 2007 paper recorded here; WAIFI 2007 published
the Bodrato solo paper. This entry records the corrected citation; the
in tree comment fix is tracked as a separate concern (bead filed at
mining time).

## Archive note

The archived papers index above carries the ISSAC 2007 entry; the self
archived preprint PDF is also snapshotted
(`http://web.archive.org/web/20260115163657/http://marco.bodrato.it/papers/WhatAboutToomCookMatricesOptimality.pdf`,
263,445 bytes, verified 2026-06-11).
