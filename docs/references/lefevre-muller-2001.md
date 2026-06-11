---
slug: lefevre-muller-2001
category: paper
citation: Lefevre, V. and Muller, J.-M. "Worst Cases for Correct Rounding of the Elementary Functions in Double Precision". ARITH-15, IEEE, 2001. The tree also cites the INRIA preprint (RR2000-35).
edition: none
canonical_url: https://doi.org/10.1109/ARITH.2001.930110
document_number: none
doi: 10.1109/ARITH.2001.930110
archived_url: none
archive_date: none
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: IEEE copyright; INRIA preprint distributed free via HAL
vendor_status: pointer-only
rot_risk: stable-publisher
provenance_class: lineage
consumers:
  - docs/lefevre-muller-corpus-provenance.md
  - tests/differential_lefevre_muller.rs
  - docs/decisions/0083-lefevre-muller-seeding-ball-generators.md
verification: the corpus lane tests/differential_lefevre_muller.rs runs the descendant data; the paper itself reaches pfloat only through CORE-MATH (see core-math-wc-corpus.md).
---

# Lefevre and Muller, Worst Cases for Correct Rounding (2001)

## Why this source

The origin of the hard to round corpus: the exhaustive search (by the
Lefevre algorithm) for binary64 inputs whose correctly rounded values
sit hardest against a rounding boundary. Everything downstream (the
maintained vinc17 database, CORE-MATH's `.wc` files, pfloat's sampled
subset, and the ball generator seeding of ADR-0083) descends from this
work.

## What it grounds

Classed lineage, not primary: pfloat transcribes no data from the paper
itself. It anchors the provenance chain recorded in
`docs/lefevre-muller-corpus-provenance.md` and the README's claim that
correctly rounded transcendentals stand on the published worst case
literature.

## Alternatives

The Stehle, Lefevre, Zimmermann lattice method
(`stehle-lefevre-zimmermann-2005.md`) is the successor search
algorithm; Ziv's strategy (`ziv-1991.md`) is the complementary route
pfloat's core actually computes by, which is why no worst case bound is
load bearing for correctness here.
