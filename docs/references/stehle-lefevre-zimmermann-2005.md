---
slug: stehle-lefevre-zimmermann-2005
category: paper
citation: Stehle, D., Lefevre, V., Zimmermann, P. "Searching Worst Cases of a One-Variable Function Using Lattice Reduction". IEEE Transactions on Computers 54(3), 2005, pp. 340-346.
edition: none
canonical_url: https://doi.org/10.1109/TC.2005.55
document_number: none
doi: 10.1109/TC.2005.55
archived_url: none
archive_date: none
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: IEEE copyright
vendor_status: pointer-only
rot_risk: stable-publisher
provenance_class: lineage
consumers:
  - docs/lefevre-muller-corpus-provenance.md
verification: none directly; the corpus produced with its descendant tooling is verified through the corpus lane (see core-math-wc-corpus.md).
---

# Stehle, Lefevre, Zimmermann, Searching Worst Cases Using Lattice Reduction (2005)

## Why this source

The SLZ algorithm: the LLL based worst case search that extended the
original Lefevre search to wider precisions and is the basis of
CORE-MATH's BaCSeL tool, whose search parameters annotate the `.wc`
files pfloat samples. Named in the corpus provenance document as the
fourth link of the lineage.

## What it grounds

Nothing in the implementation; lineage for the corpus chain of custody.

## Alternatives

The original Lefevre search (`lefevre-muller-2001.md`) is the binary64
ancestor; this is its generalization.
