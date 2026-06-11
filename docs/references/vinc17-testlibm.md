---
slug: vinc17-testlibm
category: corpus
citation: Lefevre, V. "TestLibm" maintained worst case database for correctly rounded binary64 elementary functions (successor to the 2000 Table Maker's Dilemma results pages).
edition: testlibm-data.xz as retrieved 2026-06-11 (792,872 bytes)
canonical_url: https://www.vinc17.net/research/testlibm/
document_number: none
doi: none
archived_url: http://web.archive.org/web/20260119140003/https://www.vinc17.net/research/testlibm/
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: 19ed0e6d6ec8e836a1d76818ff868bcc4fa47719ce713a342bfd8877c52a6861
vendored_path: none
license: not stated on the page; pfloat consumes this data only through CORE-MATH, whose files are MIT
vendor_status: pointer-only
rot_risk: single-maintainer
provenance_class: lineage
consumers:
  - docs/lefevre-muller-corpus-provenance.md
verification: none directly; the derived path is verified through tests/differential_lefevre_muller.rs (see core-math-wc-corpus.md).
---

# Lefevre's maintained TestLibm worst case database

## Why this source

The middle link of the hard to round corpus lineage: Lefevre and
Muller's 2001 search results, maintained and extended by Lefevre on his
personal site, are what CORE-MATH curates into the per-function `.wc`
files pfloat samples from. The corpus provenance document records the
full chain (L-M 2001 paper, this database, CORE-MATH `.wc` files, pfloat
subset).

## What it grounds

Nothing directly; pfloat never reads this database. The entry exists
because the chain of custody for the hard to round vectors runs through
it, and a single maintainer personal site is the most rot prone link in
that chain.

## Alternatives

The superseded 2000 snapshot at Muller's page (`muller-tmd-2000.md`) and
the consolidated tables in the Muller et al Handbook are the historical
and print forms of the same data; CORE-MATH is the maintained machine
readable form.

## Coverage gaps

The database covers binary64 unary elementary functions; it does not
cover the special functions (gamma, Bessel, zeta, Airy) or multi
argument kernels, so no hard to round vectors exist for those surfaces
from this lineage. See `core-math-wc-corpus.md` for the gaps in the
subset pfloat actually ships.

## Fixity note

The sha256 above is for `testlibm-data.xz` as retrieved on 2026-06-11
(792,872 bytes); the file is not vendored into this repository. The
Wayback snapshot of the data file
(`http://web.archive.org/web/20251008111500/https://www.vinc17.net/research/testlibm/testlibm-data.xz`)
was verified byte identical to the live copy (same sha256), so the hash
anchors both.
