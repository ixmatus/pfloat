---
slug: marcovecchio-log2-2009
category: paper
citation: "Raffaele Marcovecchio, The Rhin–Viola method for log 2, Acta Arithmetica 139 (2009), 147–184."
edition: Acta Arithmetica 139 (2009)
canonical_url: https://doi.org/10.4064/aa139-2-5
document_number: none
doi: 10.4064/aa139-2-5
archived_url: http://web.archive.org/web/20260416220120/https://en.wikipedia.org/wiki/Irrationality_measure
archive_date: 2026-04-16
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: publisher-restricted (IMPAN); pointer only
vendor_status: pointer-only
rot_risk: died-once
provenance_class: primary
consumers:
  - src/math/exp.rs
  - docs/decisions/0096-exp-exponent-ceiling.md
verification: the bound value 3.57455 was cross-checked against the
  Wikipedia irrationality-measure table (live and the archived
  snapshot above, both eyeballed 2026-06-11); the paper metadata
  (title, author, journal, volume, pages) was verified via the
  Crossref API (api.crossref.org/works/10.4064/aa139-2-5). The
  consuming cap derivation in ADR-0096 uses factor 4 > mu, so a
  final-decimal error in the bound is absorbed by the margin.
---

# Marcovecchio 2009: the irrationality measure of log 2

Rot note (the `died-once` class): the DOI already resolves to the
IMPAN homepage rather than the article as of the retrieval date — the
per-article redirect is broken at the publisher. Wayback SPN was down
on 2026-06-11 and no snapshot of the DOI or publisher page exists, so
the archived URL above is the verified existing snapshot of the
secondary table that carries the bound value, not the paper itself; a
save of the canonical pointer is owed at the annual re-verification
ritual. The Crossref metadata record is the surviving canonical
identifier.

## Why this source

The best published upper bound on the irrationality measure of
`log 2`: `μ(log 2) ≤ 3.57455…`, improving Rhin–Viola. pfloat needs an
effective bound, not the exact decimal: the measure caps how deeply a
dyadic input `x` can agree with `k·ln 2`, which is exactly the
cancellation depth of `exp`'s argument reduction `r = x − k·ln 2`.

## What it grounds

Two derivations in the exp exponent-rim path (ADR-0096, pf-7z66):
the realized-collapse retry cap in `exp_reduced_pinned`
(`4·(precision(x) + 64) + 1024` extra reduction bits suffices because
the agreement depth is bounded by `μ·(precision(x) + 64)` bits and
`4 > μ`), and the termination argument for the certified
`floor(x/ln2)` interval loop (`x/ln2` for a ≤64-bit-span dyadic
cannot sit closer to an integer than `~2^-230`, so `q = 1024`
certifies).

## Alternatives

Rhin–Viola's earlier bound (~3.89) would also suffice for the factor-4
margin; the cap derivation deliberately carries enough headroom that
the choice between published bounds is not load-bearing.
