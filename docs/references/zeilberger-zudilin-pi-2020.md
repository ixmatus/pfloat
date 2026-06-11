---
slug: zeilberger-zudilin-pi-2020
category: paper
citation: "Doron Zeilberger and Wadim Zudilin, The irrationality measure of pi is at most 7.103205334137..., Moscow Journal of Combinatorics and Number Theory 9 (2020), 407-419."
edition: Moscow J. Comb. Number Theory 9 (2020)
canonical_url: https://doi.org/10.2140/moscow.2020.9.407
document_number: arXiv:1912.06345
doi: 10.2140/moscow.2020.9.407
archived_url: none
archive_date: none
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: publisher-restricted (MSP); arXiv preprint openly readable; pointer only
vendor_status: pointer-only
rot_risk: stable-publisher
provenance_class: primary
consumers:
  - src/math/trig_reduce.rs
  - docs/decisions/0098-input-structure-aware-dispatch.md
verification: the bound value 7.103205334137 was cross-checked against
  the Wikipedia irrationality-measure table (retrieved 2026-06-11,
  which cites this paper) and the arXiv abstract title (1912.06345,
  retrieved 2026-06-11); paper metadata verified via the Crossref API
  (api.crossref.org/works/10.2140/moscow.2020.9.407). The consuming
  cap derivation in trig_reduce uses factor 8 > mu, so a
  final-decimal error in the bound is absorbed by the margin.
---

# Zeilberger-Zudilin 2020: the irrationality measure of pi

Wayback note: SPN was down on 2026-06-11 and no snapshot of the DOI
or arXiv page exists yet; a save is owed at the annual
re-verification ritual. The arXiv identifier (1912.06345) and the
Crossref record are the durable pointers; rot risk is classed
stable-publisher on the strength of the arXiv + MSP + Crossref
triple coverage.

## Why this source

The best published upper bound on the irrationality measure of `π`:
`μ(π) ≤ 7.103205334137…`, improving Salikhov. pfloat needs an
effective bound, not the exact decimal: the measure caps how deeply
a dyadic input `x` can agree with `q·(π/2)`, which is exactly the
collapse depth of the trig argument-reduction residual.

## What it grounds

The termination cap of the realized-collapse retry loop in
`trig_reduce::reduce` (ADR-0098, pf-k68i): the product width grows
on a collapsed residual up to `8·(precision(x) + e_x) + working +
256` bits, sufficient because the residual magnitude is bounded
below by `~2^-(μ·(precision(x) + e_x) + c)` and `8 > μ` with slack.

## Alternatives

Salikhov's earlier bound (7.6063…) would also suffice for the
factor-8 margin; the cap deliberately carries enough headroom that
the choice between published bounds is not load-bearing.
