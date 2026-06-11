---
slug: johansson-arb-2017
category: paper
citation: Johansson, F. "Arb: Efficient Arbitrary-Precision Midpoint-Radius Interval Arithmetic". IEEE Transactions on Computers 66(8), 2017, pp. 1281-1292.
edition: also arXiv:1611.02831
canonical_url: https://doi.org/10.1109/TC.2017.2690633
document_number: none
doi: 10.1109/TC.2017.2690633
archived_url: http://web.archive.org/web/20260401000258/https://arxiv.org/abs/1611.02831
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: IEEE copyright; the arXiv preprint is free
vendor_status: pointer-only
rot_risk: stable-publisher
provenance_class: primary
consumers:
  - pfloat-ball/src/spec.rs
  - pfloat-ball/src/mag.rs
  - docs/decisions/0074-pfloat-ball-crate-and-mag.md
verification: pfloat-ball property and differential lanes (including the independent Arb containment lane, pfloat-ball/tests/differential_arb.rs) verify the enclosure laws the design inherits; the deliberate divergences are tested where they bite (64 bit radius tightness).
---

# Johansson, Arb (2017)

## Why this source

The principal design reference for midpoint and radius ball arithmetic:
why balls beat inf sup intervals at high precision (one rounding
direction, half the memory, asymptotically tight radii), the `mag_t`
unsigned radius type, and the accuracy accounting. pfloat-ball names
Johansson's design by behavior and identifiers; this entry supplies the
bibliographic anchor.

## What it grounds

The ball representation and the enclosure law statement in
`pfloat-ball/src/spec.rs`, and the `Mag` radius type in
`pfloat-ball/src/mag.rs`. The deliberate divergence is recorded in
ADR-0074: a 64 bit radius significand where Arb uses 30 bits, chosen
for tightness and Kani verifiability. Arb the software is also an
oracle (separate entry, `arb.md`): design reference and oracle are
different relationships and are recorded separately.

## Alternatives

MPFI's inf sup arbitrary precision intervals (Revol and Rouillier) are
the contrasting representation, conformant to 1788 style semantics but
paying two roundings per endpoint; van der Hoeven's ball arithmetic
papers are the other lineage root. Neither is load bearing for
pfloat-ball's code, so they are named here rather than given entries.
