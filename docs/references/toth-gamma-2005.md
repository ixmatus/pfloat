---
slug: toth-gamma-2005
category: web
citation: Toth, V. T. "Programmable Calculators: The Gamma Function". rskey.org, 2005 (as cited in tree).
edition: none
canonical_url: https://www.rskey.org/gamma.htm
document_number: none
doi: none
archived_url: http://web.archive.org/web/20110927012411/http://www.rskey.org/gamma.htm
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: personal page, all rights reserved; consulted as a pattern, no code or coefficients copied
vendor_status: pointer-only
rot_risk: single-maintainer
provenance_class: primary
consumers:
  - src/math/gamma_stirling.rs
verification: tests/differential_lgamma.rs and differential_gamma.rs MPFR lanes cover the Spouge path the pattern informed; coefficients come from Spouge's formula, not this page.
---

# Toth, Programmable Calculators: The Gamma Function

## Why this source

`gamma_stirling.rs` cites this page as a reference implementation
pattern for evaluating the Spouge approximation in limited precision
settings: the practical shape of the sum (alternating term
accumulation, when to fold the leading factor) rather than any
mathematics beyond Spouge's paper.

## What it grounds

The implementation pattern of `spouge_lgamma` in
`src/math/gamma_stirling.rs`. The coefficients and the error bound come
from Spouge (1994) with Pugh's (2004) error analysis; this page
contributed engineering shape only.

## Alternatives

Numerical Recipes presents a Lanczos pattern instead; Spouge plus this
evaluation pattern was chosen for closed form coefficients with no
precomputed table (see `spouge-1994.md` when it lands and ADR-0041).

## Archive note

A hobbyist personal site maintained by one author since the 1990s;
classed single-maintainer. The recorded snapshot is from 2011 (verified
to carry the gamma article); the Wayback save endpoint was down at
mining time (2026-06-11), so a fresh save of the live page is pending
the annual re-verification pass.
