---
slug: iso-c-annex-g
category: standard
citation: ISO/IEC 9899:2011 (C11), Annex G: IEC 60559-compatible complex arithmetic. Cited via the freely available committee draft N1570.
edition: N1570 committee draft (2011-04-12), whose Annex G matches the published C11 text
canonical_url: https://www.open-std.org/jtc1/sc22/wg14/www/docs/n1570.pdf
document_number: N1570
doi: none
archived_url: http://web.archive.org/web/20260609091213/https://open-std.org/jtc1/sc22/wg14/www/docs/n1570.pdf
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: ISO copyright; the N1570 working draft is distributed free by WG14
vendor_status: pointer-only
rot_risk: stable-publisher
provenance_class: primary
consumers:
  - pfloat-complex/src/specials.rs
  - pfloat-complex/src/csqrt.rs
  - pfloat-complex/src/cexp.rs
  - pfloat-complex/src/clog.rs
  - docs/decisions/0090-complex-multiply-divide.md
  - docs/decisions/0091-complex-magnitude-phase-elementary-annex-g.md
verification: pfloat-complex/tests/annex_g_special_values.rs enumerates the Annex G special value rows; pfloat-complex/tests/dispatch_totality.rs proves the dispatch covers every input class.
---

# C11 Annex G (via N1570)

## Why this source

Annex G is the normative encoding of branch cut and special value
discipline for complex elementary functions (the standardization of
Kahan's 1987 analysis, see `kahan-branch-cuts-1987.md`): signed zero on
branch cuts, the complex infinity recovery rules of G.5.1, and the
special value tables pfloat-complex implements row by row.

## What it grounds

The special value dispatch of every pfloat-complex kernel (multiply,
divide, magnitude, phase, sqrt, exp, log) and the G.5.1 infinity
recovery behavior. ADR-0090 and ADR-0091 record the derivations; the
canonical Annex G derivation document for the v1.x tail work was built
from the verbatim N1570 text.

## Alternatives

The published ISO C11 standard is paywalled; N1570 is the standing free
draft whose Annex G text matches it, which is why the tree cites N1570
clause numbers. C99's Annex G (via N1256) is materially identical for
the surface pfloat-complex implements.
