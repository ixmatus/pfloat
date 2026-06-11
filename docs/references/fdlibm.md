---
slug: fdlibm
category: software
citation: fdlibm, the Freely Distributable LIBM (Sun Microsystems/SunSoft, developed by K.C. Ng et al), hosted at Netlib.
edition: none
canonical_url: https://www.netlib.org/fdlibm/
document_number: none
doi: none
archived_url: http://web.archive.org/web/20260608125239/https://www.netlib.org/fdlibm/
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: Sun permission notice (free use and redistribution with the notice preserved)
vendor_status: pointer-only
rot_risk: community-run
provenance_class: contextual
consumers:
  - pfloat-libm/README.md
verification: none; named as a non-template in the pfloat-libm disclosure. Its accuracy class (faithful, not correctly rounded) is the contrast row in the README comparison table via the libm crate, its Rust descendant.
---

# fdlibm

## Why this source

The ancestor of nearly every system libm (and of the Rust `libm`
crate): the canonical faithful but not correctly rounded double
precision implementation. The pfloat-libm disclosure names it, with
CRlibm, as the class of C implementations the agent treats as
behavioral references rather than templates; the comparison table
positions pfloat-libm against its descendants.

## What it grounds

Nothing in the implementation. Context for what pfloat-libm is not
(faithful rounding) and provenance hygiene (named so the non copying
instruction is checkable against a concrete artifact).

## Alternatives

CRlibm (`crlibm.md`) and CORE-MATH are the correctly rounded C
lineages; fdlibm is catalogued as the faithful one.
