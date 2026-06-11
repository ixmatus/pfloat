---
slug: kahan-branch-cuts-1987
category: paper
citation: Kahan, W. "Branch Cuts for Complex Elementary Functions, or Much Ado About Nothing's Sign Bit". In Iserles and Powell (eds), The State of the Art in Numerical Analysis, Clarendon Press, Oxford, 1987.
edition: none
canonical_url: https://people.eecs.berkeley.edu/~wkahan/
document_number: none
doi: none
archived_url: http://web.archive.org/web/20260602171616/https://people.eecs.berkeley.edu/~wkahan/
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: copyrighted book chapter (Clarendon Press); author materials on the Berkeley page carry no stated license
vendor_status: pointer-only
rot_risk: academic-personal
provenance_class: primary
consumers:
  - pfloat-complex/src/csqrt.rs
  - docs/decisions/0091-complex-magnitude-phase-elementary-annex-g.md
verification: pfloat-complex/tests/annex_g_special_values.rs rows and the acb componentwise differential lane (pfloat-complex/tests/differential_acb.rs) exercise the csqrt branch behavior.
---

# Kahan, Branch Cuts for Complex Elementary Functions (1987)

## Why this source

Kahan's paper is the origin of signed zero branch cut discipline for the
complex elementary functions and of the cancellation robust `csqrt`
reformulation. C99 Annex G encodes its conclusions normatively, but the
reformulation pfloat-complex implements in the `csqrt` interior is the
paper's, not the standard's; Annex G specifies special values, not the
stable interior formula.

## What it grounds

The interior formula of `pfloat-complex/src/csqrt.rs` (the
`kahan_brackets` decomposition: compute `w` from `|z|` and `|a|` in the
cancellation free order, then derive the components by the branch of the
real part). ADR-0091 records the adoption.

## Alternatives

Direct evaluation of `sqrt((|z| + a) / 2)` without the robust ordering
loses half the significand near the negative real axis; the naive polar
form (`sqrt(r) cis(theta/2)`) double rounds through the argument. Both
rejected in ADR-0091.

## Archive note

The paper is a book chapter with no canonical electronic home; the
canonical URL points at Kahan's Berkeley page, the author's hosting for
his floating point papers and the rot inventory item this entry insures.
The volume itself is available only in print and through libraries.
