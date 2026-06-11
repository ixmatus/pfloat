---
slug: mpc
category: software
citation: GNU MPC, a C library for multiprecision complex arithmetic with exact rounding (Enge, Gastineau, Pelissier, Zimmermann).
edition: design reference only; not a build or test dependency
canonical_url: https://www.multiprecision.org/mpc/
document_number: none
doi: none
archived_url: http://web.archive.org/web/20260501014922/https://www.multiprecision.org/mpc/
archive_date: 2026-06-11
retrieval_date: 2026-06-11
sha256: none
vendored_path: none
license: LGPL
vendor_status: pointer-only
rot_risk: community-run
provenance_class: primary
consumers:
  - pfloat-complex/README.md
  - docs/decisions/0092-complex-verification-posture.md
verification: pfloat-complex's independent differential lane runs against Arb's acb, not MPC (rug exposes no MPC binding), so the design reference and the test oracle are deliberately different implementations.
---

# GNU MPC

## Why this source

The de facto semantics reference for componentwise correctly rounded
complex arithmetic at arbitrary precision: what it means for each of
the real and imaginary components to be correctly rounded under its own
rounding mode. pfloat-complex adopted the componentwise model with MPC
as the named design reference (and not a code template), recorded in
the crate README and ADR-0092.

## What it grounds

The componentwise correct rounding contract of pfloat-complex's
arithmetic surface. Verification deliberately uses a different
implementation (Arb's acb, see `arb.md`): the design reference defines
the contract, an independent implementation checks it, so a shared
misreading cannot self confirm.

## Alternatives

C11 Annex G defines special values but not a rounding model at
arbitrary precision; Arb's acb is ball based rather than componentwise
correctly rounded, which is exactly why it works as the independent
check rather than the contract definition.
