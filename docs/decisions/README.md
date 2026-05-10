# Architecture Decision Records

This directory holds the record of *why* pfloat is the way it is.
Each significant choice (numeric representation, API shape, feature
gating, verification posture, performance tradeoffs) gets one
Architecture Decision Record. Together they form the audit log a
future reviewer would otherwise have to reconstruct from commit
messages and release notes.

The format is borrowed from ferrodec, which borrowed it from the
broader ADR community.

## Conventions

- **Filenames**: `NNNN-short-slug.md`, four-digit zero-padded sequence
  number, lowercase slug. Numbers are never re-used; superseded ADRs
  keep their slot and link forward.
- **Format**: see `template.md`. Each ADR is short. A single page is
  the target; the form matters more than the length.
- **Status lifecycle**:
  - `proposed` — drafted, not yet acted on. Avoid this for
    retroactive ADRs.
  - `accepted` — the decision is in effect.
  - `superseded by ADR-NNNN` — replaced; keep the file as a
    historical record, link forward.
  - `rejected` — considered and decided against. Document for the
    next person who wonders the same thing.
- **Plans**: planning artifacts archive under `plans/` with a date
  prefix (`YYYY-MM-DD-slug.md`). They are snapshots of the state at
  decision time, not living documents. ADRs reference the plan that
  produced them when applicable.

## Writing a new ADR

1. Pick the next available number.
2. Copy `template.md` to `NNNN-your-slug.md`.
3. Fill in: status, date, context, decision, consequences, related
   references.
4. If the decision supersedes a prior one, edit the prior ADR's
   status line to `superseded by ADR-NNNN`.

Decisions that are reversible or local in scope do not need an ADR.
These are for choices that matter to future contributors deciding
whether to revisit a path.

## Index

The Phase 0 ADRs (ratified as `accepted` on 2026-05-10):

- [0001 — `u64` limb representation, sign-magnitude, top-bit-set normalization](0001-limb-representation.md)
- [0002 — Bit-level precision granularity](0002-bit-level-precision.md)
- [0003 — Dual API: `BigFloat` (dynamic) and `FixedFloat<PREC>` (const-generic)](0003-dual-api.md)
- [0004 — Mantissa storage: `Vec<u64>` for `BigFloat`, `[u64; N]` for `FixedFloat`](0004-mantissa-storage.md)
- [0005 — Special-value encoding via tagged `Class` enum](0005-class-enum.md)
- [0006 — `i64` exponent](0006-exponent-type.md)
- [0007 — Rounding mode: function-call enum plus typestate; flags via thread-local in std, passed-context in no_std](0007-rounding-and-flags.md)
- [0008 — Differential testing oracle: `gmp-mpfr-sys` on a feature-gated CI lane](0008-differential-oracle.md)
- [0009 — Verification scaffolding: copy-paste from ferrodec, no shared crate](0009-verification-scaffolding.md)
- [0010 — Schönhage-Strassen FFT multiplication deferred to 1.x](0010-fft-deferred.md)
- [0011 — MSRV moves to nightly to use `generic_const_exprs`](0011-msrv-nightly-for-generic-const-exprs.md)
