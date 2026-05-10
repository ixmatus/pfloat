# ADR-0009: Verification scaffolding, copy-paste from ferrodec, no shared crate

- **Status**: accepted
- **Date**: 2026-05-10

## Context

ferrodec ships a verification stack: Kani harnesses under
`src/verify/`, conformance tests under `tests/`, fuzz harnesses under
`fuzz/`, ADR template and conventions under `docs/decisions/`,
KNOWN_ISSUES.md, CHANGELOG.md, and a CI workflow that ties them
together. The scaffolding works.

Two ways to share scaffolding between ferrodec and pfloat:

- **Shared dev-dependency crate** (`pfloat-verify-scaffolding` or
  similar). Each project depends on it and gets harness templates,
  CI workflow snippets, ADR conventions for free. Updates to the
  shared crate propagate.
- **Copy-paste at bootstrap time**. pfloat starts with a snapshot
  of ferrodec's scaffolding and evolves it independently.

The shared-crate approach is the more "engineered" choice. It is
also the more coupled choice. Once two projects share scaffolding,
changes to either's needs become harder to make: the shared crate
either grows abstractions to accommodate both, or one project ends
up driving the API at the other's expense.

ferrodec and pfloat differ on important axes:

- **Numeric domain**: decimal vs binary. Different conformance
  corpora (decTest vs Lefèvre–Muller and IEEE binary FP suites).
- **Type architecture**: single type (`Decimal128`) vs dual type
  (`BigFloat` + `FixedFloat<PREC>`).
- **Surface scope**: ferrodec is feature-complete at 1.x; pfloat
  has six to nine months of phases ahead, each of which will
  surface scaffolding needs that ferrodec never had.
- **Differential oracle**: ferrodec uses `astro-float` and direct
  calculation; pfloat uses MPFR via `gmp-mpfr-sys`.

The honest read: the scaffolding is more useful as a starting
point than as a binding. Coupling now would tax pfloat's evolution
through every phase.

## Decision

Copy-paste from ferrodec at bootstrap time. Maintain no shared
crate.

The bootstrap commit imports:

- ADR template (`docs/decisions/template.md`).
- ADR README and conventions (`docs/decisions/README.md`).
- Cargo workspace lints config (in `Cargo.toml`).
- `rust-toolchain.toml` (stable + rustfmt + clippy + thumbv6m
  targets).
- `.gitignore` (Rust standard plus `.claude/`).
- CI workflow shape (`.github/workflows/ci.yml`), trimmed to
  pfloat's feature surface.

The bootstrap commit does not import:

- `src/verify/` Kani harnesses. Phase 6 lands harnesses adapted to
  pfloat's surface (transcendentals, special functions, dual API).
  The directory shape will follow ferrodec's, but the proofs will
  not transfer literally.
- Conformance corpus. Different domain.
- Fuzz harnesses. Phase 6 task.
- KNOWN_ISSUES.md, CHANGELOG.md. Empty until 1.0 work begins.

When ferrodec evolves its scaffolding (a CI workflow update, a
clippy lint allowance), pfloat does not automatically inherit. If
the change is relevant, it is copied across by hand. The cost is
small at the volume of changes either project makes.

## Consequences

**Wins:**

- pfloat starts with a known-working scaffolding shape on day one.
  No bikeshedding about ADR format, lint allowances, or CI matrix
  structure.
- pfloat evolves its scaffolding to fit its phases without taxing
  ferrodec.
- ferrodec evolves its scaffolding without breaking pfloat.
- The scaffolding pattern is itself documented in this ADR, so a
  future plant-flag project (`interval-1788`, `tai-time`) inherits
  the convention by reading this ADR.

**Costs:**

- Two copies of substantially-similar scaffolding live in two repos.
  When a clippy allowance or a CI improvement applies to both, it
  has to be applied twice. The cost is paid in minutes per change
  at the volume both projects are likely to see.
- A change to the ADR template that should propagate to both will
  be on the next person to notice. Acceptable; the ADR template is
  near-stable already.
- New plant-flag projects copy from whichever earlier project
  carries the most relevant scaffolding shape; there is no canonical
  source. The user's `~/.claude/CLAUDE.md` carries the cross-project
  conventions, which is the canonical-but-not-binding answer.

## Related

- DESIGN.md, "Verification" section.
- ferrodec, the source of the imported scaffolding.
- This ADR is intended to be the model for future plant-flag
  projects' scaffolding decisions.
