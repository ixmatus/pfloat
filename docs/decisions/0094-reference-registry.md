# ADR-0094: Per-source reference registry at docs/references/

- **Status**: accepted
- **Date**: 2026-06-11

## Context

pfloat implements external standards and published algorithms, and its
README promises that implementations derive from primary sources. The
citations backing that promise were scattered: a flat catalog at
`docs/references.md`, one deep provenance document for the hard to round
corpus, per-module doc comments, and some thirty five ADRs. Nothing
recorded archive snapshots, retrieval dates, licenses (outside the corpus
document), rot risk, or fixity hashes, so the registry of what pfloat
stands on could rot silently: personal academic pages disappear, vendor
URLs move, and a citation that can no longer be checked is a provenance
claim on faith.

A cross project convention (the manuals program; pfloat slice tracked as
bead `pf-ysut`) now requires each repository to carry a self contained
per-source registry that downstream synthesis repositories copy verbatim.
The crate level CLAUDE.md records the accretion ritual; this ADR records
the schema and the decisions inside it.

## Decision

One markdown file per external source at `docs/references/<slug>.md`,
with eighteen required frontmatter keys (schema table in
`docs/references/README.md`): identity (`slug`, `category`, `citation`,
`edition`, `document_number`, `doi`), location and rot insurance
(`canonical_url`, `archived_url`, `archive_date`, `retrieval_date`,
`rot_risk`), legal posture (`license`, `vendor_status`, `sha256`,
`vendored_path`), and honesty linkage (`provenance_class`, `consumers`,
`verification`). `INDEX.md` is generated and checked by
`scripts/references-index.sh` (the `adr-index.sh` pattern); CI runs
`--check` in the conformance job, offline. Wayback saves happen at mining
time from the dev machine, and an archive date is recorded only after the
snapshot has been opened and verified.

Four provenance classes keep the registry honest about what the tree
actually does with each source: `primary` (an implementation derives from
it), `oracle` (a verification lane compares against it), `lineage`
(acknowledged through a direct source), and `contextual` (named for
context, proxy, or gap documentation). Two entries are deliberately
catalogued because pfloat does NOT use them: Berkeley TestFloat (its
fixed format vectors do not apply at arbitrary precision, which is the
structural reason third party conformance vectors do not exist for this
surface) and the Muller et al Handbook (a bibliographic anchor the corpus
provenance document explicitly does not transcribe from). They are the
canonical examples of `contextual`; class inflation toward `primary` is
this registry's named failure mode.

`vendor_status` extends the cross project three value enum (`vendored`,
`pointer-only`, `legally-cannot`) with a fourth value, `derived-subset`:
a source whose data was transcribed or sampled into a derived in tree
artifact without redistributing upstream bytes verbatim. The CORE-MATH
hard to round corpus is the motivating case: inputs are transcribed,
outputs are independently recomputed, and the upstream `.wc` files are
not vendored. A `derived-subset` entry pins the derived artifact by
`sha256`, so regenerating the artifact forces a registry update in the
same change. The owner chose widening the enum over recording
`pointer-only` with prose; the divergence is flagged to the master
program (`smil-27za`) for estate wide adoption.

Theorems get no entries (the Fundamental Theorem of Interval Arithmetic
is stated canonically in tree at `pfloat-ball/src/spec.rs`); toolchain
dependencies carry provenance in their manifests; and sources the code
does not actually cite are not added for prestige. The function major
catalog at `docs/references.md` remains, with its bibliography section
replaced by links into this registry.

## Consequences

- Every load bearing citation gains an archived snapshot, a license
  record, and a machine checked consumers list; a future maintainer can
  verify provenance without the original authors. The registry is the
  durable artifact the README's provenance paragraph points at.
- The consumer path check turns module renames into registry edits. That
  is the cost of auditability; the rule is fix forward (update the path),
  never delete the consumer line to silence CI.
- The `derived-subset` fixity pin makes corpus regeneration a two file
  change (data plus registry). Deliberate: regenerating a corpus is a
  provenance event.
- Offline CI cannot detect a rotted or stub snapshot. The mitigation is
  procedural, not mechanical: verify at save time, re-verify rot prone
  snapshots roughly annually.
- Downstream copies can go stale against this registry; the copy contract
  (importer rewrites `consumers` and `verification`, keeps all other
  fields) is prose, not tooling, until the master program builds the
  drift check sweep.

## Related

- Beads: `pf-ysut` (program), master brief `smil-27za`.
- Other ADRs: ADR-0083 (corpus seeding; the corpus this schema's
  `derived-subset` value was shaped by), ADR-0087 (the prior
  docs-posture CI gate this one follows the pattern of).
- `docs/references/README.md` (schema of record),
  `scripts/references-index.sh` (validator),
  `docs/lefevre-muller-corpus-provenance.md` (narrative depth the
  corpus registry entry summarizes).
