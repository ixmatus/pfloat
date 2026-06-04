#!/bin/sh
# scripts/adr-index.sh
#
# Render (and assert) `docs/decisions/README.md`. The directory preamble
# (Conventions, Writing a new ADR) is owned by this script; the `## Index`
# list is generated from the ADR files themselves, so the index cannot
# rot behind the directory as ADRs land. Each `docs/decisions/NNNN-slug.md`
# is the single source of truth for its own title and status.
#
# To change the preamble, edit the heredoc below, not the generated file.
# This mirrors `scripts/rounding-status-table.sh`, which likewise owns its
# whole output (so `--emit > file` is safe; the script never reads the
# file it writes).
#
# Modes:
#   --emit   Print the full `docs/decisions/README.md` to stdout.
#   --check  Diff the rendered file against the committed
#            `docs/decisions/README.md`. Exit 0 on match, 1 on drift
#            (with a unified diff on stderr). Run in CI to keep the
#            index from rotting behind the ADR directory.
#
# POSIX `sh`; no bash features. No Rust toolchain required. Run from any
# directory; the script discovers the repo root from its own location.

set -eu

LC_ALL=C
export LC_ALL

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
cd "$repo_root"

doc=docs/decisions/README.md
adr_dir=docs/decisions

usage() {
    printf 'usage: %s {--emit | --check}\n' "$0" >&2
    exit 2
}

# --- Render the index body (one row per ADR, sorted by number) --------
# Each ADR's title is its first `# ` heading (minus the marker); the
# filename is the link target. The status line comes in two historical
# shapes (`- **Status**: accepted` and `Status: Accepted (date)`); both
# are matched and normalised to the leading keyword, lowercased, so dates
# and slice notes in the file do not leak into the index. The shell glob
# sorts numerically under LC_ALL=C.
index_rows() {
    for f in "$adr_dir"/[0-9][0-9][0-9][0-9]-*.md; do
        base=$(basename "$f")
        title=$(awk '/^# / { sub(/^# +/, ""); print; exit }' "$f")
        status=$(awk '
            /[Ss]tatus(\*\*)?:/ {
                sub(/^.*[Ss]tatus(\*\*)?:[[:space:]]*/, "")
                kw = tolower($0)
                sub(/[^a-z-].*$/, "", kw)
                print kw
                exit
            }' "$f")
        printf -- '- [%s](%s) (%s)\n' "$title" "$base" "$status"
    done
}

emit() {
    # Preamble, owned here (quoted heredoc: backticks are literal).
    cat <<'PREAMBLE'
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

PREAMBLE
    # Generated index (unquoted heredoc: backticks escaped, $(...) runs).
    cat <<MARKDOWN
## Index

This list is generated from the ADR files by \`scripts/adr-index.sh\`
and checked in CI (\`scripts/adr-index.sh --check\`); each ADR file is
the source of truth for its own title and status. Edit the script (the
preamble and this note live in its heredocs), not this file by hand.

$(index_rows)
MARKDOWN
}

check() {
    if [ ! -f "$doc" ]; then
        printf 'error: %s not found; run `%s --emit > %s`\n' "$doc" "$0" "$doc" >&2
        exit 1
    fi
    emit_tmp=$(mktemp)
    trap 'rm -f "$emit_tmp"' EXIT INT TERM
    emit > "$emit_tmp"
    if diff -u "$doc" "$emit_tmp" >/dev/null 2>&1; then
        exit 0
    fi
    printf 'adr-index: %s drift detected\n' "$doc" >&2
    printf '  committed (left) vs script-emitted (right):\n' >&2
    diff -u "$doc" "$emit_tmp" >&2 || true
    exit 1
}

case "${1:-}" in
    --emit)  emit ;;
    --check) check ;;
    *)       usage ;;
esac
