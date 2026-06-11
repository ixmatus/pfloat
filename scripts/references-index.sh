#!/bin/sh
# scripts/references-index.sh
#
# Render (and assert) `docs/references/INDEX.md`, and validate every
# registry entry against the schema defined in ADR-0094 and
# `docs/references/README.md`. The index preamble is owned by this
# script; the list is generated from the entry files themselves, so the
# index cannot rot behind the directory as entries land. Each
# `docs/references/<slug>.md` is the single source of truth for its own
# metadata.
#
# This mirrors `scripts/adr-index.sh`, which likewise owns its whole
# output (so `--emit > file` is safe; the script never reads the file
# it writes).
#
# Modes:
#   --emit   Print the full `docs/references/INDEX.md` to stdout.
#   --check  Validate every entry (frontmatter keys, enums, consumer
#            paths, fixity hashes, body headings), then diff the
#            rendered index against the committed INDEX.md. Exit 0 on
#            success, 1 on any violation or drift. Run in CI; no
#            network, no Rust toolchain.
#
# POSIX `sh`; no bash features. Run from any directory; the script
# discovers the repo root from its own location.

set -eu

LC_ALL=C
export LC_ALL

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
cd "$repo_root"

ref_dir=docs/references
doc=$ref_dir/INDEX.md

usage() {
    printf 'usage: %s {--emit | --check}\n' "$0" >&2
    exit 2
}

# Files in docs/references/ that are not registry entries.
is_exempt() {
    case "$(basename "$1")" in
        README.md|TEMPLATE.md|INDEX.md|coverage-gaps.md) return 0 ;;
        *) return 1 ;;
    esac
}

entries() {
    for f in "$ref_dir"/*.md; do
        is_exempt "$f" && continue
        printf '%s\n' "$f"
    done
}

# Extract a scalar frontmatter value (first match wins).
fm() {
    awk -v key="$2" '
        NR == 1 && $0 != "---" { exit }
        NR > 1 && $0 == "---" { exit }
        $0 ~ "^" key ":" {
            sub("^" key ":[[:space:]]*", "")
            print
            exit
        }' "$1"
}

# --- Per-entry schema validation -------------------------------------
# Frontmatter keys, in the exact required order. `consumers` is a list:
# its items are the `  - path` lines between it and `verification`.
KEYS='slug category citation edition canonical_url document_number doi
archived_url archive_date retrieval_date sha256 vendored_path license
vendor_status rot_risk provenance_class consumers verification'

validate_entry() {
    f=$1
    base=$(basename "$f")
    stem=${base%.md}
    errs=0

    err() {
        printf 'references-index: %s: %s\n' "$base" "$1" >&2
        errs=$((errs + 1))
    }

    # Key sequence must match the canonical order exactly.
    want=$(printf '%s\n' $KEYS)
    got=$(awk '
        NR == 1 && $0 != "---" { exit 1 }
        NR > 1 && $0 == "---" { exit }
        NR > 1 && /^[a-z0-9_]+:/ { sub(":.*$", ""); print }' "$f") || {
        err 'missing frontmatter (file must start with ---)'
        printf '%d\n' "$errs"
        return
    }
    if [ "$want" != "$got" ]; then
        err 'frontmatter keys missing, unknown, or out of order'
    fi

    # Every scalar key non-empty.
    for k in $KEYS; do
        [ "$k" = consumers ] && continue
        v=$(fm "$f" "$k")
        if [ -z "$v" ]; then
            err "empty value for '$k'"
        fi
    done

    slug=$(fm "$f" slug)
    [ "$slug" = "$stem" ] || err "slug '$slug' != filename stem '$stem'"

    category=$(fm "$f" category)
    case "$category" in
        standard|paper|book|web|corpus|software) ;;
        *) err "bad category '$category'" ;;
    esac

    vendor_status=$(fm "$f" vendor_status)
    case "$vendor_status" in
        vendored|pointer-only|legally-cannot|derived-subset) ;;
        *) err "bad vendor_status '$vendor_status'" ;;
    esac

    rot_risk=$(fm "$f" rot_risk)
    case "$rot_risk" in
        died-once|single-maintainer|community-run|academic-personal|stable-publisher|ephemeral) ;;
        *) err "bad rot_risk '$rot_risk'" ;;
    esac

    provenance_class=$(fm "$f" provenance_class)
    case "$provenance_class" in
        primary|oracle|lineage|contextual) ;;
        *) err "bad provenance_class '$provenance_class'" ;;
    esac

    # Consumer paths: at least one, each must exist in the repo
    # (fix forward on renames; never delete the line).
    consumers=$(awk '
        NR == 1 && $0 != "---" { exit }
        NR > 1 && $0 == "---" { exit }
        /^consumers:/ { on = 1; next }
        on && /^  - / { sub(/^  - /, ""); print; next }
        on { exit }' "$f")
    if [ -z "$consumers" ]; then
        err 'consumers list is empty'
    else
        printf '%s\n' "$consumers" | while IFS= read -r p; do
            if [ ! -e "$p" ]; then
                printf 'references-index: %s: consumer path does not exist: %s\n' "$base" "$p" >&2
            fi
        done
        missing=$(printf '%s\n' "$consumers" | { c=0; while IFS= read -r p; do [ -e "$p" ] || c=$((c + 1)); done; printf '%s' "$c"; })
        [ "$missing" -eq 0 ] || errs=$((errs + missing))
    fi

    # Fixity: vendored and derived-subset entries pin an in-tree
    # artifact by hash. Regenerating the artifact must update the hash.
    sha256=$(fm "$f" sha256)
    vendored_path=$(fm "$f" vendored_path)
    case "$vendor_status" in
        vendored|derived-subset)
            [ "$sha256" != none ] || err "vendor_status $vendor_status requires sha256"
            [ "$vendored_path" != none ] || err "vendor_status $vendor_status requires vendored_path"
            if [ "$vendored_path" != none ] && [ -f "$vendored_path" ] && [ "$sha256" != none ]; then
                actual=$(shasum -a 256 "$vendored_path" | awk '{print $1}')
                [ "$actual" = "$sha256" ] || err "sha256 mismatch for $vendored_path (recorded $sha256, actual $actual)"
            fi
            ;;
    esac

    # Rot prone sources must carry a verified archive snapshot.
    archived_url=$(fm "$f" archived_url)
    case "$rot_risk" in
        died-once|single-maintainer|academic-personal|ephemeral)
            [ "$archived_url" != none ] || err "rot_risk $rot_risk requires archived_url"
            ;;
    esac

    # Required body headings.
    for h in '## Why this source' '## What it grounds' '## Alternatives'; do
        grep -q "^$h\$" "$f" || err "missing body heading '$h'"
    done
    if [ "$category" = corpus ]; then
        grep -q '^## Coverage gaps$' "$f" || err "corpus entry missing '## Coverage gaps'"
    fi

    printf '%d\n' "$errs"
}

# --- Render the index body (one row per entry, sorted by slug) --------
index_rows() {
    entries | sort | while IFS= read -r f; do
        base=$(basename "$f")
        slug=${base%.md}
        category=$(fm "$f" category)
        provenance_class=$(fm "$f" provenance_class)
        rot_risk=$(fm "$f" rot_risk)
        printf -- '- [%s](%s) (%s, %s, %s)\n' \
            "$slug" "$base" "$category" "$provenance_class" "$rot_risk"
    done
}

emit() {
    # Preamble, owned here (quoted heredoc: backticks are literal).
    cat <<'PREAMBLE'
# Reference registry index

One line per registry entry; the entry file is the source of truth for
everything else. See `README.md` in this directory for the schema and
the rules, and ADR-0094 for the decision record.

PREAMBLE
    # Generated index (unquoted heredoc: backticks escaped, $(...) runs).
    cat <<MARKDOWN
This list is generated by \`scripts/references-index.sh\` and checked in
CI (\`scripts/references-index.sh --check\`). Edit the script, not this
file by hand. Row format: slug (category, provenance class, rot risk).

$(index_rows)
MARKDOWN
}

check() {
    total=0
    for f in $(entries); do
        e=$(validate_entry "$f")
        total=$((total + e))
    done
    if [ "$total" -gt 0 ]; then
        printf 'references-index: %d schema violation(s)\n' "$total" >&2
        exit 1
    fi
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
    printf 'references-index: %s drift detected\n' "$doc" >&2
    printf '  committed (left) vs script-emitted (right):\n' >&2
    diff -u "$doc" "$emit_tmp" >&2 || true
    exit 1
}

case "${1:-}" in
    --emit)  emit ;;
    --check) check ;;
    *)       usage ;;
esac
