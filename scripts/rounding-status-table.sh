#!/bin/sh
# scripts/rounding-status-table.sh
#
# Generate (and assert) the per-function rounding-status table that
# pfloat publishes at `docs/rounding-status.md`. The verification
# oracle's per-function status records under `tests/oracle/status/` are
# the single source of truth; this script renders them to markdown.
#
# Modes:
#   --emit   Print the full `docs/rounding-status.md` to stdout.
#   --check  Diff the rendered file against the committed
#            `docs/rounding-status.md`. Exit 0 on match, 1 on drift
#            (with a unified diff on stderr). Run in CI to keep the
#            published table from rotting behind the status records.
#
# Each status record carries, per function, the correct-rounding
# verdict in each of the five IEEE 754-2019 rounding modes. Two verdicts
# occur: `correctly-rounded` (certified across the exhaustive binary32
# grid) and `unswept` (the grid's bf-to-f32 bridge carries NearestEven
# only, so the directed mode is certified instead by the five-mode
# differential lanes and the cross-check sweep). The legend in the
# emitted file records this distinction; nothing here asserts more than
# the records do.
#
# POSIX `sh`; no bash features (verified under dash and the bash 3.2
# macOS ships as /bin/sh). Avoids `case` inside `$(...)` for the same
# reason `feature-union-check.sh` does. No Rust toolchain required.
#
# Run from any directory; the script discovers the repo root from its
# own location.

set -eu

LC_ALL=C
export LC_ALL

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
cd "$repo_root"

status_dir=tests/oracle/status
doc=docs/rounding-status.md

usage() {
    printf 'usage: %s {--emit | --check}\n' "$0" >&2
    exit 2
}

# --- Render the table body (sorted markdown rows) ---------------------
# One awk per file extracts the fields; a leading sort key (lowercased
# function name plus zero-padded order) groups parametric Bessel orders
# numerically within their family. The second awk maps the verdicts to
# `CR` / `CR(d)` and formats the markdown row.
table_rows() {
    for f in "$status_dir"/*.toml; do
        awk -F' *= *' '
            /^function /  { gsub(/"/, "", $2); fn  = $2 }
            /^order /     { gsub(/"/, "", $2); ord = $2 }
            /^oracle /    { gsub(/"/, "", $2); orc = $2 }
            /^NE /        { gsub(/"/, "", $2); ne  = $2 }
            /^NA /        { gsub(/"/, "", $2); na  = $2 }
            /^TZ /        { gsub(/"/, "", $2); tz  = $2 }
            /^TP /        { gsub(/"/, "", $2); tp  = $2 }
            /^TN /        { gsub(/"/, "", $2); tn  = $2 }
            END {
                key = tolower(fn) sprintf("%05d", ord + 0)
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
                    key, fn, ord, orc, ne, na, tz, tp, tn
            }
        ' "$f"
    done | sort | cut -f2- | awk -F'\t' '
        function disp(s) {
            return s == "correctly-rounded" ? "CR" \
                 : (s == "unswept" ? "CR(d)" : s)
        }
        {
            fn = $1; ord = $2; orc = $3
            ne = $4; na = $5; tz = $6; tp = $7; tn = $8
            name = (ord == "" ? fn : fn " (n=" ord ")")
            printf "| %s | %s | %s | %s | %s | %s | %s |\n", \
                name, orc, disp(ne), disp(na), disp(tz), disp(tp), disp(tn)
        }
    '
}

# --- Worst-case correctness numbers for the summary line --------------
# Max worst_ulp and the sums of mismatch/panic counts across all rows.
# Generated from the data so a future regression changes the published
# numbers (and trips `--check`) rather than hiding under stale prose.
worst_stats() {
    for f in "$status_dir"/*.toml; do
        awk -F' *= *' '
            /^worst_ulp /      { w = $2 + 0 }
            /^mismatch_count / { m = $2 + 0 }
            /^panic_count /    { p = $2 + 0 }
            END { printf "%d %d %d\n", w, m, p }
        ' "$f"
    done | awk '
        { if ($1 > maxw) maxw = $1; sm += $2; sp += $3 }
        END { printf "%d %d %d\n", maxw, sm, sp }
    '
}

emit() {
    nfns=$(ls "$status_dir"/*.toml 2>/dev/null | wc -l | tr -d ' ')
    set -- $(worst_stats)
    maxw=$1; summ=$2; span=$3
    cat <<MARKDOWN
# pfloat rounding status

Per-function correct-rounding status across the five IEEE 754-2019
rounding modes (NE NearestEven, NA NearestAway, TZ TowardZero, TP
TowardPositive, TN TowardNegative), for the $nfns functions the
verification oracle tracks. Generated from the status records under
\`tests/oracle/status/\` by \`scripts/rounding-status-table.sh\` and
checked in CI; the records are the single source of truth.

Across all rows the worst observed error is $maxw ULP, with $summ
mismatches and $span panics over the sampled input grids.

Legend:

- \`CR\`: correctly rounded, certified across the exhaustive binary32
  input grid. Every binary32 value is computed at high working
  precision and rounded to the target, then compared bit for bit
  against the oracle.
- \`CR(d)\`: correctly rounded, certified by the five-mode differential
  lanes against MPFR (with the synthesized NearestAway oracle, which
  MPFR lacks a primitive for) and reconfirmed by the per-release
  cross-check sweep (ADR-0049). The exhaustive binary32 oracle's
  bf-to-f32 bridge carries NearestEven only, so it does not sweep this
  directed mode; the guarantee for these cells rests on the lanes.

The oracle column names the primary rigorous backend: MPFR where it has
a primitive for the function, Arb otherwise.

| Function | Oracle | NE | NA | TZ | TP | TN |
| --- | --- | --- | --- | --- | --- | --- |
$(table_rows)
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
    printf 'rounding-status-table: %s drift detected\n' "$doc" >&2
    printf '  committed (left) vs script-emitted (right):\n' >&2
    diff -u "$doc" "$emit_tmp" >&2 || true
    exit 1
}

case "${1:-}" in
    --emit)  emit ;;
    --check) check ;;
    *)       usage ;;
esac
