#!/bin/sh
# scripts/ball-enclosure-status.sh
#
# Assert that the published pfloat-ball enclosure-accuracy posture
# (`docs/ball-enclosure-status.md`) declares a posture for every public
# ball operation, and lists no operation that is not on the surface.
#
# Unlike `scripts/rounding-status-table.sh`, this posture is authored, not
# generated: the enclosure shape and the tightest/accurate class are design
# declarations the measured `differential_arb` tightness lane substantiates,
# not values read off a status record. So this script does NOT regenerate
# the doc; it gates its COMPLETENESS, the same anti-rot role
# `scripts/feature-union-check.sh` plays for the CI feature union.
#
# The source of truth for the operation set is the `pub fn` surface of
# `pfloat-ball/src/arith.rs` and `pfloat-ball/src/elem.rs` (every public fn
# there is a ball operation). The declared set is the backtick-quoted
# identifiers in the table rows of the posture doc (the lines beginning
# `|`). The two must be equal: a new operation added without a posture row,
# or a posture row left behind after a rename, fails this check.
#
# Modes:
#   (no args) / --check  Run the assertion. Exit 0 on agreement, 1 on
#                        drift (with a diagnostic on stderr). The CI mode.
#   --show               Print the two derived sets, then run the check.
#
# POSIX `sh`; no bash features. Runs from any directory; discovers the repo
# root from its own location. No Rust toolchain required, so the cheap
# `conformance` CI job runs it without building.

set -eu

LC_ALL=C
export LC_ALL

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
cd "$repo_root"

doc=docs/ball-enclosure-status.md
arith=pfloat-ball/src/arith.rs
elem=pfloat-ball/src/elem.rs

for f in "$doc" "$arith" "$elem"; do
    if [ ! -f "$f" ]; then
        printf 'ball-enclosure-status: %s not found\n' "$f" >&2
        exit 1
    fi
done

# --- Source of truth: the public ball ops -----------------------------
# Every `pub fn <name>` in arith.rs and elem.rs, sorted and deduped.
surface=$(grep -hoE 'pub fn [a-z0-9_]+' "$arith" "$elem" \
    | sed -E 's/pub fn //' | sort -u)

# --- Declared set: backtick-quoted idents in the doc's table rows ------
# Table rows begin with `|`; only operation names are backtick-quoted in
# them (the shape and posture cells are plain prose), so this selects the
# declared operations and nothing else.
declared=$(grep '^|' "$doc" \
    | grep -oE '`[a-z0-9_]+`' | tr -d '`' | sort -u)

if [ -z "$surface" ]; then
    printf 'ball-enclosure-status: no `pub fn` found in %s / %s\n' "$arith" "$elem" >&2
    exit 1
fi
if [ -z "$declared" ]; then
    printf 'ball-enclosure-status: no backtick-quoted ops in the table rows of %s\n' "$doc" >&2
    exit 1
fi

case "${1:-}" in
    --show)
        printf 'surface  : %s\n' "$(printf '%s' "$surface" | tr '\n' ' ')"
        printf 'declared : %s\n' "$(printf '%s' "$declared" | tr '\n' ' ')"
        ;;
    --check|"") : ;;
    *)
        printf 'usage: %s [--check | --show]\n' "$0" >&2
        exit 2
        ;;
esac

if [ "$surface" != "$declared" ]; then
    printf 'ball-enclosure-status: the posture doc drifted from the public ball surface\n' >&2
    surf_tmp=$(mktemp)
    decl_tmp=$(mktemp)
    trap 'rm -f "$surf_tmp" "$decl_tmp"' EXIT INT TERM
    printf '%s\n' "$surface" > "$surf_tmp"
    printf '%s\n' "$declared" > "$decl_tmp"
    missing=$(comm -23 "$surf_tmp" "$decl_tmp" | tr '\n' ' ')
    extra=$(comm -13 "$surf_tmp" "$decl_tmp" | tr '\n' ' ')
    [ -n "$missing" ] && printf '  public ops with NO posture row: %s\n' "$missing" >&2
    [ -n "$extra" ]   && printf '  posture rows for a non-existent op: %s\n' "$extra" >&2
    printf '  add the operation to the table in %s (or fix the name).\n' "$doc" >&2
    exit 1
fi

exit 0
