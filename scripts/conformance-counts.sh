#!/bin/sh
# scripts/conformance-counts.sh
#
# Emit and assert pfloat's per-bucket conformance harness counts.
#
# Modes:
#   --emit   Print the canonical `## Conformance evidence` markdown
#            block to stdout. The same numbers live in README.md;
#            this script is the single source of truth.
#
#   --check  Diff the emitted block against README.md's
#            `## Conformance evidence` section. Exit 0 on match,
#            1 on drift (with a unified diff on stderr). Run in CI
#            to prevent the README counts from rotting.
#
# The buckets (Kani harnesses, differential lanes, fuzz targets,
# property files) are asserted independently. A per-bucket gate is
# the failure-mode the plan calls out: an aggregate floor would let
# one bucket silently shrink while another grew. The differential
# sweep size is extracted from ADR-0014 rather than hardcoded so
# the gate and the design record cannot drift apart.
#
# POSIX `sh` — no bash features. Runs anywhere `sh`, `awk`, `sed`,
# `grep`, `find`, `diff`, and `mktemp` exist (the standard Unix
# toolset). No Rust toolchain required, so CI can run it without
# building anything.
#
# Run from any directory; the script discovers the repo root from
# its own location.

set -eu

# Byte-wise locale so `sort` and `grep` are deterministic across
# machines (a French-locale CI sorts `é` differently from a C-locale
# laptop, which would surface as spurious drift).
LC_ALL=C
export LC_ALL

usage() {
    printf 'usage: %s {--emit | --check}\n' "$0" >&2
    exit 2
}

# Discover the repo root from this script's location.
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
cd "$repo_root"

# --- Bucket: Kani proof harnesses -------------------------------------
# Each `#[kani::proof]` attribute is one harness. `helpers.rs` is the
# one verify file without a proof attribute (utility module).
kani_proofs=$(grep -h '#\[kani::proof\]' src/verify/*.rs | wc -l | tr -d ' ')
kani_files=$(grep -l '#\[kani::proof\]' src/verify/*.rs | wc -l | tr -d ' ')

# --- Bucket: differential lanes ---------------------------------------
# One toplevel `tests/differential_*.rs` file per kernel or op. The
# shared harness module under `tests/differential/` is intentionally
# not counted (it is scaffolding, not a lane).
differential_files=$(ls tests/differential_*.rs 2>/dev/null | wc -l | tr -d ' ')

# Sweep sizes extracted from ADR-0014 (the design-record source of
# truth). Both lines are stable bolded list-items; the sed pulls the
# token immediately after the colon and before "inputs".
adr14=docs/decisions/0014-mpfr-differential-ci-gating.md
sweep_ci=$(grep -F '**CI default**' "$adr14" \
    | sed -n 's/^[^:]*: *\([^ ]*\) inputs.*/\1/p')
sweep_deep=$(grep -F 'Local deep' "$adr14" \
    | sed -n 's/^[^:]*: *\([^ ]*\) inputs.*/\1/p')

# --- Bucket: fuzz targets ---------------------------------------------
fuzz_targets=$(ls fuzz/fuzz_targets/*.rs 2>/dev/null | wc -l | tr -d ' ')

# --- Bucket: property tests -------------------------------------------
# `tests/property_*.rs` files. `.proptest-regressions` artifacts are
# excluded by the `*.rs` glob — they are generated regression seeds,
# not tests.
property_files=$(ls tests/property_*.rs 2>/dev/null | wc -l | tr -d ' ')

# Per-op listing: file basenames stripped of the `property_` prefix
# and `.rs` suffix, sorted, space-joined. Concrete enumeration is
# what the per-bucket discipline asks for — a count alone hides
# which op the regression is in.
property_ops=$(ls tests/property_*.rs 2>/dev/null \
    | sed -e 's|^tests/property_||' -e 's|\.rs$||' \
    | sort \
    | tr '\n' ' ' \
    | sed 's/ $//')

# --- Emit the markdown block ------------------------------------------
emit() {
    cat <<MARKDOWN
## Conformance evidence

Per-bucket counts of pfloat's verification harnesses. Each bucket
is asserted independently by \`scripts/conformance-counts.sh\` and
gated in CI; one bucket shrinking while another grows cannot hide
under an aggregate floor.

- **Kani proof harnesses:** ${kani_proofs} \`#[kani::proof]\`
  attributes across ${kani_files} files in \`src/verify/\`.
- **Differential lanes:** ${differential_files} \`tests/differential_*.rs\`
  files. CI sweep ${sweep_ci} inputs per (op × precision × rounding
  mode); \`PFLOAT_DEEP=1\` escalates to ${sweep_deep} (ADR-0014).
- **Fuzz targets:** ${fuzz_targets} targets under \`fuzz/fuzz_targets/\`.
- **Property tests:** ${property_files} \`tests/property_*.rs\` files
  (${property_ops}).
MARKDOWN
}

# --- Check mode: diff emitted vs README -------------------------------
check() {
    if [ ! -f README.md ]; then
        printf 'error: README.md not found at %s/README.md\n' "$repo_root" >&2
        exit 2
    fi
    emit_tmp=$(mktemp)
    readme_tmp=$(mktemp)
    # POSIX `trap` cleanup; the script exits via `exit` only, so
    # EXIT covers normal and error paths.
    trap 'rm -f "$emit_tmp" "$readme_tmp"' EXIT INT TERM

    # `strip_trailing_blanks` normalizes both blocks before diffing.
    # The README naturally carries a blank line before the next `##`
    # heading (markdown convention); `--emit` does not. Without this
    # step `--check` would always flag a one-line difference.
    strip_trailing_blanks() {
        awk '
            /./ { for (i = 0; i < blanks; i++) print ""; blanks = 0; print; next }
            /^$/ { blanks++ }
        ' "$1"
    }

    emit | strip_trailing_blanks /dev/stdin > "$emit_tmp"

    # Extract README's `## Conformance evidence` section: from the
    # heading line up to (but not including) the next `## ` heading.
    awk '
        /^## Conformance evidence$/ { in_block = 1 }
        in_block && /^## / && !/^## Conformance evidence$/ { in_block = 0 }
        in_block { print }
    ' README.md | strip_trailing_blanks /dev/stdin > "$readme_tmp"

    if [ ! -s "$readme_tmp" ]; then
        printf 'error: README.md has no "## Conformance evidence" section\n' >&2
        exit 1
    fi

    if diff -u "$readme_tmp" "$emit_tmp" >/dev/null 2>&1; then
        exit 0
    fi

    printf 'conformance-counts: README.md drift detected\n' >&2
    printf '  expected (left) vs script-emitted (right):\n' >&2
    diff -u "$readme_tmp" "$emit_tmp" >&2 || true
    exit 1
}

case "${1:-}" in
    --emit)  emit ;;
    --check) check ;;
    *)       usage ;;
esac
