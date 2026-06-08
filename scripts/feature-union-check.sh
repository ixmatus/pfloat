#!/bin/sh
# scripts/feature-union-check.sh
#
# Assert that CI's full-feature-union list cannot silently drift from
# the crate's actual feature surface.
#
# Two facts are gated:
#
#   1. The test-matrix union entry and the clippy job's `--features=`
#      list in `.github/workflows/ci.yml` are byte-identical. They must
#      stay in lockstep: clippy lints the surface the matrix tests, and
#      a kernel feature present in one but not the other is a hole.
#
#   2. Every `[features]` key in `Cargo.toml` either appears in that
#      union list or is on the explicit exclusion allowlist below. A
#      new kernel feature then fails this check until it is added to the
#      union (or deliberately excluded), so the per-push gate cannot
#      drift behind the surface it is meant to cover.
#
# This is the runnable form of a lesson that used to live in reviewer
# memory: the per-push matrix once ran a 6-combo subset that did not
# enable the high-tier kernels, and defects in those kernels' tests
# sat unsurfaced because no combo turned the feature on. Entry 7 (the
# union) closed that, but the union is a hardcoded list; without this
# guard, the next kernel feature added to `Cargo.toml` but forgotten in
# CI would silently reopen the hole. ADR-0053.
#
# Modes:
#   (no args) / --check  Run the assertions. Exit 0 on agreement, 1 on
#                        drift (with a diagnostic on stderr). The CI
#                        mode.
#   --show               Print the three derived lists (matrix union,
#                        clippy union, expected-from-Cargo.toml) for
#                        debugging, then run the check.
#
# POSIX `sh`; no bash features (verified under both dash and the
# bash 3.2 that macOS ships as /bin/sh). Runs anywhere `sh`, `awk`,
# `grep`, `comm`, `sort`, and `mktemp` exist. No Rust toolchain
# required, so the cheap `conformance` CI job runs it without building.
#
# Run from any directory; the script discovers the repo root from its
# own location.

set -eu

# Byte-wise locale so `sort` is deterministic across machines.
LC_ALL=C
export LC_ALL

# Discover the repo root from this script's location.
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
cd "$repo_root"

ci=.github/workflows/ci.yml
cargo=Cargo.toml

# Features that legitimately do NOT belong in the per-push union list:
#   default            the default *set*, not a feature; its members
#                      (std, fmt, big) are listed explicitly instead.
#   alloc              transitively pulled by `std`/`big`; never named
#                      on its own in the union line.
#   kani               compiles proof harnesses; a separate manual job.
#   ziv-instrumented   test/cross-check instrumentation only.
#   differential-mpfr  differential/oracle tiers; their own CI job with
#   differential-arb   the gmp-mpfr-sys / Arb system dependencies.
# Adding a new name here is a deliberate edit, reviewed like any other.
exclude="default alloc kani ziv-instrumented differential-mpfr differential-arb"

# --- Extract the two CI union lists -----------------------------------
# Matrix union entry: the sole matrix line carrying `--features=` with
# no `--no-default-features` prefix (entries 3-6 carry both, so the
# anchored `- "--features=` pattern selects entry 7 alone).
matrix_union=$(grep -E '^[[:space:]]*-[[:space:]]*"--features=' "$ci" \
    | sed -E 's/.*--features=([^"]*)".*/\1/')

# Clippy union: the `cargo clippy ... --features=<list> ...` line for the
# root `pfloat` crate. The workspace members (`pfloat-libm`, `pfloat-ball`)
# have their own clippy lines carrying their own features (`differential-mpfr`,
# `big`, ...); those are different crates' surfaces, not the root's
# kernel-feature union, so exclude them here.
clippy_union=$(grep -E 'cargo clippy .*--features=' "$ci" \
    | grep -v -e 'pfloat-libm' -e 'pfloat-ball' -e 'pfloat-complex' \
    | sed -E 's/.*--features=([^ ]*).*/\1/')

if [ -z "$matrix_union" ]; then
    printf 'feature-union-check: no `- "--features=..."` matrix entry found in %s\n' "$ci" >&2
    exit 1
fi
if [ -z "$clippy_union" ]; then
    printf 'feature-union-check: no `cargo clippy ... --features=` line found in %s\n' "$ci" >&2
    exit 1
fi

# --- Derive the expected union from Cargo.toml ------------------------
# Every `[features]` key, minus the exclusion allowlist, sorted.
all_features=$(awk '
    /^\[features\]/ { in_f = 1; next }
    in_f && /^\[/   { in_f = 0 }
    in_f && /^[A-Za-z0-9_-]+[[:space:]]*=/ { print $1 }
' "$cargo")

# Drop the exclusion allowlist from the feature set, then sort. Done
# with `grep -vxF` against a one-per-line allowlist file rather than a
# `case` inside this command substitution: bash 3.2 (macOS /bin/sh)
# mis-parses `case ... esac` nested in `$(...)`, so the filter stays
# out of the substitution.
exclude_file=$(mktemp)
printf '%s\n' $exclude > "$exclude_file"
expected=$(printf '%s\n' "$all_features" | grep -vxF -f "$exclude_file" | sort)
rm -f "$exclude_file"

actual=$(printf '%s' "$matrix_union" | tr ',' '\n' | sort)

# --- --show: print the derived lists ----------------------------------
case "${1:-}" in
    --show)
        printf 'matrix union : %s\n' "$matrix_union"
        printf 'clippy union : %s\n' "$clippy_union"
        printf 'expected set : %s\n' "$(printf '%s' "$expected" | tr '\n' ' ')"
        ;;
    --check|"") : ;;
    *)
        printf 'usage: %s [--check | --show]\n' "$0" >&2
        exit 2
        ;;
esac

# --- Assertion 1: matrix union == clippy union ------------------------
if [ "$matrix_union" != "$clippy_union" ]; then
    printf 'feature-union-check: matrix and clippy feature lists differ\n' >&2
    printf '  matrix: %s\n' "$matrix_union" >&2
    printf '  clippy: %s\n' "$clippy_union" >&2
    printf '  the test matrix and the clippy gate must lint/test the same surface.\n' >&2
    exit 1
fi

# --- Assertion 2: Cargo.toml surface == CI union (modulo allowlist) ---
if [ "$expected" != "$actual" ]; then
    printf 'feature-union-check: CI union drifted from the Cargo.toml feature surface\n' >&2
    exp_tmp=$(mktemp)
    act_tmp=$(mktemp)
    trap 'rm -f "$exp_tmp" "$act_tmp"' EXIT INT TERM
    printf '%s\n' "$expected" > "$exp_tmp"
    printf '%s\n' "$actual" > "$act_tmp"
    # `comm` on the two sorted lists: -23 = only-in-expected (missing
    # from CI), -13 = only-in-actual (named in CI but not a feature).
    missing=$(comm -23 "$exp_tmp" "$act_tmp" | tr '\n' ' ')
    extra=$(comm -13 "$exp_tmp" "$act_tmp" | tr '\n' ' ')
    [ -n "$missing" ] && printf '  in Cargo.toml but NOT in the CI union: %s\n' "$missing" >&2
    [ -n "$extra" ]   && printf '  in the CI union but NOT a Cargo.toml feature: %s\n' "$extra" >&2
    printf '  add the feature to the union in %s, or to the exclusion allowlist in this script.\n' "$ci" >&2
    exit 1
fi

exit 0
