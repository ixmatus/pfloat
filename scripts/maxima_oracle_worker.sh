#!/usr/bin/env nix-shell
#!nix-shell -i bash -p maxima python3
#
# Maxima oracle worker (ADR-0035 Tier 6 sampling oracle).
#
# This is the wrapper: nix-shell provides the Maxima interpreter
# without requiring a system install. The wrapper forwards "$@" to
# the inner script that does the per-request dispatch via
# `maxima --very-quiet --batch-string=...`.
#
# Per-request cost is dominated by Maxima startup (~500ms-1s);
# Tier 6 use mode is sampling (hand-derived corpus + tie-breakers
# + N-sample per release), not full f32 sweep. ADR-0035 records
# the use-mode choice.
#
# Function coverage caveats (probed at slice p1.10):
#
# - `bessel_i(n, x)` at very small subnormals (|x| < ~2^-200) can
#   trigger Maxima's "Exceeded maximum allowed fpprec" because its
#   internal hypergeometric falls through to a symbolic form. The
#   worker treats these as INC for the affected inputs rather than
#   ERR; the sampling layer's coverage gap on the smallest
#   subnormals is acceptable since Arb + mpmath two-oracle
#   agreement already covers that class.
# - Maxima does not have a direct logarithmic-integral primitive;
#   `li[1](x)` is the polylogarithm, not the logarithmic integral.
#   The worker composes via `li_via_ei(x) = ei(log(x))` for the
#   `li` FnId. Maxima has `gamma_incomplete` and `expintegral_ei`
#   suitable for this.
# - `airy_dai(x)` and `airy_dbi(x)` are the Airy derivatives in
#   Maxima's naming.
#
# Protocol: same as the Arb / mpmath workers. Request line on
# stdin: `<fn_id> <order_or_dash> <input_bits_hex> <mode>`.
# Response on stdout: `OK <f32_bits_hex>` | `INC` | `ERR <msg>`.

set -euo pipefail

# Locate the helper script: same directory as this wrapper.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HELPER="${SCRIPT_DIR}/maxima_oracle_worker.py"

if [ ! -x "$(command -v python3)" ]; then
    echo "ERR python3 not found in nix-shell environment" >&2
    exit 1
fi

if [ ! -f "$HELPER" ]; then
    echo "ERR helper script not found at $HELPER" >&2
    exit 1
fi

# The helper handles the protocol; maxima is invoked per request
# via subprocess. We forward "$@" so the helper can take optional
# flags (e.g. --once for single-request mode in slice p1.11+
# tests).
exec python3 "$HELPER" "$@"
