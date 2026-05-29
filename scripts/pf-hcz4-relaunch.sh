#!/usr/bin/env bash
#
# pf-hcz4 single-shard relaunch (spot reclaim recovery). ADR-0049.
#
# Reclaim rate on spot is empirically ~10% per ferrodec. This script
# relaunches a single shard, optionally forcing on-demand market.
#
# Usage:
#   pf-hcz4-relaunch.sh <RUN_ID> <shard>
#   MARKET=on-demand pf-hcz4-relaunch.sh <RUN_ID> Yn:5
#
# Just a thin wrapper that re-invokes pf-hcz4-launch.sh with --only.

set -euo pipefail

RUN_ID="${1:?usage: $0 <RUN_ID> <shard>}"
SHARD="${2:?usage: $0 <RUN_ID> <shard>}"

export MARKET="${MARKET:-on-demand}"   # default to on-demand for reliability
export RUN_ID

# Delegate.
exec "$(dirname "$0")/pf-hcz4-launch.sh" --only "${SHARD}" --run-id "${RUN_ID}"
