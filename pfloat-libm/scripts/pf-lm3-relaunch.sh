#!/usr/bin/env bash
#
# Relaunch one pf-lm3 shard that failed or was reclaimed, on-demand, into
# the same run prefix. Adapts pfloat's pf-hcz4-relaunch.sh.
#
# Usage: pf-lm3-relaunch.sh <RUN_ID> <fn>_<K>of<M>
#   e.g. pf-lm3-relaunch.sh pf-lm3-20260601-120000 exp_3of16

set -euo pipefail

RUN_ID="${1:?usage: $0 <RUN_ID> <fn>_<K>of<M>}"
SHARD="${2:?usage: $0 <RUN_ID> <fn>_<K>of<M>}"

export MARKET="${MARKET:-on-demand}"   # force on-demand for the retry
export RUN_ID

exec "$(dirname "$0")/pf-lm3-launch.sh" --yes-i-will-pay --only "${SHARD}" --run-id "${RUN_ID}"
