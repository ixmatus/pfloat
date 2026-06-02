#!/usr/bin/env bash
#
# pf-lm3 sweep status poller. Reports running instances and S3 sentinels
# for a run; when nothing is running it syncs results and suggests the
# aggregate step. Adapts pfloat's pf-hcz4-status.sh.
#
# Usage: pf-lm3-status.sh <RUN_ID>
#   S3_BUCKET, AWS_REGION as for the launcher.

set -euo pipefail

RUN_ID="${1:?usage: $0 <RUN_ID>}"
AWS_REGION="${AWS_REGION:-us-east-1}"
S3_BUCKET="${S3_BUCKET:?S3_BUCKET unset}"

echo "[status] run_id=${RUN_ID}"

# Running / pending instances for this run.
mapfile -t INSTANCES < <(aws ec2 describe-instances --region "${AWS_REGION}" \
    --filters "Name=tag:RunId,Values=${RUN_ID}" \
              "Name=instance-state-name,Values=running,pending" \
    --query 'Reservations[*].Instances[*].[InstanceId,Tags[?Key==`Shard`].Value|[0]]' \
    --output text 2>/dev/null || true)
RUNNING="${#INSTANCES[@]}"
[[ "${RUNNING}" -eq 1 && -z "${INSTANCES[0]}" ]] && RUNNING=0

DONE_COUNT="$(aws s3 ls "s3://${S3_BUCKET}/${RUN_ID}/" --recursive --region "${AWS_REGION}" 2>/dev/null | grep -c '/_DONE$' || true)"
FAILED_COUNT="$(aws s3 ls "s3://${S3_BUCKET}/${RUN_ID}/" --recursive --region "${AWS_REGION}" 2>/dev/null | grep -c '/_FAILED_PREFLIGHT$' || true)"

echo "[status] running=${RUNNING} done=${DONE_COUNT} failed=${FAILED_COUNT}"
if [[ "${RUNNING}" -gt 0 ]]; then
    printf '  %s\n' "${INSTANCES[@]}"
fi
if [[ "${FAILED_COUNT}" -gt 0 ]]; then
    echo "[status] failed shards:"
    aws s3 ls "s3://${S3_BUCKET}/${RUN_ID}/" --recursive --region "${AWS_REGION}" \
        | grep '/_FAILED_PREFLIGHT$' | sed 's#.*/'"${RUN_ID}"'/\([^/]*\)/.*#  \1#'
fi

if [[ "${RUNNING}" -eq 0 ]]; then
    mkdir -p "/tmp/${RUN_ID}"
    aws s3 sync "s3://${S3_BUCKET}/${RUN_ID}/" "/tmp/${RUN_ID}/" --region "${AWS_REGION}" >/dev/null
    echo "[status] no instances running; synced to /tmp/${RUN_ID}/"
    echo "[status] aggregate: python3 pfloat-libm/scripts/pf-lm3-aggregate.py /tmp/${RUN_ID}/"
    if [[ "${FAILED_COUNT}" -gt 0 ]]; then
        echo "[status] relaunch failed: pfloat-libm/scripts/pf-lm3-relaunch.sh ${RUN_ID} <fn>_<K>of<M>"
    fi
fi
