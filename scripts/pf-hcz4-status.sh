#!/usr/bin/env bash
#
# pf-hcz4 status poller. ADR-0049.
#
# Lists running instances tagged with the given RunId, fetches one
# `tail -n 1 /var/log/sweep.log` line each via SSM, and checks S3 for
# _DONE / _FAILED_PREFLIGHT sentinels. When every shard has a
# sentinel, auto-`s3 sync` results to /tmp/<RUN_ID>/ and print the
# next-step command.
#
# Usage:
#   pf-hcz4-status.sh <RUN_ID>
#
# Environment:
#   S3_BUCKET        (required)
#   AWS_REGION       default us-east-1

set -euo pipefail

RUN_ID="${1:?usage: $0 <RUN_ID>}"
AWS_REGION="${AWS_REGION:-us-east-1}"
S3_BUCKET="${S3_BUCKET:?S3_BUCKET unset}"

# Running instances by tag.
mapfile -t INSTANCES < <(
    aws ec2 describe-instances --region "${AWS_REGION}" \
        --filters "Name=tag:RunId,Values=${RUN_ID}" \
                  "Name=instance-state-name,Values=running,pending" \
        --query 'Reservations[*].Instances[*].[InstanceId,Tags[?Key==`Shard`].Value|[0]]' \
        --output text
)

echo "[status] running/pending instances: ${#INSTANCES[@]}"

# Sentinel counts from S3.
DONE_COUNT="$(aws s3 ls "s3://${S3_BUCKET}/${RUN_ID}/" --recursive \
    --region "${AWS_REGION}" 2>/dev/null | grep -c '/_DONE$' || true)"
FAILED_COUNT="$(aws s3 ls "s3://${S3_BUCKET}/${RUN_ID}/" --recursive \
    --region "${AWS_REGION}" 2>/dev/null | grep -c '/_FAILED_PREFLIGHT$' || true)"

echo "[status] _DONE sentinels: ${DONE_COUNT}"
echo "[status] _FAILED_PREFLIGHT sentinels: ${FAILED_COUNT}"

# Per-running-instance log tail via SSM (bulk).
if [[ ${#INSTANCES[@]} -gt 0 ]]; then
    # Build the instance-id list (skip the shard-name column).
    INSTANCE_IDS=()
    for line in "${INSTANCES[@]}"; do
        INSTANCE_IDS+=("$(echo "${line}" | awk '{print $1}')")
    done

    CMD_ID="$(aws ssm send-command \
        --region "${AWS_REGION}" \
        --document-name AWS-RunShellScript \
        --instance-ids "${INSTANCE_IDS[@]}" \
        --parameters 'commands=["tail -n 1 /var/log/sweep.log 2>/dev/null || echo no-log"]' \
        --query 'Command.CommandId' --output text 2>/dev/null || true)"

    if [[ -n "${CMD_ID}" ]]; then
        sleep 3
        echo "[status] live tails:"
        for iid in "${INSTANCE_IDS[@]}"; do
            tail_line="$(aws ssm get-command-invocation \
                --region "${AWS_REGION}" \
                --command-id "${CMD_ID}" \
                --instance-id "${iid}" \
                --query 'StandardOutputContent' --output text 2>/dev/null || true)"
            shard="$(for line in "${INSTANCES[@]}"; do
                if [[ "${line}" == ${iid}* ]]; then echo "${line}" | awk '{print $2}'; fi
            done)"
            printf "  %-12s %-20s %s\n" "${iid}" "${shard}" "${tail_line:0:100}"
        done
    fi
fi

TOTAL_EXPECTED=63
SENTINELS=$((DONE_COUNT + FAILED_COUNT))
echo "[status] ${SENTINELS}/${TOTAL_EXPECTED} shards complete (DONE + FAILED)"

if [[ "${SENTINELS}" -ge "${TOTAL_EXPECTED}" && "${#INSTANCES[@]}" -eq 0 ]]; then
    echo "[status] all shards complete. Syncing s3://${S3_BUCKET}/${RUN_ID}/ → /tmp/${RUN_ID}/"
    mkdir -p "/tmp/${RUN_ID}"
    aws s3 sync "s3://${S3_BUCKET}/${RUN_ID}/" "/tmp/${RUN_ID}/" --region "${AWS_REGION}"
    echo "[status] next: python3 scripts/pf-hcz4-aggregate.py /tmp/${RUN_ID}/"
fi
