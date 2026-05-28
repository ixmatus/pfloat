#!/usr/bin/env bash
#
# pf-hcz4 — pf-tqzz full per-release cross-check sweep launcher.
# One EC2 instance per (FnId, order) shard; 63 shards total. ADR-0049.
#
# Models the ferrodec d32-sweep pattern: per-shard cloud-init bakes
# FN_ID + ORDER + RUN_ID + S3_BUCKET; instance self-terminates on
# script end via --instance-initiated-shutdown-behavior terminate;
# results land in s3://${S3_BUCKET}/${RUN_ID}/<shard>/{result.json,
# sweep.log,_DONE}.
#
# Usage:
#   pf-hcz4-launch.sh [--only <shard>] [--smoke] [--run-id <id>]
#                     [--instance-type <type>] [--market on-demand|spot]
#                     [--dry-run]
#
# Required environment:
#   S3_BUCKET            S3 bucket name (you create + own).
#   AWS_REGION           Default us-east-1.
#   IAM_INSTANCE_PROFILE Name of an instance profile with
#                        s3:PutObject on s3://${S3_BUCKET}/* and
#                        AmazonSSMManagedInstanceCore.
#   GITHUB_REPO          The git URL the cloud-init clones from.
#                        Default: derived from `git remote get-url`.
#
# Quota check (per ferrodec lesson): the launcher inspects spot and
# on-demand vCPU quotas before fanning out. Bump via
#   aws service-quotas request-service-quota-increase
#       --service-code ec2 --quota-code L-34B43A08 \
#       --desired-value <N>
# (synchronous; may require multi-hour wait for AWS approval).

set -euo pipefail

AWS_REGION="${AWS_REGION:-us-east-1}"
INSTANCE_TYPE="${INSTANCE_TYPE:-c8g.large}"
MARKET="${MARKET:-spot}"           # "spot" or "on-demand"
RUN_ID="${RUN_ID:-pf-hcz4-$(date +%Y%m%d-%H%M%S)}"
DRY_RUN=0
ONLY_SHARD=""
SMOKE=0
S3_BUCKET="${S3_BUCKET:-}"
IAM_INSTANCE_PROFILE="${IAM_INSTANCE_PROFILE:-}"

# Default the GitHub repo URL from the local remote if not set.
GITHUB_REPO="${GITHUB_REPO:-$(git remote get-url origin 2>/dev/null || echo '')}"
GIT_SHA="${GIT_SHA:-$(git rev-parse HEAD 2>/dev/null || echo '')}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --only) ONLY_SHARD="$2"; shift 2;;
        --smoke) SMOKE=1; shift;;
        --run-id) RUN_ID="$2"; shift 2;;
        --instance-type) INSTANCE_TYPE="$2"; shift 2;;
        --market) MARKET="$2"; shift 2;;
        --dry-run) DRY_RUN=1; shift;;
        --help|-h)
            sed -n '2,/^set -e/p' "$0" | sed 's/^# \?//'
            exit 0
            ;;
        *) echo "unknown flag $1" >&2; exit 1;;
    esac
done

# ===== Prereqs =====

[[ -z "$S3_BUCKET" ]] && { echo "S3_BUCKET unset" >&2; exit 1; }
[[ -z "$IAM_INSTANCE_PROFILE" ]] && { echo "IAM_INSTANCE_PROFILE unset" >&2; exit 1; }
[[ -z "$GITHUB_REPO" ]] && { echo "GITHUB_REPO unset (couldn't derive from git remote)" >&2; exit 1; }
[[ -z "$GIT_SHA" ]] && { echo "GIT_SHA unset (couldn't derive from git rev-parse)" >&2; exit 1; }

echo "[launch] run_id=${RUN_ID}"
echo "[launch] bucket=s3://${S3_BUCKET}"
echo "[launch] repo=${GITHUB_REPO} sha=${GIT_SHA}"
echo "[launch] instance=${INSTANCE_TYPE} market=${MARKET} region=${AWS_REGION}"

# Idempotent bucket check.
if ! aws s3api head-bucket --bucket "${S3_BUCKET}" --region "${AWS_REGION}" 2>/dev/null; then
    echo "[launch] creating bucket s3://${S3_BUCKET}"
    [[ "${DRY_RUN}" -eq 1 ]] || aws s3 mb "s3://${S3_BUCKET}" --region "${AWS_REGION}"
fi

# Default VPC presence.
if ! aws ec2 describe-vpcs --region "${AWS_REGION}" \
    --filters Name=is-default,Values=true \
    --query 'Vpcs[0].VpcId' --output text 2>/dev/null | grep -q '^vpc-'; then
    echo "[launch] creating default VPC in ${AWS_REGION}"
    [[ "${DRY_RUN}" -eq 1 ]] || aws ec2 create-default-vpc --region "${AWS_REGION}"
fi

# AMI lookup: latest Ubuntu 24.04 ARM64.
AMI_ID="$(aws ssm get-parameters \
    --region "${AWS_REGION}" \
    --names /aws/service/canonical/ubuntu/server/24.04/stable/current/arm64/hvm/ebs-gp3/ami-id \
    --query 'Parameters[0].Value' --output text)"
echo "[launch] ami=${AMI_ID}"

# Quota check. Each c8g.large is 2 vCPU; 63 shards × 2 = 126 vCPU.
if [[ "${MARKET}" == "spot" ]]; then
    QUOTA="$(aws service-quotas get-service-quota --service-code ec2 \
        --quota-code L-34B43A08 --region "${AWS_REGION}" \
        --query 'Quota.Value' --output text)"
    echo "[launch] spot quota: ${QUOTA} vCPU"
    if (( $(printf '%.0f' "${QUOTA}") < 126 )); then
        echo "[launch] WARNING: spot quota ${QUOTA} < 126 vCPU needed; expect throttling" >&2
    fi
fi

# ===== Shard enumeration =====

SHARDS_FILE="$(mktemp)"
if [[ -n "${ONLY_SHARD}" ]]; then
    echo "${ONLY_SHARD}" > "${SHARDS_FILE}"
elif [[ "${SMOKE}" -eq 1 ]]; then
    echo "exp" > "${SHARDS_FILE}"
else
    # 63 shards from tests/oracle/status/*.toml. Map filename to
    # --fn-id arg: replace trailing _N on J|Y|I|K with :N.
    for f in tests/oracle/status/*.toml; do
        name="$(basename "$f" .toml)"
        case "${name}" in
            Jn_*|Yn_*|In_*|Kn_*)
                family="${name%%_*}"
                order="${name#*_}"
                echo "${family}:${order}"
                ;;
            *)
                echo "${name}"
                ;;
        esac
    done > "${SHARDS_FILE}"
fi
N_SHARDS="$(wc -l < "${SHARDS_FILE}")"
echo "[launch] ${N_SHARDS} shards to launch"

# ===== Per-shard launch =====

launch_one() {
    local fn_arg="$1"
    # Filename-safe slug for S3 prefix; replace ':' with '_'.
    local slug="${fn_arg//:/_}"
    local user_data="$(mktemp)"
    cat > "${user_data}" <<USERDATA
#!/bin/bash
# Cloud-init for pf-hcz4 shard ${slug} of run ${RUN_ID}.
set -euo pipefail
# tee to keep aws ec2 get-console-output useful (ferrodec gotcha).
exec > >(tee -a /var/log/sweep.log) 2>&1

export DEBIAN_FRONTEND=noninteractive
apt-get update -q
apt-get install -y -q git build-essential pkg-config curl \\
    python3 python3-venv python3-pip libgmp-dev libmpfr-dev

# awscli via pip (ferrodec gotcha: apt awscli flaky on noble).
pip3 install --quiet awscli

# Rust toolchain via rustup. rustup respects the repo's
# rust-toolchain.toml after the clone.
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | \\
    sh -s -- -y --default-toolchain none
. "/root/.cargo/env"

git clone "${GITHUB_REPO}" /opt/pfloat
cd /opt/pfloat
git checkout "${GIT_SHA}"

# Arb venv bootstrap (idempotent).
PFLOAT_ARB_ORACLE_VENV=/opt/pfloat-arb-oracle/venv \\
    bash scripts/setup_arb_oracle.sh

# Build release binary.
PFLOAT_ARB_ORACLE_VENV=/opt/pfloat-arb-oracle/venv \\
    cargo build --release --features differential-arb,ziv-instrumented \\
    --example pf_tqzz_sweep

# Pre-flight smoke. Bail out if the existing per-push smoke fails;
# the full shard cannot pass if the smoke does not.
if ! PFLOAT_ARB_ORACLE_VENV=/opt/pfloat-arb-oracle/venv \\
    cargo test --release --test oracle_cross_check_smoke \\
    --features differential-arb,ziv-instrumented -- --nocapture; then
    aws s3 cp /var/log/sweep.log \\
        "s3://${S3_BUCKET}/${RUN_ID}/${slug}/sweep.log"
    aws s3 cp /dev/null \\
        "s3://${S3_BUCKET}/${RUN_ID}/${slug}/_FAILED_PREFLIGHT"
    shutdown -h now
    exit 1
fi

# Run the shard. Note: --fn-id takes the dispatchable form (Yn:5 etc.)
PFLOAT_ARB_ORACLE_VENV=/opt/pfloat-arb-oracle/venv \\
    PFLOAT_GIT_SHA="${GIT_SHA}" \\
    ./target/release/examples/pf_tqzz_sweep \\
    --fn-id "${fn_arg}" \\
    --modes all \\
    --instance-type "${INSTANCE_TYPE}" \\
    --output /tmp/result.json

# Upload artifacts.
aws s3 cp /tmp/result.json \\
    "s3://${S3_BUCKET}/${RUN_ID}/${slug}/result.json"
aws s3 cp /var/log/sweep.log \\
    "s3://${S3_BUCKET}/${RUN_ID}/${slug}/sweep.log"
aws s3 cp /dev/null \\
    "s3://${S3_BUCKET}/${RUN_ID}/${slug}/_DONE"

shutdown -h now
USERDATA

    local market_args=""
    if [[ "${MARKET}" == "spot" ]]; then
        market_args="--instance-market-options MarketType=spot"
    fi

    echo "[launch] shard=${slug} fn_arg=${fn_arg}"
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        echo "[launch] (dry-run) would launch with user-data $(wc -l < "${user_data}") lines"
        rm -f "${user_data}"
        return 0
    fi

    local instance_id
    instance_id="$(aws ec2 run-instances \
        --region "${AWS_REGION}" \
        --image-id "${AMI_ID}" \
        --instance-type "${INSTANCE_TYPE}" \
        --iam-instance-profile "Name=${IAM_INSTANCE_PROFILE}" \
        --instance-initiated-shutdown-behavior terminate \
        --user-data "file://${user_data}" \
        --tag-specifications "ResourceType=instance,Tags=[{Key=RunId,Value=${RUN_ID}},{Key=Shard,Value=${slug}},{Key=Project,Value=pfloat}]" \
        ${market_args} \
        --query 'Instances[0].InstanceId' --output text)"
    echo "[launch] shard=${slug} instance=${instance_id}"

    rm -f "${user_data}"
}

while IFS= read -r shard; do
    launch_one "${shard}"
done < "${SHARDS_FILE}"

rm -f "${SHARDS_FILE}"

# Defensive shutdown fallback (ferrodec lesson B). Bounds billing if
# user-data hangs and the in-script shutdown never fires.
if [[ "${DRY_RUN}" -eq 0 ]]; then
    sleep 60   # let instances boot before sending the SSM command
    aws ssm send-command \
        --region "${AWS_REGION}" \
        --document-name AWS-RunShellScript \
        --targets "Key=tag:RunId,Values=${RUN_ID}" \
        --parameters 'commands=["shutdown -h +180 || true"]' \
        --comment "pf-hcz4 ${RUN_ID} 3h fallback" \
        > /dev/null || echo "[launch] SSM fallback scheduling skipped (no instances yet?)"
fi

echo "[launch] done. Poll with: scripts/pf-hcz4-status.sh ${RUN_ID}"
