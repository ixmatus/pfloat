#!/usr/bin/env bash
#
# pf-lm3 — pfloat-libm exhaustive f32 verification sweep launcher.
# One EC2 instance per function; each instance splits its function's 2^32
# input space across all its vCPUs and runs the sub-shards in parallel.
# A natural evolution of pfloat's pf-hcz4 "one FnId per instance" (ADR-0049)
# for the exhaustive-2^32 scale: the per-core sharding is the capability
# pf-hcz4 lacked. MPFR-only, no Arb venv (ADR-0058).
#
# Per-instance cloud-init bakes FN + RUN_ID + S3_BUCKET; the instance
# self-terminates on completion; each function's NPROC sub-shard results
# land in s3://${S3_BUCKET}/${RUN_ID}/<fn>/result_<k>.json. The aggregator
# merges them and validates [0, 2^32) coverage.
#
# Usage:
#   pf-lm3-launch.sh --yes-i-will-pay [--instance-type <type>]
#                    [--directed-sample N] [--only <fn>] [--smoke]
#                    [--run-id <id>] [--market on-demand|spot]
#   pf-lm3-launch.sh --dry-run [...]      # print the fan-out, launch nothing
#
# A real launch REQUIRES --yes-i-will-pay (the EC2-spend gate). ALWAYS run
# --smoke (one function on one instance) FIRST: it validates the cloud-init
# end to end and gives the true per-function wall time to size the spend.
#
# Required environment:
#   S3_BUCKET            S3 bucket name (you create + own).
#   IAM_INSTANCE_PROFILE Instance profile with s3:PutObject on
#                        s3://${S3_BUCKET}/* and AmazonSSMManagedInstanceCore.
#   AWS_REGION           Default us-east-1.
#   GITHUB_REPO          Git URL cloud-init clones. Default: git remote.
#
# Instance sizing: cost is ~constant per core-hour across c8g sizes, so
# pick the size for wall time and your vCPU quota. The default c8g.8xlarge
# (32 vCPU) finishes the slowest function (atanh, ~120-300 core-hours) in
# ~4-9h; c8g.16xlarge (64 vCPU) halves that. 25 instances at 32 vCPU need
# an 800-vCPU on-demand quota; request a bump if needed.

set -euo pipefail

AWS_REGION="${AWS_REGION:-us-east-1}"
INSTANCE_TYPE="${INSTANCE_TYPE:-c8g.8xlarge}"
MARKET="${MARKET:-on-demand}"      # on-demand: spot reclaim wastes a long sweep
DIRECTED_SAMPLE="${DIRECTED_SAMPLE:-1048576}"
# Wall-clock billing guard (minutes), issued post-launch via SSM. The slowest
# saturating functions (atanh, tanh) legitimately run many hours; size this
# above their measured wall time. 720 (12h) was too tight for a 32-vCPU atanh
# (the in-domain shards finished right at the guard); a 64-vCPU re-run halves
# that, but keep generous headroom.
FALLBACK_MINUTES="${FALLBACK_MINUTES:-720}"
# Targeted re-run knobs. Default empty = full-space sweep, instance-nproc
# shard-count, per-function S3 dir (the normal path). Set all three to re-sweep
# a subset of shard indices at a finer shard-count into a SEPARATE S3 prefix —
# used to finish a load-imbalanced tail (a few slow contiguous shards spread
# across all cores) without recomputing the banked shards or colliding their
# result_<k>.json names. The aggregator groups by the JSON `function` field, so
# the separate prefix still merges into the same function row.
SHARD_COUNT_OVERRIDE="${SHARD_COUNT_OVERRIDE:-}"
SHARD_INDEX_LIST="${SHARD_INDEX_LIST:-}"
OUTPUT_SUBDIR_OVERRIDE="${OUTPUT_SUBDIR:-}"
RUN_ID="${RUN_ID:-pf-lm3-$(date +%Y%m%d-%H%M%S)}"
DRY_RUN=0
PAY_OK=0
ONLY_FN=""
SMOKE=0
S3_BUCKET="${S3_BUCKET:-}"
IAM_INSTANCE_PROFILE="${IAM_INSTANCE_PROFILE:-}"

GITHUB_REPO="${GITHUB_REPO:-$(git remote get-url origin 2>/dev/null || echo '')}"
GIT_SHA="${GIT_SHA:-$(git rev-parse HEAD 2>/dev/null || echo '')}"

# The 25 unary functions (must match LibmFnId::UNARY).
FUNCTIONS=(
    exp exp2 exp10 expm1
    ln log2 log10 log1p
    sqrt cbrt
    sin cos tan cot sec csc
    asin acos atan
    sinh cosh tanh
    asinh acosh atanh
)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --yes-i-will-pay) PAY_OK=1; shift;;
        --dry-run) DRY_RUN=1; shift;;
        --only) ONLY_FN="$2"; shift 2;;
        --smoke) SMOKE=1; shift;;
        --run-id) RUN_ID="$2"; shift 2;;
        --directed-sample) DIRECTED_SAMPLE="$2"; shift 2;;
        --instance-type) INSTANCE_TYPE="$2"; shift 2;;
        --market) MARKET="$2"; shift 2;;
        --help|-h) sed -n '2,/^set -e/p' "$0" | sed 's/^# \?//'; exit 0;;
        *) echo "unknown flag $1" >&2; exit 1;;
    esac
done

[[ -z "$S3_BUCKET" ]] && { echo "S3_BUCKET unset" >&2; exit 1; }
[[ -z "$IAM_INSTANCE_PROFILE" ]] && { echo "IAM_INSTANCE_PROFILE unset" >&2; exit 1; }
[[ -z "$GITHUB_REPO" ]] && { echo "GITHUB_REPO unset (couldn't derive from git remote)" >&2; exit 1; }
[[ -z "$GIT_SHA" ]] && { echo "GIT_SHA unset (couldn't derive from git rev-parse)" >&2; exit 1; }

# ===== Function (instance) enumeration =====

if [[ -n "${ONLY_FN}" ]]; then
    # Accept a comma-separated list so a targeted re-run can relaunch several
    # functions (e.g. --only atanh,tanh) under one RUN_ID and one fallback.
    IFS=',' read -r -a LAUNCH <<< "${ONLY_FN}"
elif [[ "${SMOKE}" -eq 1 ]]; then
    LAUNCH=("exp")
else
    LAUNCH=("${FUNCTIONS[@]}")
fi

echo "[launch] run_id=${RUN_ID}"
echo "[launch] bucket=s3://${S3_BUCKET} repo=${GITHUB_REPO} sha=${GIT_SHA}"
echo "[launch] instance=${INSTANCE_TYPE} market=${MARKET} region=${AWS_REGION}"
echo "[launch] ${#LAUNCH[@]} instances (one per function): ${LAUNCH[*]}"
echo "[launch] each instance self-shards across its vCPUs; directed_sample=${DIRECTED_SAMPLE}/sub-shard"

if [[ "${DRY_RUN}" -eq 0 && "${PAY_OK}" -eq 0 ]]; then
    echo "[launch] REFUSING to launch: pass --yes-i-will-pay to spend, or --dry-run to preview" >&2
    echo "[launch] strongly recommended: run --smoke (one function) first to calibrate + validate" >&2
    exit 2
fi

# ===== AWS prereqs =====

if ! aws s3api head-bucket --bucket "${S3_BUCKET}" --region "${AWS_REGION}" 2>/dev/null; then
    echo "[launch] creating bucket s3://${S3_BUCKET}"
    [[ "${DRY_RUN}" -eq 1 ]] || aws s3 mb "s3://${S3_BUCKET}" --region "${AWS_REGION}"
fi
AMI_ID="$(aws ssm get-parameters --region "${AWS_REGION}" \
    --names /aws/service/canonical/ubuntu/server/24.04/stable/current/arm64/hvm/ebs-gp3/ami-id \
    --query 'Parameters[0].Value' --output text 2>/dev/null || echo 'DRY')"
echo "[launch] ami=${AMI_ID}"

# ===== Per-function launch =====

launch_one() {
    local fn="$1"
    local user_data; user_data="$(mktemp)"
    # Unquoted heredoc: ${...} expands now (baked into user-data); \${...}
    # and \$(...) stay literal for instance-time evaluation.
    cat > "${user_data}" <<USERDATA
#!/bin/bash
# Cloud-init for pf-lm3 function ${fn} of run ${RUN_ID}.
set -euo pipefail
exec > >(tee -a /var/log/sweep.log) 2>&1
export HOME=/root
export DEBIAN_FRONTEND=noninteractive
# Targeted re-run params baked from the launcher (empty => instance-time
# defaults computed below: full nproc sharding into the per-function dir).
SHARD_COUNT_OVERRIDE="${SHARD_COUNT_OVERRIDE}"
SHARD_INDEX_LIST="${SHARD_INDEX_LIST}"
OUTPUT_SUBDIR="${OUTPUT_SUBDIR_OVERRIDE}"
apt-get update -q
# m4: gmp-mpfr-sys builds GMP/MPFR from vendored source (ADR-0049). No
# python/Arb: the libm harness is MPFR-only (ADR-0058).
apt-get install -y -q git build-essential m4 pkg-config curl libgmp-dev libmpfr-dev python3-pip
pip3 install --quiet --break-system-packages awscli

trap 'rc=\$?; [ \$rc -eq 0 ] && exit 0; aws s3 cp /var/log/sweep.log "s3://${S3_BUCKET}/${RUN_ID}/${fn}/sweep.log" 2>/dev/null || true; echo done | aws s3 cp - "s3://${S3_BUCKET}/${RUN_ID}/${fn}/_FAILED_PREFLIGHT" 2>/dev/null || true; shutdown -h now' EXIT

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain none
. "/root/.cargo/env"

git clone "${GITHUB_REPO}" /opt/pfloat
cd /opt/pfloat
git checkout "${GIT_SHA}"

# Pinned nightly with retries (lazy install flakes under concurrency).
TC=\$(grep -oE 'nightly-[0-9-]+' rust-toolchain.toml | head -1)
for i in 1 2 3 4 5; do
    rustup toolchain install "\$TC" --profile minimal && break
    echo "[pf-lm3] rustup install attempt \$i failed; retry in 20s"; sleep 20
done

cargo build -p pfloat-libm --release --features differential-mpfr --example libm_sweep

# Pre-flight: the per-push smoke must pass or the shard cannot.
if ! cargo test -p pfloat-libm --release --features differential-mpfr \\
    --test libm_smoke_gate -- --nocapture; then
    aws s3 cp /var/log/sweep.log "s3://${S3_BUCKET}/${RUN_ID}/${fn}/sweep.log"
    echo done | aws s3 cp - "s3://${S3_BUCKET}/${RUN_ID}/${fn}/_FAILED_PREFLIGHT"
    shutdown -h now; exit 1
fi

# Split this function's 2^32 space across all vCPUs; run in parallel. A
# sub-shard exit 1 means it FOUND value mismatches (a real finding, not a
# crash), so do not abort the batch; the aggregator reads them from the
# uploaded JSON.
NPROC=\$(nproc)
# Normal path: shard-count = nproc, all indices 0..nproc-1, output -> <fn>/.
# Targeted-tail path: a finer SHARD_COUNT over an explicit index subset, into a
# SEPARATE OUTSUB prefix so result_<k>.json names never collide with the banked
# full-run shards (the aggregator merges by the JSON function-field, not dir).
SHARD_COUNT="\${SHARD_COUNT_OVERRIDE:-\$NPROC}"
INDICES="\${SHARD_INDEX_LIST:-\$(seq 0 \$((SHARD_COUNT - 1)))}"
OUTSUB="\${OUTPUT_SUBDIR:-${fn}}"
echo "[pf-lm3] ${fn}: shard-count=\$SHARD_COUNT, \$(echo \$INDICES | wc -w) indices -> ${RUN_ID}/\$OUTSUB"
set +e
# Each sub-shard uploads its OWN result the instant it finishes (not in one
# batch after all NPROC complete). A wall-clock guard that fires mid-run then
# loses only the still-running shards; every completed shard is already in S3,
# so the aggregator sees partial coverage and a targeted re-run fills only the
# gap. The earlier batch-at-end design lost an entire 24/32 function when the
# guard fired inside the wait loop. run_and_upload is backgrounded per shard;
# PIDs are waited on EXPLICITLY (a bare \`wait\` would also block on the \`tee\`
# from the \`exec > >(tee ...)\` redirection, which never exits).
run_and_upload() {
    local k="\$1"
    PFLOAT_GIT_SHA="${GIT_SHA}" ./target/release/examples/libm_sweep \\
        --function "${fn}" --width f32 --exhaustive \\
        --shard-index "\$k" --shard-count "\$SHARD_COUNT" \\
        --directed-sample "${DIRECTED_SAMPLE}" \\
        --instance-type "${INSTANCE_TYPE}" \\
        --output-json "/tmp/result_\$k.json" > "/tmp/sub_\$k.log" 2>&1
    cat "/tmp/sub_\$k.log" >> /var/log/sweep.log 2>/dev/null || true
    aws s3 cp "/tmp/result_\$k.json" "s3://${S3_BUCKET}/${RUN_ID}/\$OUTSUB/result_\$k.json" || true
}
sweep_pids=""
for k in \$INDICES; do
    run_and_upload "\$k" &
    sweep_pids="\$sweep_pids \$!"
done
for p in \$sweep_pids; do wait "\$p"; done
set -e

aws s3 cp /var/log/sweep.log "s3://${S3_BUCKET}/${RUN_ID}/\$OUTSUB/sweep.log"
echo done | aws s3 cp - "s3://${S3_BUCKET}/${RUN_ID}/\$OUTSUB/_DONE"
shutdown -h now
USERDATA

    local market_args=""
    [[ "${MARKET}" == "spot" ]] && market_args="--instance-market-options MarketType=spot"

    if [[ "${DRY_RUN}" -eq 1 ]]; then
        echo "[launch] (dry-run) ${fn}: $(wc -l < "${user_data}") user-data lines"
        rm -f "${user_data}"; return 0
    fi

    local instance_id
    instance_id="$(aws ec2 run-instances --region "${AWS_REGION}" \
        --image-id "${AMI_ID}" --instance-type "${INSTANCE_TYPE}" \
        --iam-instance-profile "Name=${IAM_INSTANCE_PROFILE}" \
        --instance-initiated-shutdown-behavior terminate \
        --user-data "file://${user_data}" \
        --tag-specifications "ResourceType=instance,Tags=[{Key=RunId,Value=${RUN_ID}},{Key=Shard,Value=${fn}},{Key=Project,Value=pfloat-libm}]" \
        ${market_args} \
        --query 'Instances[0].InstanceId' --output text)"
    echo "[launch] ${fn} instance=${instance_id}"
    rm -f "${user_data}"
}

for fn in "${LAUNCH[@]}"; do
    launch_one "${fn}"
done

# Defensive shutdown fallback: bound billing if a cloud-init hangs. The
# slowest function may legitimately run many hours (set FALLBACK_MINUTES
# above its measured wall time); a hung build is caught far sooner by the
# absence of a _DONE sentinel in the status poller. Per-shard incremental
# upload (above) means a guard firing mid-run loses only the running shards,
# not the whole function.
if [[ "${DRY_RUN}" -eq 0 ]]; then
    sleep 60
    aws ssm send-command --region "${AWS_REGION}" \
        --document-name AWS-RunShellScript \
        --targets "Key=tag:RunId,Values=${RUN_ID}" \
        --parameters "commands=[\"shutdown -h +${FALLBACK_MINUTES} || true\"]" \
        --comment "pf-lm3 ${RUN_ID} ${FALLBACK_MINUTES}m fallback" \
        > /dev/null || echo "[launch] SSM fallback scheduling skipped"
fi

echo "[launch] done. Poll with: pfloat-libm/scripts/pf-lm3-status.sh ${RUN_ID}"
