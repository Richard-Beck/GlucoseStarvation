#!/bin/bash

set -euo pipefail

CONFIG_PATH=${1:-"scripts/selection_strategy_config.R"}
DETAILED=${2:-0}
CPUS_PER_TASK=${3:-4}
MEMORY=${4:-16G}
TIME_LIMIT=${5:-08:00:00}
QOS=${6:-normal}

if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "Config file not found: $CONFIG_PATH" >&2
  exit 1
fi

mkdir -p slurm/logs
find slurm/logs -maxdepth 1 -type f -delete

BATCH_LABEL=$(date +%Y%m%d_%H%M%S)
RUN_DIR="slurm/runs/gpath_selection_strategy/${BATCH_LABEL}"
mkdir -p "$RUN_DIR"

CONFIG_COPY="${RUN_DIR}/selection_strategy_config.R"
cp "$CONFIG_PATH" "$CONFIG_COPY"

MANIFEST_PATH="${RUN_DIR}/manifest.tsv"
Rscript scripts/make_gpath_selection_strategy_manifest.R "$CONFIG_COPY" "$MANIFEST_PATH"

N_TASKS=$(awk 'END{print NR-1}' "$MANIFEST_PATH")
if [[ -z "$N_TASKS" || "$N_TASKS" -le 0 ]]; then
  echo "No tasks found in manifest: $MANIFEST_PATH" >&2
  exit 1
fi

JOB_SCRIPT="slurm/jobs/gpath_selection_strategy_array.sh"
JOB_ID=$(sbatch \
  --array=1-"${N_TASKS}" \
  --cpus-per-task="${CPUS_PER_TASK}" \
  --mem="${MEMORY}" \
  --time="${TIME_LIMIT}" \
  --qos="${QOS}" \
  "$JOB_SCRIPT" \
  "$CONFIG_COPY" \
  "$MANIFEST_PATH" \
  "$DETAILED" | awk '{print $4}')

JOB_TABLE="${RUN_DIR}/job_ids.txt"
{
  echo "array_job_id=${JOB_ID}"
  echo "run_dir=${RUN_DIR}"
  echo "config=${CONFIG_COPY}"
  echo "manifest=${MANIFEST_PATH}"
  echo "n_tasks=${N_TASKS}"
} > "$JOB_TABLE"

echo "Selection-strategy array job submitted: ${JOB_ID}"
echo "Run metadata written to ${JOB_TABLE}"
echo "run_dir=${RUN_DIR}"
