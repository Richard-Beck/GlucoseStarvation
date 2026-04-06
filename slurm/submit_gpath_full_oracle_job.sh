#!/bin/bash

set -euo pipefail

MANIFEST_PATH=${1:?Usage: bash slurm/submit_gpath_full_oracle_job.sh <config_manifest.tsv> [qos] [clear_logs]}
QOS=${2:-"normal"}
CLEAR_LOGS=${3:-1}

ARRAY_SCRIPT="slurm/jobs/gpath_full_oracle_array.sh"
RUN_DIR=$(dirname "$MANIFEST_PATH")
JOB_INFO_PATH="${RUN_DIR}/job_ids.txt"

mkdir -p "$RUN_DIR" "slurm/logs"
if [ "$CLEAR_LOGS" -eq 1 ]; then
  find slurm/logs -maxdepth 1 -type f -delete
fi

N_TASKS=$(awk 'END { print NR - 1 }' "$MANIFEST_PATH")
if [ "$N_TASKS" -le 0 ]; then
  echo "ERROR: manifest contains no tasks"
  exit 1
fi

ARRAY_JOB_ID=$(sbatch --parsable \
  --array="1-${N_TASKS}" \
  --qos="${QOS}" \
  "$ARRAY_SCRIPT" \
  "$MANIFEST_PATH" \
  "scripts/run_job.R")

cat > "$JOB_INFO_PATH" <<EOF
run_dir=${RUN_DIR}
manifest_path=${MANIFEST_PATH}
array_job_id=${ARRAY_JOB_ID}
qos=${QOS}
EOF

echo "Array job submitted: ${ARRAY_JOB_ID}"
echo "Run metadata written to ${JOB_INFO_PATH}"
echo "run_dir=${RUN_DIR}"
