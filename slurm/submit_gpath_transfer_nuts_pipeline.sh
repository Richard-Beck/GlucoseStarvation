#!/bin/bash

set -euo pipefail

MODEL_NAME=${1:-"gpath"}
R_VAL=${2:-1}
P_VAL=${3:-1}
W_VAL=${4:-0}
CONSTRAINT_FLAG=${5:-0}
WASTE_MECH_FLAG=${6:-1}

LINE_SPEC=${7:-"all"}
DIRECTION_SPEC=${8:-"low_to_high,high_to_low"}
FIT_SPEC=${9:-"null,transfer,oracle"}
CHAIN_SPEC=${10:-"1,2,3,4"}

ITER_WARMUP=${11:-500}
ITER_SAMPLING=${12:-1000}
ADAPT_DELTA=${13:-0.99}
MAX_TREED=${14:-12}
NUM_THREADS=${15:-16}
STAN_DATA_PATH=${16:-"data/inputs/stan/gstarvation_v1/stan_ready_data.Rds"}
OUTPUT_ROOT=${17:-"data/gpath_transfer_cv_nuts"}

CONTAINER_URI=${CONTAINER_URI:-"docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"}
BINDS=${BINDS:-"-B /home/$USER,/share,/etc/passwd,/etc/group"}

RUN_ID="${R_VAL}R_${P_VAL}P_${W_VAL}W_C${CONSTRAINT_FLAG}_M${WASTE_MECH_FLAG}"
STAMP=$(date +"%Y%m%d_%H%M%S")
RUN_DIR="slurm/runs/gpath_transfer_nuts/${RUN_ID}/${STAMP}"
MANIFEST_PATH="${RUN_DIR}/task_manifest.tsv"
JOB_INFO_PATH="${RUN_DIR}/job_ids.txt"

ARRAY_SCRIPT="slurm/jobs/gpath_transfer_nuts_array.sh"
SUMMARY_SCRIPT="slurm/jobs/gpath_transfer_nuts_summarize.sh"
MANIFEST_SCRIPT="scripts/make_gpath_transfer_nuts_manifest.R"

mkdir -p "$RUN_DIR" "slurm/logs"
find slurm/logs -maxdepth 1 -type f -delete
export APPTAINER_NOHTTPS=1

echo "PWD: $(pwd)"
echo "Preparing transfer CV NUTS manifest in ${RUN_DIR}"
apptainer exec $BINDS $CONTAINER_URI \
  Rscript "$MANIFEST_SCRIPT" \
  "$MANIFEST_PATH" \
  "$STAN_DATA_PATH" \
  "$LINE_SPEC" \
  "$DIRECTION_SPEC" \
  "$FIT_SPEC" \
  "$CHAIN_SPEC"

N_TASKS=$(awk 'END { print NR - 1 }' "$MANIFEST_PATH")
if [ "$N_TASKS" -le 0 ]; then
  echo "ERROR: manifest contains no tasks"
  exit 1
fi

ARRAY_JOB_ID=$(sbatch --parsable \
  --array="1-${N_TASKS}" \
  "$ARRAY_SCRIPT" \
  "$MODEL_NAME" \
  "$R_VAL" \
  "$P_VAL" \
  "$W_VAL" \
  "$CONSTRAINT_FLAG" \
  "$WASTE_MECH_FLAG" \
  "$MANIFEST_PATH" \
  "scripts/run_gpath_transfer_nuts.R" \
  "$ITER_WARMUP" \
  "$ITER_SAMPLING" \
  "$ADAPT_DELTA" \
  "$MAX_TREED" \
  "$NUM_THREADS" \
  "$OUTPUT_ROOT")

SUMMARY_JOB_ID=$(sbatch --parsable \
  --dependency="afterany:${ARRAY_JOB_ID}" \
  "$SUMMARY_SCRIPT" \
  "$MODEL_NAME" \
  "$RUN_ID" \
  "$OUTPUT_ROOT" \
  "scripts/summarize_gpath_transfer_nuts.R")

cat > "$JOB_INFO_PATH" <<EOF
run_id=${RUN_ID}
run_dir=${RUN_DIR}
manifest_path=${MANIFEST_PATH}
array_job_id=${ARRAY_JOB_ID}
summary_job_id=${SUMMARY_JOB_ID}
output_root=${OUTPUT_ROOT}
EOF

echo "Array job submitted: ${ARRAY_JOB_ID}"
echo "Summary job submitted: ${SUMMARY_JOB_ID}"
echo "Run metadata written to ${JOB_INFO_PATH}"
