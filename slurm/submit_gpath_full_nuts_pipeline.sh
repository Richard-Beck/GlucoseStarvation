#!/bin/bash

set -euo pipefail

RUN_ID_PATH=${1:-"models/gpath/v1/assessment_run_ids.txt"}
CHAIN_SPEC=${2:-"1,2,3,4"}
ITER_WARMUP=${3:-500}
ITER_SAMPLING=${4:-1000}
ADAPT_DELTA=${5:-0.99}
MAX_TREED=${6:-12}
NUM_THREADS=${7:-16}
STAN_DATA_PATH=${8:-"data/inputs/stan/gstarvation_v1/stan_ready_data.Rds"}
OUTPUT_ROOT=${9:-"data/gpath_transfer_cv_nuts"}
INIT_MODE=${10:-"optim"}
OPTIM_DATASET_LABEL=${11:-"gstarvation_v1"}
OPTIM_RUN_LABEL=${12:?Usage: bash slurm/submit_gpath_full_nuts_pipeline.sh [assessment_run_ids] [chain_spec] [warmup] [sampling] [adapt_delta] [max_treedepth] [threads] [stan_data_path] [output_root] [init_mode] [optim_dataset_label] <optim_run_label> [qos] [line_id] [direction] [clear_logs] [optim_root_override]}
QOS=${13:-"normal"}
LINE_ID=${14:-1}
DIRECTION=${15:-"low_to_high"}
CLEAR_LOGS=${16:-1}
OPTIM_ROOT_OVERRIDE=${17:-""}

CONTAINER_URI=${CONTAINER_URI:-"docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"}
BINDS=${BINDS:-"-B /home/$USER,/share,/etc/passwd,/etc/group"}

STAMP=$(date +"%Y%m%d_%H%M%S")
RUN_DIR="slurm/runs/gpath_full_nuts/${STAMP}"
MANIFEST_PATH="${RUN_DIR}/task_manifest.tsv"
JOB_INFO_PATH="${RUN_DIR}/job_ids.txt"

ARRAY_SCRIPT="slurm/jobs/gpath_full_nuts_array.sh"
MANIFEST_SCRIPT="scripts/make_gpath_full_nuts_manifest.R"

mkdir -p "$RUN_DIR" "slurm/logs"
if [ "$CLEAR_LOGS" -eq 1 ]; then
  find slurm/logs -maxdepth 1 -type f -delete
fi
export APPTAINER_NOHTTPS=1

echo "PWD: $(pwd)"
echo "Preparing full-data NUTS manifest in ${RUN_DIR}"
apptainer exec $BINDS $CONTAINER_URI \
  Rscript "$MANIFEST_SCRIPT" \
  "$MANIFEST_PATH" \
  "$RUN_ID_PATH" \
  "$CHAIN_SPEC" \
  "$LINE_ID" \
  "$DIRECTION" \
  "$STAN_DATA_PATH" \
  "$INIT_MODE" \
  "$OPTIM_DATASET_LABEL" \
  "$OPTIM_RUN_LABEL" \
  "$OPTIM_ROOT_OVERRIDE"

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
  "scripts/run_gpath_transfer_nuts.R" \
  "$ITER_WARMUP" \
  "$ITER_SAMPLING" \
  "$ADAPT_DELTA" \
  "$MAX_TREED" \
  "$NUM_THREADS" \
  "$OUTPUT_ROOT" \
  "all")

cat > "$JOB_INFO_PATH" <<EOF
run_dir=${RUN_DIR}
manifest_path=${MANIFEST_PATH}
array_job_id=${ARRAY_JOB_ID}
run_id_path=${RUN_ID_PATH}
chain_spec=${CHAIN_SPEC}
stan_data_path=${STAN_DATA_PATH}
output_root=${OUTPUT_ROOT}
init_mode=${INIT_MODE}
optim_dataset_label=${OPTIM_DATASET_LABEL}
optim_run_label=${OPTIM_RUN_LABEL}
optim_root_override=${OPTIM_ROOT_OVERRIDE}
line_id=${LINE_ID}
direction=${DIRECTION}
qos=${QOS}
EOF

echo "Array job submitted: ${ARRAY_JOB_ID}"
echo "Run metadata written to ${JOB_INFO_PATH}"
echo "run_dir=${RUN_DIR}"
