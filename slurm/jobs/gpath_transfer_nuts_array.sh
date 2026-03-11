#!/bin/bash
#SBATCH --job-name=gpath_xnuts
#SBATCH --output=slurm/logs/%x_%A_%a.out
#SBATCH --error=slurm/logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --qos=normal

set -euo pipefail

MODEL_NAME=${1:-"gpath"}
R_VAL=${2:-1}
P_VAL=${3:-1}
W_VAL=${4:-0}
CONSTRAINT_FLAG=${5:-0}
WASTE_MECH_FLAG=${6:-1}
MANIFEST_PATH=${7:?manifest path required}
TASK_SCRIPT=${8:-"scripts/run_gpath_transfer_nuts.R"}
ITER_WARMUP=${9:-500}
ITER_SAMPLING=${10:-1000}
ADAPT_DELTA=${11:-0.99}
MAX_TREED=${12:-12}
NUM_THREADS=${13:-16}
OUTPUT_ROOT=${14:-"data/gpath_transfer_cv_nuts"}

TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
TASK_LINE=$(awk -v task_id="$TASK_ID" 'BEGIN { FS="\t" } NR == (task_id + 1) { print; exit }' "$MANIFEST_PATH")
if [ -z "$TASK_LINE" ]; then
  echo "ERROR: no manifest row found for task $TASK_ID in $MANIFEST_PATH"
  exit 1
fi

IFS=$'\t' read -r _TASK LINE_ID DIRECTION FIT_TYPE CHAIN_ID OBSERVED_VALUE HOLDOUT_VALUE HOLDOUT_WELLS HOLDOUT_OBS STAN_DATA_PATH <<< "$TASK_LINE"

CONTAINER_URI=${CONTAINER_URI:-"docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"}
BINDS=${BINDS:-"-B /home/$USER,/share,/etc/passwd,/etc/group"}

export OMP_NUM_THREADS=1
export TZ="US/Eastern"
export APPTAINER_NOHTTPS=1

mkdir -p slurm/logs

echo "PWD: $(pwd)"
echo "Starting transfer CV NUTS task on $(hostname)"
echo "Slurm job: ${SLURM_JOB_ID} | array task: ${TASK_ID}"
echo "Task: line=${LINE_ID} direction=${DIRECTION} fit=${FIT_TYPE} chain=${CHAIN_ID}"

id "$USER" > /dev/null 2>&1
sleep 2

apptainer exec $BINDS $CONTAINER_URI \
  Rscript "$TASK_SCRIPT" \
  "$MODEL_NAME" \
  "$R_VAL" \
  "$P_VAL" \
  "$W_VAL" \
  "$CONSTRAINT_FLAG" \
  "$WASTE_MECH_FLAG" \
  "$LINE_ID" \
  "$DIRECTION" \
  "$FIT_TYPE" \
  "$CHAIN_ID" \
  "$ITER_WARMUP" \
  "$ITER_SAMPLING" \
  "$ADAPT_DELTA" \
  "$MAX_TREED" \
  "$NUM_THREADS" \
  "$STAN_DATA_PATH" \
  "$OUTPUT_ROOT"
