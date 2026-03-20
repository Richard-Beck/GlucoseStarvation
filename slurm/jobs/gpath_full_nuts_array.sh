#!/bin/bash
#SBATCH --job-name=gpath_fnuts
#SBATCH --output=slurm/logs/%x_%A_%a.out
#SBATCH --error=slurm/logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --qos=normal

set -euo pipefail

MANIFEST_PATH=${1:?manifest path required}
TASK_SCRIPT=${2:-"scripts/run_gpath_transfer_nuts.R"}
ITER_WARMUP=${3:-500}
ITER_SAMPLING=${4:-1000}
ADAPT_DELTA=${5:-0.99}
MAX_TREED=${6:-12}
NUM_THREADS=${7:-16}
OUTPUT_ROOT=${8:-"data/gpath_transfer_cv_nuts"}
PLOIDY_EFFECT_MASK_SPEC=${9:-"all"}

TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
TASK_LINE=$(awk -v task_id="$TASK_ID" 'BEGIN { FS="\t" } NR == (task_id + 1) { print; exit }' "$MANIFEST_PATH")
if [ -z "$TASK_LINE" ]; then
  echo "ERROR: no manifest row found for task $TASK_ID in $MANIFEST_PATH"
  exit 1
fi

IFS=$'\t' read -r _TASK RUN_ID CHAIN_ID LINE_ID DIRECTION STAN_DATA_PATH INIT_MODE INIT_SOURCE <<< "$TASK_LINE"

IFS='_' read -r r_part p_part w_part c_part m_part <<< "$RUN_ID"
R_VAL=${r_part//R/}
P_VAL=${p_part//P/}
W_VAL=${w_part//W/}
CONSTRAINT_FLAG=${c_part//C/}
WASTE_MECH_FLAG=${m_part//M/}

CONTAINER_URI=${CONTAINER_URI:-"docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"}
BINDS=${BINDS:-"-B /home/$USER,/share,/etc/passwd,/etc/group"}

export OMP_NUM_THREADS=1
export TZ="US/Eastern"
export APPTAINER_NOHTTPS=1

mkdir -p slurm/logs

echo "PWD: $(pwd)"
echo "Starting full-data NUTS task on $(hostname)"
echo "Slurm job: ${SLURM_JOB_ID} | array task: ${TASK_ID}"
echo "Run: ${RUN_ID} chain=${CHAIN_ID}"
echo "Init mode=${INIT_MODE}"
echo "Init source=${INIT_SOURCE}"

id "$USER" > /dev/null 2>&1
sleep 2

apptainer exec $BINDS $CONTAINER_URI \
  Rscript "$TASK_SCRIPT" \
  "gpath" \
  "$R_VAL" \
  "$P_VAL" \
  "$W_VAL" \
  "$CONSTRAINT_FLAG" \
  "$WASTE_MECH_FLAG" \
  "$LINE_ID" \
  "$DIRECTION" \
  "oracle" \
  "$CHAIN_ID" \
  "$ITER_WARMUP" \
  "$ITER_SAMPLING" \
  "$ADAPT_DELTA" \
  "$MAX_TREED" \
  "$NUM_THREADS" \
  "$STAN_DATA_PATH" \
  "$OUTPUT_ROOT" \
  "$PLOIDY_EFFECT_MASK_SPEC" \
  "$INIT_MODE" \
  "$INIT_SOURCE"
