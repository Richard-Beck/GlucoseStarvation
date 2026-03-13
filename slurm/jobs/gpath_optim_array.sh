#!/bin/bash
#SBATCH --job-name=gpath_opt
#SBATCH --output=slurm/logs/%x_%A_%a.out
#SBATCH --error=slurm/logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --qos=normal

set -euo pipefail

MODEL_NAME=${1:-"gpath"}
MODEL_VERSION=${2:-"v1"}
RUN_ID=${3:?run_id required}
NUM_THREADS=${4:-8}
STAN_DATA_PATH=${5:-"data/inputs/stan/gstarvation_v1/stan_ready_data.Rds"}
DATASET_LABEL=${6:-"gstarvation_v1"}
OUTPUT_DIR=${7:?output_dir required}

CONTAINER_URI=${CONTAINER_URI:-"docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"}
BINDS=${BINDS:-"-B /home/$USER,/share,/etc/passwd,/etc/group"}

export OMP_NUM_THREADS=1
export TZ="US/Eastern"
export APPTAINER_NOHTTPS=1

mkdir -p slurm/logs

echo "PWD: $(pwd)"
echo "Starting optimisation task ${SLURM_ARRAY_TASK_ID} for ${RUN_ID} on $(hostname)"

id "$USER" > /dev/null 2>&1
sleep 2

apptainer exec $BINDS $CONTAINER_URI \
  Rscript pipelines/gpath/optim/run_start.R \
  "$MODEL_NAME" \
  "$MODEL_VERSION" \
  "$RUN_ID" \
  "$NUM_THREADS" \
  "$SLURM_ARRAY_TASK_ID" \
  "$STAN_DATA_PATH" \
  "$DATASET_LABEL" \
  "$OUTPUT_DIR"
