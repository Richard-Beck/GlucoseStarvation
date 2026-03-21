#!/bin/bash
#SBATCH --job-name=gpath_select
#SBATCH --output=slurm/logs/%x_%A_%a.out
#SBATCH --error=slurm/logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --qos=normal

set -euo pipefail

CONFIG_PATH=${1:?config path required}
MANIFEST_PATH=${2:?manifest path required}
DETAILED=${3:-0}
TASK_SCRIPT=${4:-"scripts/run_gpath_selection_strategy_task.R"}

TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}

CONTAINER_URI=${CONTAINER_URI:-"docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"}
BINDS=${BINDS:-"-B /home/$USER,/share,/etc/passwd,/etc/group"}

export OMP_NUM_THREADS=1
export TZ="US/Eastern"
export APPTAINER_NOHTTPS=1

mkdir -p slurm/logs

echo "PWD: $(pwd)"
echo "Starting selection-strategy task on $(hostname)"
echo "Slurm job: ${SLURM_JOB_ID} | array task: ${TASK_ID}"
echo "Config: ${CONFIG_PATH}"
echo "Manifest: ${MANIFEST_PATH}"
echo "Detailed trajectories: ${DETAILED}"

id "$USER" > /dev/null 2>&1
sleep 2

apptainer exec $BINDS $CONTAINER_URI \
  Rscript "$TASK_SCRIPT" \
  "config=${CONFIG_PATH}" \
  "manifest=${MANIFEST_PATH}" \
  "task_id=${TASK_ID}" \
  "detailed=${DETAILED}"
