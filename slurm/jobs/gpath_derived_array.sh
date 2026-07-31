#!/bin/bash
#SBATCH --job-name=gpath_derived
#SBATCH --output=slurm/logs/%x_%A_%a.out
#SBATCH --error=slurm/logs/%x_%A_%a.err

set -euo pipefail

SCRIPT=${1:?R builder script required}
CONFIG=${2:?release config required}
CONTAINER_URI=${CONTAINER_URI:-"docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"}
BINDS=${BINDS:-"-B /home/$USER,/share,/etc/passwd,/etc/group"}
export TZ="US/Eastern"
export APPTAINER_NOHTTPS=1
mkdir -p slurm/logs

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  apptainer exec $BINDS "$CONTAINER_URI" Rscript "$SCRIPT" "$CONFIG" "$SLURM_ARRAY_TASK_ID"
else
  apptainer exec $BINDS "$CONTAINER_URI" Rscript "$SCRIPT" "$CONFIG"
fi
