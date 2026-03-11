#!/bin/bash
#SBATCH --job-name=gpath_xsum
#SBATCH --output=slurm/logs/%x_%j.out
#SBATCH --error=slurm/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --qos=normal

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "$REPO_ROOT"

MODEL_NAME=${1:-"gpath"}
RUN_ID=${2:?run id required}
OUTPUT_ROOT=${3:-"data/gpath_transfer_cv"}
SUMMARY_SCRIPT=${4:-"scripts/summarize_gpath_transfer_cv.R"}

CONTAINER_URI=${CONTAINER_URI:-"docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"}
BINDS=${BINDS:-"-B /home/$USER,/share,/etc/passwd,/etc/group"}

export TZ="US/Eastern"
export APPTAINER_NOHTTPS=1

mkdir -p slurm/logs

echo "Summarizing transfer CV outputs for ${MODEL_NAME} ${RUN_ID}"

apptainer exec $BINDS $CONTAINER_URI \
  Rscript "$SUMMARY_SCRIPT" \
  "$MODEL_NAME" \
  "$RUN_ID" \
  "$OUTPUT_ROOT"
