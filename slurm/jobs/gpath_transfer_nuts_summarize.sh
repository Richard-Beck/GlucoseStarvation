#!/bin/bash
#SBATCH --job-name=gpath_xnsum
#SBATCH --output=slurm/logs/%x_%j.out
#SBATCH --error=slurm/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --qos=normal

set -euo pipefail

MODEL_NAME=${1:-"gpath"}
RUN_ID=${2:?run id required}
OUTPUT_ROOT=${3:-"data/gpath_transfer_cv_nuts"}
SUMMARY_SCRIPT=${4:-"scripts/summarize_gpath_transfer_nuts.R"}

CONTAINER_URI=${CONTAINER_URI:-"docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"}
BINDS=${BINDS:-"-B /home/$USER,/share,/etc/passwd,/etc/group"}

export TZ="US/Eastern"
export APPTAINER_NOHTTPS=1

mkdir -p slurm/logs

echo "PWD: $(pwd)"
echo "Summarizing transfer CV NUTS outputs for ${MODEL_NAME} ${RUN_ID}"

apptainer exec $BINDS $CONTAINER_URI \
  Rscript "$SUMMARY_SCRIPT" \
  "$MODEL_NAME" \
  "$RUN_ID" \
  "$OUTPUT_ROOT"
