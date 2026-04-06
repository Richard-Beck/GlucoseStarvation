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
TASK_SCRIPT=${2:-"scripts/run_job.R"}

TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
TASK_LINE=$(awk -v task_id="$TASK_ID" 'BEGIN { FS="\t" } NR == (task_id + 1) { print; exit }' "$MANIFEST_PATH")
if [ -z "$TASK_LINE" ]; then
  echo "ERROR: no manifest row found for task $TASK_ID in $MANIFEST_PATH"
  exit 1
fi

IFS=$'\t' read -r _TASK RUN_ID CHAIN_ID CONFIG_PATH OUTPUT_DIR RUN_TAG <<< "$TASK_LINE"

export OMP_NUM_THREADS=1
export TZ="US/Eastern"

mkdir -p slurm/logs

echo "PWD: $(pwd)"
echo "Starting full-oracle NUTS task on $(hostname)"
echo "Slurm job: ${SLURM_JOB_ID} | array task: ${TASK_ID}"
echo "Run: ${RUN_ID} chain=${CHAIN_ID} tag=${RUN_TAG}"
echo "Config: ${CONFIG_PATH}"
echo "Output dir: ${OUTPUT_DIR}"

id "$USER" > /dev/null 2>&1
sleep 2

scripts/agentRrunner.sh "$TASK_SCRIPT" "$CONFIG_PATH"
