#!/bin/bash

#SBATCH --job-name=sum159_comp_reseg
#SBATCH --time=02:00:00
#SBATCH --qos=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --gres=gpu:a30:1

set -euo pipefail

MANIFEST_TSV=$1
RUN_DIR=$2
CLASSIFIER_PATH=$3
PROJECT_ROOT=$4

cd "${PROJECT_ROOT}"

set +u
source /app/eb/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh
conda activate cpose
set -u

python "${RUN_DIR}/code/batch_segment_and_classify_competition.py" \
  --manifest_tsv "${MANIFEST_TSV}" \
  --run_dir "${RUN_DIR}" \
  --series_index "${SLURM_ARRAY_TASK_ID}" \
  --classifier_path "${CLASSIFIER_PATH}" \
  --cellpose_batch_size 8
