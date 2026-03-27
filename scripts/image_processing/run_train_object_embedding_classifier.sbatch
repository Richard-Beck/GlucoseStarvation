#!/bin/bash

#----------------------------------------------------------------#
# SBATCH - SLURM Resource Directives
#----------------------------------------------------------------#
#SBATCH --job-name=obj_clf_train
#SBATCH --output=logs/obj-clf-train-%j.out
#SBATCH --error=logs/obj-clf-train-%j.err
#SBATCH --time=02:00:00
#SBATCH --qos=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

echo "=========================================================="
echo "Job Started at $(date)"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Running on node: ${SLURMD_NODENAME:-NA}"
echo "=========================================================="

set +u
source /app/eb/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh
conda activate cpose
set -u

RUN_DIR=${1:-data/image_processing_runs/run_20260324_233122}
INPUT_CSV=${2:-${RUN_DIR}/object_train_features_90.csv}
OUTPUT_DIR=${3:-${RUN_DIR}/classifier_first_pass}
TEST_SIZE=${4:-0.2}
SEED=${5:-1}
MAX_ITER=${6:-2000}

mkdir -p logs "${OUTPUT_DIR}"

echo "Input CSV: ${INPUT_CSV}"
echo "Output dir: ${OUTPUT_DIR}"
echo "Test size: ${TEST_SIZE}"
echo "Seed: ${SEED}"
echo "Max iter: ${MAX_ITER}"

/home/4473331/.conda/envs/cpose/bin/python scripts/image_processing/train_object_embedding_classifier.py \
  --input_csv "${INPUT_CSV}" \
  --output_dir "${OUTPUT_DIR}" \
  --test_size "${TEST_SIZE}" \
  --seed "${SEED}" \
  --max_iter "${MAX_ITER}"

echo "Rendering R summary..."
scripts/agentRrunner.sh scripts/image_processing/plot_object_embedding_classifier_results.R "${OUTPUT_DIR}"

echo "=========================================================="
echo "Job Finished at $(date)"
echo "=========================================================="
