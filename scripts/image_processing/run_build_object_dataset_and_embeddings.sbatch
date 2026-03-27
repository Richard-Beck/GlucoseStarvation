#!/bin/bash

#----------------------------------------------------------------#
# SBATCH - SLURM Resource Directives
#----------------------------------------------------------------#
#SBATCH --job-name=obj_embed_90
#SBATCH --output=logs/obj-embed-90-%j.out
#SBATCH --error=logs/obj-embed-90-%j.err
#SBATCH --time=01:00:00
#SBATCH --qos=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --gres=gpu:a30:1

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
CELLPOSE_BATCH_SIZE=${2:-32}
EMBED_BATCH_SIZE=${3:-512}

STACK_TIF="${RUN_DIR}/annotation_stack_90.tif"
ANNOTATION_CSV="${RUN_DIR}/annotation_stack_90_Results.csv"
MANIFEST_CSV="${RUN_DIR}/annotation_stack_90_manifest.csv"
OUT_LABELS_TIF="${RUN_DIR}/objects_labels_90.tif"
OUT_FEATURES="${RUN_DIR}/object_train_features_90.csv"
OUT_META="${RUN_DIR}/object_train_features_90_meta.json"

mkdir -p logs "${RUN_DIR}"

echo "Conda env: cpose"
echo "Run dir: ${RUN_DIR}"
echo "Cellpose batch size: ${CELLPOSE_BATCH_SIZE}"
echo "Embedding batch size: ${EMBED_BATCH_SIZE}"

python scripts/image_processing/build_object_dataset_and_embeddings.py \
  --input_stack_tif "${STACK_TIF}" \
  --annotation_results_csv "${ANNOTATION_CSV}" \
  --manifest_csv "${MANIFEST_CSV}" \
  --output_labels_tif "${OUT_LABELS_TIF}" \
  --output_feature_table "${OUT_FEATURES}" \
  --output_meta_json "${OUT_META}" \
  --cellpose_batch_size "${CELLPOSE_BATCH_SIZE}" \
  --embed_batch_size "${EMBED_BATCH_SIZE}"

echo "=========================================================="
echo "Job Finished at $(date)"
echo "=========================================================="
