#!/bin/bash

#SBATCH --job-name=wp3_cpsam_nuclear
#SBATCH --output=logs/wp3-cpsam-nuclear-%j.out
#SBATCH --error=logs/wp3-cpsam-nuclear-%j.err
#SBATCH --time=12:00:00
#SBATCH --qos=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=96G
#SBATCH --gres=gpu:a30:1

set -euo pipefail

MANIFEST_CSV=${1:-data/report_exports/wp3_nuclear_size/pilot_manifest.csv}
OUTPUT_ROOT=${2:-data/image_processing_runs/wp3_nuclear_size_pilot}
RAW_DATA_DIR=${3:-all_raw}
MAX_IMAGES=${MAX_IMAGES:-1000}
QC_COUNT=${QC_COUNT:-12}
CELLPOSE_BATCH_SIZE=${CELLPOSE_BATCH_SIZE:-8}

echo "=========================================================="
echo "WP3 CPSAM nuclear-size pilot started at $(date)"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Node: ${SLURMD_NODENAME:-NA}"
echo "Manifest: ${MANIFEST_CSV}"
echo "Output root: ${OUTPUT_ROOT}"
echo "Raw data dir: ${RAW_DATA_DIR}"
echo "Max images: ${MAX_IMAGES}"
echo "=========================================================="

mkdir -p logs
nvidia-smi || true

set +u
source /app/eb/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh
conda activate cpose
set -u

python scripts/image_processing/wp3_nuclear_size_pilot.py \
  --manifest_csv "${MANIFEST_CSV}" \
  --raw_data_dir "${RAW_DATA_DIR}" \
  --output_root "${OUTPUT_ROOT}" \
  --max_images "${MAX_IMAGES}" \
  --cellpose_batch_size "${CELLPOSE_BATCH_SIZE}" \
  --qc_count "${QC_COUNT}" \
  --save_masks

echo "=========================================================="
echo "WP3 CPSAM nuclear-size pilot finished at $(date)"
echo "=========================================================="
