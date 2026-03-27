#!/bin/bash

#----------------------------------------------------------------#
# SBATCH - SLURM Resource Directives
#----------------------------------------------------------------#
#SBATCH --job-name=cellpose_object_features
#SBATCH --output=logs/object-features-%j.out
#SBATCH --error=logs/object-features-%j.err

#SBATCH --time=16:00:00
#SBATCH --qos=large
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
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

OUTPUT_ROOT=${1:-data/image_processing_runs}
RAW_DATA_DIR=${2:-all_raw}

echo "Python environment activated. Running the script..."
echo "Output root: ${OUTPUT_ROOT}"
echo "Raw data dir: ${RAW_DATA_DIR}"

python scripts/image_processing/batch_extract_object_features.py \
  --raw_data_dir "${RAW_DATA_DIR}" \
  --output_root "${OUTPUT_ROOT}" \
  --max_workers "${SLURM_CPUS_PER_TASK:-6}"

echo "=========================================================="
echo "Job Finished at $(date)"
echo "=========================================================="
