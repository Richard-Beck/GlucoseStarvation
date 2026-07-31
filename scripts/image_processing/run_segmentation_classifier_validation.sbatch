#!/bin/bash
#SBATCH --job-name=segclass_val90
#SBATCH --partition=red
#SBATCH --qos=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:a30:1
#SBATCH --time=00:30:00

set -euo pipefail

PROJECT_ROOT=${1:?project root required}
RELEASE_DIR=${2:?release directory required}
FROZEN_RUNNER=${3:?frozen run_shard.py required}

cd "${PROJECT_ROOT}"
set +u
source /app/eb/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh
conda activate cpose
set -u

nvidia-smi
python "${FROZEN_RUNNER}" \
  --manifest "${RELEASE_DIR}/pipeline_manifest.csv" \
  --output-dir "${RELEASE_DIR}/pipeline" \
  --classifier data/train/classifier_training_outputs/object_classifier.pkl \
  --expected-classifier-sha256 e20436c7a5ef141f06996b47982bcfaf5f17bb1a799f601d297008da6d4b3aed \
  --backbone /home/4473331/.cache/torch/hub/checkpoints/resnet18-f37072fd.pth \
  --expected-backbone-sha256 f37072fd47e89c5e827621c5baffa7500819f7896bbacec160b1a16c560e07ec \
  --cellpose-batch-size 8 \
  --embedding-batch-size 64 \
  --min-cell-area-px 50 \
  --min-nucleus-area-px 10 \
  --qc-count 8 \
  --save-masks

python scripts/image_processing/build_segmentation_classifier_validation.py summarize \
  --release-dir "${RELEASE_DIR}" \
  --frozen-runner "${FROZEN_RUNNER}" \
  --point-snap-radius 3
