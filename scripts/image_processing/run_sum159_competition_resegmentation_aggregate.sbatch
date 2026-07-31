#!/bin/bash

#SBATCH --job-name=sum159_comp_agg
#SBATCH --time=00:30:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

MANIFEST_TSV=$1
RUN_DIR=$2
PROJECT_ROOT=$3

cd "${PROJECT_ROOT}"

python3 "${RUN_DIR}/code/aggregate_sum159_competition_resegmentation.py" \
  --manifest_tsv "${MANIFEST_TSV}" \
  --run_dir "${RUN_DIR}"
