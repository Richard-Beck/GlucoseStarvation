#!/usr/bin/env bash
set -euo pipefail

PRELIM_COUNTS_CSV=${1:-data/image_processing_runs/run_20260324_233122/prelim_counts_by_image.csv}
RAW_DATA_DIR=${2:-all_raw}
OUTPUT_ROOT=${3:-data/image_processing_runs/wp3_nuclear_size_pilot}
MAX_IMAGES=${MAX_IMAGES:-90}
IMAGES_PER_STRATUM=${IMAGES_PER_STRATUM:-1}
GLUCOSE_BINS=${GLUCOSE_BINS:-low,mid,high}

RUN_ID=$(date +"%Y%m%d_%H%M%S")
SLURM_RUN_DIR="slurm/runs/wp3_nuclear_size/${RUN_ID}"
MANIFEST_CSV="data/report_exports/wp3_nuclear_size/pilot_manifest_${RUN_ID}.csv"
MANIFEST_META_JSON="data/report_exports/wp3_nuclear_size/pilot_manifest_${RUN_ID}.json"

mkdir -p "${SLURM_RUN_DIR}" logs data/report_exports/wp3_nuclear_size "${OUTPUT_ROOT}"

python scripts/image_processing/make_wp3_nuclear_pilot_manifest.py \
  --prelim_counts_csv "${PRELIM_COUNTS_CSV}" \
  --output_csv "${MANIFEST_CSV}" \
  --output_meta_json "${MANIFEST_META_JSON}" \
  --images_per_stratum "${IMAGES_PER_STRATUM}" \
  --max_images "${MAX_IMAGES}" \
  --glucose_bins "${GLUCOSE_BINS}"

cp "${MANIFEST_CSV}" "${SLURM_RUN_DIR}/pilot_manifest.csv"
cp "${MANIFEST_META_JSON}" "${SLURM_RUN_DIR}/pilot_manifest_meta.json"

export MAX_IMAGES
JOB_ID=$(sbatch --parsable \
  scripts/image_processing/run_wp3_nuclear_size_pilot.sbatch \
  "${MANIFEST_CSV}" \
  "${OUTPUT_ROOT}" \
  "${RAW_DATA_DIR}")

{
  echo "run_id=${RUN_ID}"
  echo "job_id=${JOB_ID}"
  echo "submitted_at=$(date --iso-8601=seconds)"
  echo "prelim_counts_csv=${PRELIM_COUNTS_CSV}"
  echo "raw_data_dir=${RAW_DATA_DIR}"
  echo "output_root=${OUTPUT_ROOT}"
  echo "manifest_csv=${MANIFEST_CSV}"
  echo "manifest_meta_json=${MANIFEST_META_JSON}"
  echo "max_images=${MAX_IMAGES}"
  echo "images_per_stratum=${IMAGES_PER_STRATUM}"
  echo "glucose_bins=${GLUCOSE_BINS}"
} > "${SLURM_RUN_DIR}/submission.env"

echo "${JOB_ID}" > "${SLURM_RUN_DIR}/job_ids.txt"
echo "Submitted WP3 CPSAM nuclear-size pilot job ${JOB_ID}"
echo "Manifest: ${MANIFEST_CSV}"
echo "Submission metadata: ${SLURM_RUN_DIR}"
