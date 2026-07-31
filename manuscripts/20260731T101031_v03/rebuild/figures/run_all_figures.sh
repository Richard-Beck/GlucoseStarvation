#!/usr/bin/env bash
set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ASSEMBLY_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PROJECT_ROOT="$(cd "${ASSEMBLY_ROOT}/../.." && pwd)"

if (( $# > 1 )); then
  echo "usage: $0 [OUTPUT_DIR]" >&2
  exit 2
fi

OUTPUT_DIR="${1:-${SCRIPT_DIR}/rebuilt_figures}"
OUTPUT_DIR="$(mkdir -p "${OUTPUT_DIR}" && cd "${OUTPUT_DIR}" && pwd)"
if find "${OUTPUT_DIR}" -mindepth 1 -print -quit | grep -q .; then
  echo "output directory must be empty: ${OUTPUT_DIR}" >&2
  exit 2
fi

WORK_ROOT="$(mktemp -d "${SCRIPT_DIR}/.validation_work.XXXXXX")"
trap 'rm -rf -- "${WORK_ROOT}"' EXIT
WORK_PACKAGES="${WORK_ROOT}/packages"
mkdir -p "${WORK_PACKAGES}"
cp -a "${SCRIPT_DIR}/package_scripts/." "${WORK_PACKAGES}/"

cd "${PROJECT_ROOT}"

scripts/agentRrunner.sh \
  "${WORK_PACKAGES}/resegmentation_core/scripts/polish_figures.R" \
  --phase final --project-root "${PROJECT_ROOT}"

scripts/agentRrunner.sh \
  "${WORK_PACKAGES}/mechanistic_diagnostics/polish_figures.R" \
  --output-root "${WORK_PACKAGES}/mechanistic_diagnostics"

scripts/agentRrunner.sh \
  "${WORK_PACKAGES}/posterior_size_context/polish_figures.R" \
  --output-root "${WORK_PACKAGES}/posterior_size_context"

scripts/agentRrunner.sh \
  "${WORK_PACKAGES}/posterior_strategy_F5/scripts/polish_figures.R" \
  --phase final

scripts/agentRrunner.sh \
  "${WORK_PACKAGES}/posterior_strategy_context/scripts/polish_figures.R" \
  --phase final

scripts/agentRrunner.sh \
  "${WORK_PACKAGES}/morphology_metrics/scripts/polish_figures.R" \
  --phase final --project-root "${PROJECT_ROOT}"

scripts/agentRrunner.sh \
  "${WORK_PACKAGES}/sum159_label_swap/scripts/polish_figures.R" \
  --phase final

"${ASSEMBLY_ROOT}/rebuild/python_runner.sh" "${SCRIPT_DIR}/normalize_rebuild_outputs.py" \
  "${WORK_PACKAGES}" "${OUTPUT_DIR}"

echo "Rebuilt and byte-validated 23 figures in ${OUTPUT_DIR}"
