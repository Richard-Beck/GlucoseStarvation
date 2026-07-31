#!/usr/bin/env bash
set -euo pipefail

OPTIM_DIR=${1:?optimization directory required}

for ((attempt=1; attempt<=2880; attempt++)); do
  if [[ -f "${OPTIM_DIR}/optim_draws_all.Rds" && \
        -f "${OPTIM_DIR}/optim_lp_all.Rds" && \
        -f "${OPTIM_DIR}/optim_rc_all.Rds" ]]; then
    echo "optimization_combined=${OPTIM_DIR}"
    exit 0
  fi
  sleep 30
done

echo "Timed out waiting for combined optimization files: ${OPTIM_DIR}" >&2
exit 1
