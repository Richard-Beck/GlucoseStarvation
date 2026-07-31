#!/usr/bin/env bash
set -euo pipefail

if [[ -n "${MANUSCRIPT_PYTHON:-}" ]]; then
  PYTHON_BIN="${MANUSCRIPT_PYTHON}"
elif command -v python3 >/dev/null 2>&1 && python3 -c 'import sys; raise SystemExit(sys.version_info < (3, 9))'; then
  PYTHON_BIN="$(command -v python3)"
elif [[ -x /home/4473331/.conda/envs/cpose/bin/python ]]; then
  PYTHON_BIN="/home/4473331/.conda/envs/cpose/bin/python"
else
  echo "A Python >=3.9 interpreter is required. Set MANUSCRIPT_PYTHON to its path." >&2
  exit 2
fi

exec "${PYTHON_BIN}" "$@"
