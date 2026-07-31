#!/usr/bin/env bash
set -euo pipefail

ARRAY_JOB_ID=${1:?array job id required}
N_STARTS=${2:?number of starts required}
ORIGINAL_TIME=${3:?original wall time required}
OUTPUT_DIR=${4:?output directory required}
LOG_DIR=${5:?log directory required}
MODEL_NAME=${6:?model name required}
MODEL_VERSION=${7:?model version required}
MODEL_ID=${8:?model id required}
THREADS=${9:?threads required}
STAN_DATA_PATH=${10:?stan data path required}
DATASET_ID=${11:?dataset id required}
CPUS=${12:?cpus required}
MEM_GB=${13:?memory required}
QOS=${14:?qos required}
RECORD_PATH=${15:?record path required}

mkdir -p "$LOG_DIR" "$(dirname "$RECORD_PATH")"

time_to_seconds() {
  local value=$1 days=0 clock=$1
  if [[ "$value" == *-* ]]; then
    days=${value%%-*}
    clock=${value#*-}
  fi
  IFS=: read -r hours minutes seconds <<< "$clock"
  seconds=${seconds:-0}
  echo $((10#$days * 86400 + 10#$hours * 3600 + 10#$minutes * 60 + 10#$seconds))
}

seconds_to_time() {
  local total=$1
  local days=$((total / 86400)) remainder=$((total % 86400))
  local hours=$((remainder / 3600)) minutes=$(((remainder % 3600) / 60)) seconds=$((remainder % 60))
  if ((days > 0)); then
    printf '%d-%02d:%02d:%02d' "$days" "$hours" "$minutes" "$seconds"
  else
    printf '%02d:%02d:%02d' "$hours" "$minutes" "$seconds"
  fi
}

declare -A states=()
for attempt in 1 2 3 4 5 6; do
  while IFS='|' read -r job_id state; do
    [[ "$job_id" =~ ^${ARRAY_JOB_ID}_([0-9]+)$ ]] || continue
    states["${BASH_REMATCH[1]}"]=$state
  # On this cluster JobIDRaw is the internal numeric record for array tasks;
  # JobID retains the submit-time <array>_<task> identity needed below.
  done < <(sacct -n -P -X -j "$ARRAY_JOB_ID" --format=JobID,State)
  ((${#states[@]} >= N_STARTS)) && break
  sleep 5
done

timeouts=()
other_failures=()
for ((task=1; task<=N_STARTS; task++)); do
  complete=1
  for prefix in optim_draws optim_lp optim_rc; do
    [[ -f "${OUTPUT_DIR}/${prefix}_${task}.Rds" ]] || complete=0
  done
  ((complete == 1)) && continue
  state=${states[$task]:-UNKNOWN}
  if [[ "$state" == TIMEOUT* ]]; then
    timeouts+=("$task")
  else
    other_failures+=("${task}:${state}")
  fi
done

if ((${#other_failures[@]})); then
  printf 'Non-timeout failures prevent automatic recovery: %s\n' "${other_failures[*]}" >&2
  printf 'initial_array_job_id=%s\nstatus=FAILED_NON_TIMEOUT\nfailures=%s\n' \
    "$ARRAY_JOB_ID" "${other_failures[*]}" > "$RECORD_PATH"
  exit 1
fi

retry_job_id=""
retry_time=""
if ((${#timeouts[@]})); then
  retry_time=$(seconds_to_time $((3 * $(time_to_seconds "$ORIGINAL_TIME"))))
  retry_spec=$(IFS=,; echo "${timeouts[*]}")
  retry_job_id=$(sbatch --parsable \
    --job-name="pe_opt_retry" --array="$retry_spec" \
    --cpus-per-task="$CPUS" --mem="${MEM_GB}G" --time="$retry_time" --qos="$QOS" \
    --output="${LOG_DIR}/%x_%A_%a.out" --error="${LOG_DIR}/%x_%A_%a.err" \
    slurm/jobs/gpath_optim_array.sh \
    "$MODEL_NAME" "$MODEL_VERSION" "$MODEL_ID" "$THREADS" \
    "$STAN_DATA_PATH" "$DATASET_ID" "$OUTPUT_DIR")
  combine_dependency="afterok:${retry_job_id}"
else
  combine_dependency=""
fi

combine_args=(
  --parsable --job-name="pe_opt_combine"
  --cpus-per-task=2 --mem=4G --time=00:30:00 --qos="$QOS"
  --output="${LOG_DIR}/%x_%j.out" --error="${LOG_DIR}/%x_%j.err"
)
[[ -n "$combine_dependency" ]] && combine_args+=(--dependency="$combine_dependency")
combine_job_id=$(sbatch "${combine_args[@]}" slurm/jobs/gpath_optim_combine.sh "$OUTPUT_DIR" 0)

{
  echo "initial_array_job_id=${ARRAY_JOB_ID}"
  echo "timeout_tasks=$(IFS=,; echo "${timeouts[*]-}")"
  echo "retry_job_id=${retry_job_id}"
  echo "retry_time=${retry_time}"
  echo "combine_job_id=${combine_job_id}"
  echo "status=RECOVERY_SUBMITTED"
} > "$RECORD_PATH"

echo "timeout_tasks=${#timeouts[@]} retry_job_id=${retry_job_id:-none} combine_job_id=${combine_job_id}"
